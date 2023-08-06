#include <cfloat>
#include "bcMap.hpp"
#include "neknek.hpp"
#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "pointInterpolation.hpp"
#include <vector>
#include <algorithm>

namespace {

std::vector<std::string> neknekSolveFields(){
  auto fields = fieldsToSolve(platform->options);
  
  // mesh velocity is "derived" from fluid velocity, so it is not relevant here
  fields.erase(std::remove(fields.begin(), fields.end(), "mesh"), fields.end());
  return fields;
}

bool isIntBc(int bcType, std::string field){
  bool isInt = bcType == bcMap::bcTypeINT;
  
  if(field.find("scalar") != std::string::npos)
    isInt = bcType == bcMap::bcTypeINTS;
  
  return isInt;
}

void reserveAllocation(nrs_t *nrs, dlong npt)
{
  neknek_t *neknek = nrs->neknek;

  if (neknek->pointMap.size() && neknek->npt == npt)
    return;

  if (neknek->o_U.size()) {
    neknek->o_U.free();
  }

  if (neknek->o_S.size()) {
    neknek->o_S.free();
  }

  if (neknek->o_pointMap.size()) {
    neknek->o_pointMap.free();
  }

  // compute page-aligned fieldOffset
  neknek->fieldOffset = npt;
  const int pageW = ALIGN_SIZE / sizeof(dfloat);
  if (neknek->fieldOffset % pageW)
    neknek->fieldOffset = (neknek->fieldOffset / pageW + 1) * pageW;

  neknek->pointMap.resize(nrs->fieldOffset + 1);
  neknek->o_pointMap = platform->device.malloc((nrs->fieldOffset + 1) * sizeof(dlong));

  if (npt) {
    neknek->o_U =
        platform->device.malloc(nrs->NVfields * neknek->fieldOffset * (neknek->nEXT + 1) * sizeof(dfloat));
    if (neknek->Nscalar) {
      neknek->o_S = platform->device.malloc(neknek->Nscalar * neknek->fieldOffset * (neknek->nEXT + 1) *
                                            sizeof(dfloat));
    }
    else {
      neknek->o_S = platform->device.malloc(1 * sizeof(dfloat));
    }
  }
  else {
    neknek->o_U = platform->device.malloc(1 * sizeof(dfloat));
    neknek->o_S = platform->device.malloc(1 * sizeof(dfloat));
  }
  neknek->npt = npt;
}

void updateInterpPoints(nrs_t *nrs)
{
  // called in case of moving mesh ONLY
  if (!nrs->neknek->globalMovingMesh)
    return;

  auto *neknek = nrs->neknek;
  const dlong nsessions = neknek->nsessions;
  const dlong sessionID = neknek->sessionID;

  auto *mesh = nrs->meshV;

  // Setup findpts
  const dfloat tol = 5e-13;
  constexpr dlong npt_max = 128;
  const dfloat bb_tol = 0.01;

  auto &device = platform->device.occaDevice();

  // TODO: possible to cache this in moving mesh case?
  std::vector<std::shared_ptr<pointInterpolation_t>> sessionInterpolators(nsessions);
  for (dlong i = 0; i < nsessions; ++i) {
    sessionInterpolators[i] = std::make_shared<pointInterpolation_t>(nrs, bb_tol, tol, i == sessionID);
    sessionInterpolators[i]->setTimerLevel(TimerLevel::Basic);
    sessionInterpolators[i]->setTimerName("neknek_t::");
  }

  neknek->interpolator.reset();
  neknek->interpolator = std::make_shared<pointInterpolation_t>(nrs, bb_tol, tol);
  neknek->interpolator->setTimerLevel(TimerLevel::Basic);
  neknek->interpolator->setTimerName("neknek_t::");

  // neknekX[:] = mesh->x[pointMap[:]]
  neknek->copyNekNekPointsKernel(mesh->Nlocal,
                                 neknek->o_pointMap,
                                 mesh->o_x,
                                 mesh->o_y,
                                 mesh->o_z,
                                 neknek->o_x,
                                 neknek->o_y,
                                 neknek->o_z);

  // add points (use GPU version)
  for (dlong sess = 0; sess < nsessions; ++sess) {
    const auto nPoint = (sess == sessionID) ? 0 : neknek->npt;
    sessionInterpolators[sess]->setPoints(nPoint, neknek->o_x, neknek->o_y, neknek->o_z);
  }

  neknek->interpolator->setPoints(neknek->npt, neknek->o_x, neknek->o_y, neknek->o_z);

  const auto warningLevel = pointInterpolation_t::VerbosityLevel::Detailed;
  for (dlong sess = 0; sess < nsessions; ++sess) {
    sessionInterpolators[sess]->find(warningLevel);
  }

  auto &sessionData = neknek->interpolator->data();

  // TODO: possible to move to GPU?
  // copy results from other session into the point interpolator
  for (dlong sess = 0; sess < nsessions; ++sess) {
    auto data = sessionInterpolators[sess]->data();
    const auto nPoint = (sess == sessionID) ? 0 : neknek->npt;
    for (dlong pt = 0; pt < nPoint; ++pt) {
      sessionData.code[pt] = data.code[pt];
      sessionData.proc[pt] = data.proc[pt];
      sessionData.el[pt] = data.el[pt];

      sessionData.r[3 * pt + 0] = data.r[3 * pt + 0];
      sessionData.r[3 * pt + 1] = data.r[3 * pt + 1];
      sessionData.r[3 * pt + 2] = data.r[3 * pt + 2];

      sessionData.dist2[pt] = data.dist2[pt];
    }
  }

  neknek->interpolator->update();
}

dlong computeNumInterpPoints(nrs_t *nrs)
{
  auto *mesh = nrs->meshV;
  auto fields = neknekSolveFields();

  if(fields.size() == 0) return 0;

  // number of interpolated faces is constant across all solved fields
  dlong numInterpFaces = 0;
  for (dlong e = 0; e < mesh->Nelements; ++e) {
    for (dlong f = 0; f < mesh->Nfaces; ++f) {
      auto bID = mesh->EToB[f + mesh->Nfaces * e];
      auto bcType = bcMap::id(bID, fields[0]);
      numInterpFaces += (isIntBc(bcType, fields[0]));
    }
  }
  return numInterpFaces * mesh->Nfp;
}

void findInterpPoints(nrs_t *nrs)
{
  auto *neknek = nrs->neknek;
  const dlong nsessions = neknek->nsessions;
  const dlong sessionID = neknek->sessionID;

  auto *mesh = nrs->meshV;

  // Setup findpts
  const dfloat tol = 5e-13;
  constexpr dlong npt_max = 128;
  const dfloat bb_tol = 0.01;

  auto &device = platform->device.occaDevice();

  // TODO: possible to cache this in moving mesh case?
  std::vector<std::shared_ptr<pointInterpolation_t>> sessionInterpolators(nsessions);
  for (dlong i = 0; i < nsessions; ++i) {
    sessionInterpolators[i] = std::make_shared<pointInterpolation_t>(nrs, bb_tol, tol, i == sessionID);
    sessionInterpolators[i]->setTimerLevel(TimerLevel::Basic);
    sessionInterpolators[i]->setTimerName("neknek_t::");
  }

  neknek->interpolator.reset();
  neknek->interpolator = std::make_shared<pointInterpolation_t>(nrs, bb_tol, tol);
  neknek->interpolator->setTimerLevel(TimerLevel::Basic);
  neknek->interpolator->setTimerName("neknek_t::");

  auto numPoints = computeNumInterpPoints(nrs);
  reserveAllocation(nrs, numPoints);

  std::vector<dfloat> neknekX(numPoints, 0.0);
  std::vector<dfloat> neknekY(numPoints, 0.0);
  std::vector<dfloat> neknekZ(numPoints, 0.0);

  auto fields = neknekSolveFields();

  std::fill(neknek->pointMap.begin(), neknek->pointMap.end(), -1);

  if(fields.size()){
    dlong ip = 0;
    for (dlong e = 0; e < mesh->Nelements; ++e) {
      for (dlong f = 0; f < mesh->Nfaces; ++f) {

        for (dlong m = 0; m < mesh->Nfp; ++m) {
          dlong id = mesh->Nfaces * mesh->Nfp * e + mesh->Nfp * f + m;
          dlong idM = mesh->vmapM[id];

          auto bID = mesh->EToB[f + mesh->Nfaces * e];
          auto bcType = bcMap::id(bID, fields[0]);

          if (isIntBc(bcType, fields[0])) {
            neknekX[ip] = mesh->x[idM];
            neknekY[ip] = mesh->y[idM];
            neknekZ[ip] = mesh->z[idM];

            neknek->pointMap[idM] = ip;
            ++ip;
          }
        }
      }
    }
  }

  neknek->pointMap[nrs->fieldOffset] = neknek->fieldOffset;
  neknek->o_pointMap.copyFrom(neknek->pointMap.data());

  // add points
  for (dlong sess = 0; sess < nsessions; ++sess) {
    const auto nPoint = (sess == sessionID) ? 0 : numPoints;
    sessionInterpolators[sess]->setPoints(nPoint, neknekX.data(), neknekY.data(), neknekZ.data());
  }

  neknek->interpolator->setPoints(numPoints, neknekX.data(), neknekY.data(), neknekZ.data());

  const auto warningLevel = pointInterpolation_t::VerbosityLevel::Detailed;
  for (dlong sess = 0; sess < nsessions; ++sess) {
    sessionInterpolators[sess]->find(warningLevel);
  }

  auto &sessionData = neknek->interpolator->data();

  // copy results from other session into the point interpolator
  for (dlong sess = 0; sess < nsessions; ++sess) {
    auto data = sessionInterpolators[sess]->data();
    const auto nPoint = (sess == sessionID) ? 0 : numPoints;
    for (dlong pt = 0; pt < nPoint; ++pt) {
      sessionData.code[pt] = data.code[pt];
      sessionData.proc[pt] = data.proc[pt];
      sessionData.el[pt] = data.el[pt];

      sessionData.r[3 * pt + 0] = data.r[3 * pt + 0];
      sessionData.r[3 * pt + 1] = data.r[3 * pt + 1];
      sessionData.r[3 * pt + 2] = data.r[3 * pt + 2];

      sessionData.dist2[pt] = data.dist2[pt];
    }
  }

  neknek->interpolator->update();

  // allocate device coordinates for later use
  if (neknek->globalMovingMesh) {
    neknek->o_x = platform->device.malloc(neknek->npt * sizeof(dfloat), neknekX.data());
    neknek->o_y = platform->device.malloc(neknek->npt * sizeof(dfloat), neknekY.data());
    neknek->o_z = platform->device.malloc(neknek->npt * sizeof(dfloat), neknekZ.data());
  }
}

void neknekSetup(nrs_t *nrs)
{
  neknek_t *neknek = nrs->neknek;

  nrsCheck(platform->options.compareArgs("CONSTANT FLOW RATE", "TRUE"), platform->comm.mpiComm, EXIT_FAILURE,
           "%s\n", "constant flow rate support not supported");

  const dlong nsessions = neknek->nsessions;

  dlong globalRank;
  MPI_Comm_rank(platform->comm.mpiCommParent, &globalRank);

  if(globalRank == 0) printf("configuring neknek with %d sessions\n", nsessions);

  std::vector<int> Nscalars(nsessions, -1);
  Nscalars[neknek->sessionID] = nrs->Nscalar;

  MPI_Allreduce(MPI_IN_PLACE, Nscalars.data(), nsessions, MPI_DLONG, MPI_MAX, platform->comm.mpiCommParent);
  auto minNscalar = *std::min_element(Nscalars.begin(),Nscalars.end());

  bool allSame = std::all_of(Nscalars.begin(), Nscalars.end(), [minNscalar] (auto v) {return v == minNscalar;});

  if(globalRank == 0 && !allSame){
    std::cout << "WARNING: Nscalar is not the same across all sessions. Using the minimum value: " << minNscalar << "\n";
  }

  neknek->Nscalar = minNscalar;
  
  const dlong movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
  dlong globalMovingMesh;
  MPI_Allreduce(&movingMesh, &globalMovingMesh, 1, MPI_DLONG, MPI_MAX, platform->comm.mpiCommParent);
  neknek->globalMovingMesh = globalMovingMesh;

  findInterpPoints(nrs);
}

} // namespace

bool neknekCoupled()
{
  auto solverFields = neknekSolveFields();

  std::set<int> intBIDs;
  std::ostringstream errorLogger;

  int coupled = 0;
  for(auto&& field : solverFields){
    for(int bID = 1; bID <= bcMap::size(field); ++bID)
    {
      auto bcType = bcMap::id(bID, field);
      bool isInt = isIntBc(bcType, field);

      if(isInt){
        intBIDs.insert(bID);
        coupled += 1;
      }

      // INT bid in one field -> INT bid in all fields
      if((intBIDs.find(bID) != intBIDs.end()) && !isInt){
        errorLogger << "ERROR: expected INT boundary condition on boundary id "
                    << bID << " for field " << field << "\n";
      }
    }
  }

  int errorLength = errorLogger.str().length();
  MPI_Allreduce(MPI_IN_PLACE, &errorLength, 1, MPI_INT, MPI_MAX, platform->comm.mpiCommParent);

  nrsCheck(errorLength > 0, platform->comm.mpiCommParent, EXIT_FAILURE,
           "%s\n", errorLogger.str().c_str());

  MPI_Allreduce(MPI_IN_PLACE, &coupled, 1, MPI_INT, MPI_MAX, platform->comm.mpiCommParent);
  return coupled > 0;
}

neknek_t::neknek_t(nrs_t *nrs, dlong _nsessions, dlong _sessionID)
    : nsessions(_nsessions), sessionID(_sessionID)
{

  nrs->neknek = this;
  if (nrs->cds) {
    nrs->cds->neknek = this;
  }

  this->nEXT = 1;
  if(!platform->options.getArgs("NEKNEK BOUNDARY EXT ORDER").empty())
    platform->options.getArgs("NEKNEK BOUNDARY EXT ORDER", this->nEXT);

  // set boundary ext order to report to user, if not specified
  platform->options.setArgs("NEKNEK BOUNDARY EXT ORDER", std::to_string(this->nEXT));

  this->coeffEXT.resize(this->nEXT);
  this->o_coeffEXT = platform->device.malloc(this->nEXT * sizeof(dfloat));

  neknekSetup(nrs);

  // variable p0th + nek-nek is not supported
  int issueError = 0;
  if(nrs->pSolver){
    if (nrs->pSolver->allNeumann && platform->options.compareArgs("LOWMACH", "TRUE")) {
      issueError = 1;
    }
  }

  nrsCheck(issueError, platform->comm.mpiCommParent, EXIT_FAILURE, "%s\n", "variable p0th is not supported!");

  this->copyNekNekPointsKernel = platform->kernels.get("copyNekNekPoints");
}

void neknek_t::updateBoundary(nrs_t *nrs, int tstep, int stage)
{
  // do not invoke barrier -- this is performed later
  platform->timer.tic("neknek update boundary", 0);

  // do not invoke barrier in timer_t::tic
  platform->timer.tic("neknek sync", 0);
  MPI_Barrier(platform->comm.mpiCommParent);
  platform->timer.toc("neknek sync");
  this->tSync = platform->timer.query("neknek sync", "HOST:MAX");

  if (this->globalMovingMesh) {
    platform->timer.tic("neknek updateInterpPoints", 1);
    updateInterpPoints(nrs);
    platform->timer.toc("neknek updateInterpPoints");
  }

  platform->timer.tic("neknek exchange", 1);

  this->interpolator->eval(nrs->NVfields, nrs->fieldOffset, nrs->o_U, this->fieldOffset, this->o_U);

  if (this->Nscalar) {
    this->interpolator->eval(this->Nscalar,
      nrs->fieldOffset,
      nrs->cds->o_S,
      this->fieldOffset,
      this->o_S);
  }

  // lag state, update timestepper coefficients and compute extrapolated state
  if (stage == 1) {
    auto *mesh = nrs->meshV;
    int extOrder = std::min(tstep, this->nEXT);
    int bdfOrder = std::min(tstep, nrs->nBDF);
    nek::extCoeff(this->coeffEXT.data(), nrs->dt, extOrder, bdfOrder);

    for (int i = this->nEXT; i > extOrder; i--)
      this->coeffEXT[i - 1] = 0.0;

    this->o_coeffEXT.copyFrom(this->coeffEXT.data(), this->nEXT * sizeof(dfloat));

    for (int s = this->nEXT + 1; s > 1; s--) {
      auto Nbyte = nrs->NVfields * this->fieldOffset * sizeof(dfloat);
      this->o_U.copyFrom(this->o_U, Nbyte, (s - 1) * Nbyte, (s - 2) * Nbyte);

      Nbyte = this->Nscalar * this->fieldOffset * sizeof(dfloat);
      this->o_S.copyFrom(this->o_S, Nbyte, (s - 1) * Nbyte, (s - 2) * Nbyte);
    }

    auto o_Uold = this->o_U + this->fieldOffset * nrs->NVfields * sizeof(dfloat);
    auto o_Sold = this->o_S + this->fieldOffset * this->Nscalar * sizeof(dfloat);

    if(this->npt){
      nrs->extrapolateKernel(this->npt,
                             nrs->NVfields,
                             this->nEXT,
                             this->fieldOffset,
                             this->o_coeffEXT,
                             o_Uold,
                             this->o_U);
    }

    if (this->Nscalar && this->npt) {
      nrs->extrapolateKernel(this->npt,
                             this->Nscalar,
                             this->nEXT,
                             this->fieldOffset,
                             this->o_coeffEXT,
                             o_Sold,
                             this->o_S);
    }
  }

  platform->timer.toc("neknek exchange");

  this->tExch = platform->timer.query("neknek exchange", "DEVICE:MAX");
  this->ratio = this->tSync / this->tExch;

  platform->timer.toc("neknek update boundary");
}

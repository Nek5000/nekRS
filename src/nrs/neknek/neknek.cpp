#include <cfloat>
#include "platform.hpp"
#include "bcMap.hpp"
#include "neknek.hpp"
#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "pointInterpolation.hpp"

#include "sha1.hpp"

namespace
{

bool isIntBc(int bcType, std::string field)
{
  bool isInt = bcType == bcMap::bcTypeINT;

  if (field.find("scalar") != std::string::npos) {
    isInt = bcType == bcMap::bcTypeINTS;
  }

  return isInt;
}

} // namespace

bool neknekCoupled()
{
  int intFound = 0;
  for (auto &&field : nrsFieldsToSolve(platform->options)) {
    for (int bID = 1; bID <= bcMap::size(field); ++bID) {
      auto bcType = bcMap::id(bID, field);
      bool isInt = isIntBc(bcType, field);

      if (isInt) {
        intFound = 1;
      }
    }
  }

  // findpts functions have to be called collectively across all sessions
  MPI_Allreduce(MPI_IN_PLACE, &intFound, 1, MPI_INT, MPI_MAX, platform->comm.mpiCommParent);
  return intFound > 0;
}

void neknek_t::reserveAllocation()
{
  auto mesh = (nrs->cht) ? nrs->cds->mesh[0] : nrs->mesh;
  this->fieldOffset_ = alignStride<dfloat>(this->npt_);
  this->o_pointMap_ = platform->device.malloc<dlong>(mesh->Nlocal);

  int nStates = this->nEXT_ + 1;
  if (this->multirate()) {
    nStates += 1;
  }

  if (this->npt_) {
    if (std::find(this->fields.begin(), this->fields.end(), "velocity") != this->fields.end()) {
      this->o_U_ = platform->device.malloc<dfloat>(nrs->NVfields * this->fieldOffset_ * nStates);
    }
    if (this->Nscalar_) {
      this->o_S_ = platform->device.malloc<dfloat>(this->Nscalar_ * this->fieldOffset_ * nStates);
    }
    if (this->multirate()) {
      this->o_time_ = platform->device.malloc<dfloat>(this->fieldOffset_ * (maxOrd + 1));
    }
  }
}

void neknek_t::updateInterpPoints()
{
  // called in case of moving mesh ONLY
  if (!this->globalMovingMesh) {
    return;
  }

  auto mesh = (nrs->cht) ? nrs->cds->mesh[0] : nrs->mesh;

  this->interpolator.reset();
  this->interpolator =
      std::make_shared<pointInterpolation_t>(mesh, platform->comm.mpiCommParent, true, intBIDs);
  this->interpolator->setTimerName("neknek_t::");

  // neknekX[:] = mesh->x[pointMap[:]]
  this->copyNekNekPointsKernel(mesh->Nlocal,
                               this->o_pointMap_,
                               mesh->o_x,
                               mesh->o_y,
                               mesh->o_z,
                               this->o_x_,
                               this->o_y_,
                               this->o_z_);

  this->interpolator->setPoints(this->o_x_, this->o_y_, this->o_z_, this->o_session_);

  const auto verboseLevel = pointInterpolation_t::VerbosityLevel::Detailed;
  this->interpolator->find(verboseLevel);
}

void neknek_t::findIntPoints()
{
  const dlong sessionID = this->sessionID_;

  auto mesh = (nrs->cht) ? nrs->cds->mesh[0] : nrs->mesh;

  this->interpolator.reset();
  this->interpolator =
      std::make_shared<pointInterpolation_t>(mesh, platform->comm.mpiCommParent, true, intBIDs);
  this->interpolator->setTimerName("neknek_t::");

  // int points are the same for all neknek fields
  dlong numInterpFaces = 0;
  for (dlong e = 0; e < mesh->Nelements; ++e) {
    for (dlong f = 0; f < mesh->Nfaces; ++f) {
      auto bID = mesh->EToB[f + mesh->Nfaces * e];
      auto bcType = bcMap::id(bID, this->fields[0]);
      numInterpFaces += (isIntBc(bcType, this->fields[0]));
    }
  }
  const auto numPoints = numInterpFaces * mesh->Nfp;

  this->npt_ = numPoints;
  this->reserveAllocation();

  std::vector<dfloat> neknekX(numPoints, 0.0);
  std::vector<dfloat> neknekY(numPoints, 0.0);
  std::vector<dfloat> neknekZ(numPoints, 0.0);
  std::vector<dlong> session(numPoints, -1);

  std::vector<dlong> pointMap(mesh->Nlocal, -1);

  if (this->fields.size()) {
    auto [x, y, z] = mesh->xyzHost();

    dlong ip = 0;
    for (dlong e = 0; e < mesh->Nelements; ++e) {
      for (dlong f = 0; f < mesh->Nfaces; ++f) {

        for (dlong m = 0; m < mesh->Nfp; ++m) {
          dlong id = mesh->Nfaces * mesh->Nfp * e + mesh->Nfp * f + m;
          dlong idM = mesh->vmapM[id];

          auto bID = mesh->EToB[f + mesh->Nfaces * e];
          auto bcType = bcMap::id(bID, this->fields[0]);

          if (isIntBc(bcType, this->fields[0])) {
            neknekX[ip] = x[idM];
            neknekY[ip] = y[idM];
            neknekZ[ip] = z[idM];
            session[ip] = sessionID;
            pointMap[idM] = ip;
            ++ip;
          }
        }
      }
    }
  }

  this->o_pointMap_.copyFrom(pointMap.data());

  this->interpolator->setPoints(neknekX, neknekY, neknekZ, session);

  const auto verboseLevel = pointInterpolation_t::VerbosityLevel::Detailed;
  this->interpolator->find(verboseLevel);

  this->o_x_ = platform->device.malloc<dfloat>(this->npt_);
  this->o_x_.copyFrom(neknekX.data());
  this->o_y_ = platform->device.malloc<dfloat>(this->npt_);
  this->o_y_.copyFrom(neknekY.data());
  this->o_z_ = platform->device.malloc<dfloat>(this->npt_);
  this->o_z_.copyFrom(neknekZ.data());
  this->o_session_ = platform->device.malloc<dlong>(this->npt_);
  this->o_session_.copyFrom(session.data());
}

void neknek_t::setup()
{
  dlong globalRank;
  MPI_Comm_rank(platform->comm.mpiCommParent, &globalRank);

  const int nsessions = this->nsessions_;
  if (platform->comm.mpiRank == 0) {
    printf("initializing neknek with %d sessions\n", nsessions);
    std::fflush(stdout);
  }

  nekrsCheck(static_cast<bool>(platform->options.compareArgs("CONSTANT FLOW RATE", "TRUE")),
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "constant flow rate support not supported");

  if (nrs->pSolver) {
    nekrsCheck(
        static_cast<bool>(nrs->pSolver->nullSpace() && platform->options.compareArgs("LOWMACH", "TRUE")),
        platform->comm.mpiComm,
        EXIT_FAILURE,
        "%s\n",
        "variable p0th is not supported!");
  }

  int movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
  MPI_Allreduce(MPI_IN_PLACE, &movingMesh, 1, MPI_INT, MPI_MAX, platform->comm.mpiCommParent);
  this->globalMovingMesh = movingMesh;

  this->fields = [&]() {
    std::vector<std::string> list;
    for (auto &&field : nrsFieldsToSolve(platform->options)) {
      auto mesh = (field == "scalar00") ?  nrs->cds->mesh[0] : nrs->mesh;

      int intFound = 0;
      for (dlong e = 0; e < mesh->Nelements; ++e) {
        for (dlong f = 0; f < mesh->Nfaces; ++f) {
          if (isIntBc(bcMap::id(mesh->EToB[f + mesh->Nfaces * e], field), field)) {
            intFound = 1;
          }
        }
      }
      MPI_Allreduce(MPI_IN_PLACE, &intFound, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);
      if (intFound) {
        list.push_back(field);
      }
    }

    return list;
  }();

  // check if all exchanged fields (within the same session) share the same INT boundaries
  std::ostringstream errorLogger;
  std::set<int> intBIDFields;
  for (auto &&field : this->fields) {
    for (int bID = 1; bID <= bcMap::size(field); ++bID) {
      const auto isInt = isIntBc(bcMap::id(bID, field), field);

      if (isInt) {
        intBIDFields.insert(bID);
      }

      if ((intBIDFields.find(bID) != intBIDFields.end()) && !isInt) {
        errorLogger << "ERROR: expected INT boundary condition on boundary id " << bID << " for field "
                    << field << "\n";
      }
    }
  }
  int errorLength = errorLogger.str().length();
  MPI_Allreduce(MPI_IN_PLACE, &errorLength, 1, MPI_INT, MPI_MAX, platform->comm.mpiCommParent);
  nekrsCheck(errorLength > 0, platform->comm.mpiCommParent, EXIT_FAILURE, "%s\n", errorLogger.str().c_str());

  std::vector<int> scalarIndices(nrs->Nscalar, -1);
  this->Nscalar_ = 0;
  for (auto &&field : this->fields) {
    if (field.find("scalar") != std::string::npos) {
      const auto id = std::stoi(field.substr(std::string("scalar").length()));
      scalarIndices[id] = this->Nscalar_;
      this->Nscalar_++;
    }
  }

  this->o_scalarIndices_ = platform->device.malloc<int>(nrs->Nscalar, scalarIndices.data());

  for (int bID = 1; bID <= bcMap::size(this->fields[0]); ++bID) {
    if (isIntBc(bcMap::id(bID, this->fields[0]), this->fields[0])) {
      intBIDs.push_back(bID);
    }
  }

  this->findIntPoints();

  if (platform->comm.mpiRank == 0) {
    std::cout << "exchanged fields: ";
    for (auto &&field : this->fields) {
      std::cout << field << "  ";
    }
    std::cout << "\n";
  }
  std::fflush(stdout);

  // check if fields across all sessions match
  {
    std::string s;
    for (auto &&field : this->fields) {
      s += field;
    }

    SHA1 sha;
    sha.update(s);
    const auto hash = sha.final();
    const auto hashTruncated = hash.substr(hash.length() - 8);
    unsigned long intHash = std::stoul("0x" + hashTruncated, nullptr, 0);

    unsigned long intHashMin;
    unsigned long intHashMax;

    MPI_Allreduce(&intHash, &intHashMin, 1, MPI_UNSIGNED_LONG, MPI_MIN, platform->comm.mpiCommParent);
    MPI_Allreduce(&intHash, &intHashMax, 1, MPI_UNSIGNED_LONG, MPI_MAX, platform->comm.mpiCommParent);

    nekrsCheck(intHashMin != intHashMax,
               platform->comm.mpiCommParent,
               EXIT_FAILURE,
               "%s\n",
               "neknek fields do not match across all sessions");
  }

  if (platform->comm.mpiRank == 0) {
    std::cout << "done\n";
  }
}

neknek_t::neknek_t(nrs_t *nrs_, dlong nsessions, dlong sessionID)
    : nsessions_(nsessions), sessionID_(sessionID), nrs(nrs_)
{
  nrs->neknek = this;

  this->nEXT_ = 1;
  if (!platform->options.getArgs("NEKNEK BOUNDARY EXT ORDER").empty()) {
    platform->options.getArgs("NEKNEK BOUNDARY EXT ORDER", this->nEXT_);
  }

  // set boundary ext order to report to user, if not specified
  platform->options.setArgs("NEKNEK BOUNDARY EXT ORDER", std::to_string(this->nEXT_));

  this->multirate_ = platform->options.compareArgs("NEKNEK MULTIRATE TIMESTEPPER", "TRUE");

  this->coeffEXT.resize(this->nEXT_);
  this->o_coeffEXT = platform->device.malloc<dfloat>(this->nEXT_);

  this->setup();

  this->copyNekNekPointsKernel = platform->kernelRequests.load("copyNekNekPoints");
  this->computeFluxKernel = platform->kernelRequests.load("computeFlux");
  this->fixSurfaceFluxKernel = platform->kernelRequests.load("fixSurfaceFlux");
  this->extrapolateBoundaryKernel = platform->kernelRequests.load("extrapolateBoundary");
  this->mapScalarKernel = platform->kernelRequests.load("mapScalar");
}

void neknek_t::updateBoundary(int tstep, int stage, double time)
{
  if (multirate()) {
    extrapolateBoundary(tstep, time, predictorStep);
    return;
  }

  // do not invoke barrier -- this is performed later
  platform->timer.tic("neknek update boundary");

  const bool exchangeAllTimes = false;
  const bool lagState = (stage == 1);
  exchange(exchangeAllTimes, lagState);

  // lag state, update timestepper coefficients and compute extrapolated state
  if (stage == 1) {
    extrapolate(tstep);
  }

  platform->timer.toc("neknek update boundary");
}

occa::memory neknek_t::partitionOfUnity()
{
  if (!this->o_partition_.isInitialized()) {
    this->o_partition_ = platform->device.malloc<dfloat>(nrs->fieldOffset);
  }

  if (!recomputePartition) {
    return this->o_partition_;
  }
  recomputePartition = false;

  auto mesh = (nrs->cht) ? nrs->cds->mesh[0] : nrs->mesh;

  auto pointInterp = pointInterpolation_t(mesh, platform->comm.mpiCommParent, true, intBIDs);

  auto o_dist = pointInterp.distanceINT();

  auto o_sess = platform->deviceMemoryPool.reserve<dlong>(nrs->fieldOffset);
  auto o_sumDist = platform->deviceMemoryPool.reserve<dfloat>(nrs->fieldOffset);
  auto o_found = platform->deviceMemoryPool.reserve<dfloat>(nrs->fieldOffset);
  auto o_interpDist = platform->deviceMemoryPool.reserve<dfloat>(nrs->fieldOffset);
  o_sumDist.copyFrom(o_dist, mesh->Nlocal);

  std::vector<dfloat> found(mesh->Nlocal);
  std::vector<dlong> sessions(mesh->Nlocal);

  for (int sess = 0; sess < this->nsessions_; ++sess) {
    auto id = (sess + this->sessionID_) % this->nsessions_;
    if (id == this->sessionID_) {
      continue;
    }
    std::fill(sessions.begin(), sessions.end(), id);
    o_sess.copyFrom(sessions.data(), mesh->Nlocal);

    pointInterp.setPoints(mesh->o_x, mesh->o_y, mesh->o_z, o_sess);
    pointInterp.find(pointInterpolation_t::VerbosityLevel::None, true);

    auto &data = pointInterp.data();
    for (int n = 0; n < mesh->Nlocal; ++n) {
      found[n] = (data.code[n] == pointInterpolation_t::CODE_NOT_FOUND) ? 0.0 : 1.0;
    }

    o_found.copyFrom(found.data());
    pointInterp.eval(1, nrs->fieldOffset, o_dist, nrs->fieldOffset, o_interpDist);

    platform->linAlg->axmy(mesh->Nlocal, 1.0, o_found, o_interpDist);
    platform->linAlg->axpby(mesh->Nlocal, 1.0, o_interpDist, 1.0, o_sumDist);
  }

  // \Xi(x) = \dfrac{\delta^s(x)}{\sum_{i=1}^S \delta^s(x_i)}
  this->o_partition_.copyFrom(o_dist, mesh->Nlocal);
  platform->linAlg->aydx(mesh->Nlocal, 1.0, o_sumDist, this->o_partition_);

  o_sess.free();
  o_sumDist.free();
  o_found.free();
  o_interpDist.free();

  return this->o_partition_;
}

void neknek_t::lag()
{
  int nStates = this->nEXT_ + 1;
  if (this->multirate()) {
    nStates += 1;
  }

  for (int s = nStates; s > 1; s--) {
    auto N = nrs->NVfields * this->fieldOffset_;
    if (std::find(this->fields.begin(), this->fields.end(), "velocity") != this->fields.end()) {
      this->o_U_.copyFrom(this->o_U_, N, (s - 1) * N, (s - 2) * N);
    }

    N = this->Nscalar_ * this->fieldOffset_;
    this->o_S_.copyFrom(this->o_S_, N, (s - 1) * N, (s - 2) * N);
  }
}

void neknek_t::extrapolate(int tstep)
{
  int extOrder = std::min(tstep, this->nEXT_);
  int bdfOrder = std::min(tstep, nrs->nBDF);
  nek::extCoeff(this->coeffEXT.data(), nrs->dt, extOrder, bdfOrder);

  for (int i = this->nEXT_; i > extOrder; i--) {
    this->coeffEXT[i - 1] = 0.0;
  }

  this->o_coeffEXT.copyFrom(this->coeffEXT.data(), this->nEXT_);

  if (this->npt_) {
    if (std::find(this->fields.begin(), this->fields.end(), "velocity") != this->fields.end()) {
      auto o_Uold = this->o_U_ + this->fieldOffset_ * nrs->NVfields;
      nrs->extrapolateKernel(this->npt_,
                             nrs->NVfields,
                             this->nEXT_,
                             this->fieldOffset_,
                             this->o_coeffEXT,
                             o_Uold,
                             this->o_U_);
    }
  }

  if (this->Nscalar_ && this->npt_) {
    auto o_Sold = this->o_S_ + this->fieldOffset_ * this->Nscalar_;
    nrs->extrapolateKernel(this->npt_,
                           this->Nscalar_,
                           this->nEXT_,
                           this->fieldOffset_,
                           this->o_coeffEXT,
                           o_Sold,
                           this->o_S_);
  }
}

void neknek_t::exchange(bool allTimeStates, bool lagState)
{
  // do not invoke barrier in timer_t::tic
  platform->timer.tic("neknek sync");
  MPI_Barrier(platform->comm.mpiCommParent);
  platform->timer.toc("neknek sync");
  this->tSync_ = platform->timer.query("neknek sync", "HOST:MAX");

  if (this->globalMovingMesh) {
    platform->timer.tic("neknek updateInterpPoints");
    this->updateInterpPoints();
    platform->timer.toc("neknek updateInterpPoints");

    this->recomputePartition = true;
  }

  if (allTimeStates) {
    auto nrsOrder = std::max(nrs->nBDF, nrs->nEXT);
    nekrsCheck(nrsOrder < this->nEXT_,
               platform->comm.mpiComm,
               EXIT_FAILURE,
               "neknek extrapolation order (%d) exceeds nekRS order (%d)\n",
               this->nEXT_,
               nrsOrder);
  }

  const auto nStates = allTimeStates ? nEXT_ : 1;

  platform->timer.tic("neknek exchange");

  if (std::find(this->fields.begin(), this->fields.end(), "velocity") != this->fields.end()) {
    this->interpolator->eval(nStates * nrs->NVfields,
                             nrs->fieldOffset,
                             nrs->o_U,
                             this->fieldOffset_,
                             this->o_U_);
  }

  if (this->Nscalar_) {
    auto o_S = nrs->cds->o_S;
    if (this->Nscalar_ != nrs->Nscalar) {
      o_S = platform->deviceMemoryPool.reserve<dfloat>(nStates * this->Nscalar_ * nrs->fieldOffset);
      this->mapScalarKernel(nrs->cds->mesh[0]->Nlocal,
                            nrs->Nscalar,
                            nrs->fieldOffset,
                            nStates,
                            this->Nscalar_,
                            this->o_scalarIndices_,
                            nrs->cds->o_S,
                            o_S);
    }
    this->interpolator->eval(nStates * this->Nscalar_, nrs->fieldOffset, o_S, this->fieldOffset_, this->o_S_);
  }

  platform->timer.toc("neknek exchange");

  this->tExch_ = platform->timer.query("neknek exchange", "DEVICE:MAX");
  this->ratio_ = this->tSync_ / this->tExch_;

  if (lagState) {
    lag();
  }
}

double neknek_t::adjustDt(double dt)
{
  if (!this->multirate()) {
    double minDt = dt;
    MPI_Allreduce(MPI_IN_PLACE, &minDt, 1, MPI_DOUBLE, MPI_MIN, platform->comm.mpiCommParent);
    double maxDt = dt;
    MPI_Allreduce(MPI_IN_PLACE, &maxDt, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiCommParent);

    const auto relErr = std::abs(maxDt - minDt) / maxDt;
    nekrsCheck(relErr > 100 * std::numeric_limits<double>::epsilon(),
               platform->comm.mpiComm,
               EXIT_FAILURE,
               "Time step size needs to be the same across all sessions.\n"
               "Max dt = %e, Min dt = %e\n",
               maxDt,
               minDt);

    return dt;
  }

  double maxDt = dt;
  MPI_Allreduce(MPI_IN_PLACE, &maxDt, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiCommParent);

  double ratio = maxDt / dt;
  int timeStepRatio = std::floor(ratio);
  double maxErr = std::abs(ratio - timeStepRatio);

  MPI_Allreduce(MPI_IN_PLACE, &maxErr, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiCommParent);
  nekrsCheck(maxErr > 1e-4,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "Multirate time stepping requires a fixed integer time step size ratio\n"
             "Max dt = %e, dt = %e, ratio = %e, ratioErr = %e\n",
             maxDt,
             dt,
             ratio,
             maxErr);

  // rescale dt to be an _exact_ integer multiple of minDt
  dt = maxDt / timeStepRatio;
  platform->options.setArgs("NEKNEK MULTIRATE STEPS", std::to_string(timeStepRatio));
  return dt;
}

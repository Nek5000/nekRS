#include <cfloat>
#include "bcMap.hpp"
#include "findpts.hpp"
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

dlong computeNumInterpPoints(nrs_t *nrs)
{
  auto *mesh = nrs->meshV;
  auto fields = neknekSolveFields();

  if (fields.size() == 0) {
    return 0;
  }

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

} // namespace

bool neknekCoupled()
{
  auto solverFields = neknekSolveFields();

  std::set<int> intBIDs;
  std::ostringstream errorLogger;

  int coupled = 0;
  for (auto &&field : solverFields) {
    for (int bID = 1; bID <= bcMap::size(field); ++bID) {
      auto bcType = bcMap::id(bID, field);
      bool isInt = isIntBc(bcType, field);

      if (isInt) {
        intBIDs.insert(bID);
        coupled += 1;
      }

      // INT bid in one field -> INT bid in all fields
      if ((intBIDs.find(bID) != intBIDs.end()) && !isInt) {
        errorLogger << "ERROR: expected INT boundary condition on boundary id " << bID << " for field "
                    << field << "\n";
      }
    }
  }

  int errorLength = errorLogger.str().length();
  MPI_Allreduce(MPI_IN_PLACE, &errorLength, 1, MPI_INT, MPI_MAX, platform->comm.mpiCommParent);

  nrsCheck(errorLength > 0, platform->comm.mpiCommParent, EXIT_FAILURE, "%s\n", errorLogger.str().c_str());

  MPI_Allreduce(MPI_IN_PLACE, &coupled, 1, MPI_INT, MPI_MAX, platform->comm.mpiCommParent);
  return coupled > 0;
}

void neknek_t::reserveAllocation()
{
  this->fieldOffset_ = alignStride<dfloat>(this->npt_);
  this->o_pointMap_ = platform->device.malloc<dlong>(nrs->meshV->Nlocal);

  if (this->npt_) {
    this->o_U_ = platform->device.malloc<dfloat>(nrs->NVfields * this->fieldOffset_ * (this->nEXT_ + 1));
    if (this->Nscalar_) {
      this->o_S_ = platform->device.malloc<dfloat>(this->Nscalar_ * this->fieldOffset_ * (this->nEXT_ + 1));
    }
  }
}

void neknek_t::updateInterpPoints()
{
  // called in case of moving mesh ONLY
  if (!this->globalMovingMesh) {
    return;
  }

  auto *neknek = nrs->neknek;
  const dlong nsessions = this->nsessions_;
  const dlong sessionID = this->sessionID_;

  auto *mesh = nrs->meshV;

  // Setup findpts
  const dfloat tol = (sizeof(dfloat) == sizeof(double)) ? 5e-13 : 1e-6;
  constexpr dlong npt_max = 1;
  const dfloat bb_tol = 0.01;

  auto &device = platform->device.occaDevice();

  this->interpolator.reset();
  this->interpolator = std::make_shared<pointInterpolation_t>(nrs, bb_tol, tol, true, sessionID_, true);
  this->interpolator->setTimerLevel(TimerLevel::Basic);
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

  this->interpolator->setPoints(this->npt_, this->o_x_, this->o_y_, this->o_z_, this->o_session_);

  const auto verboseLevel = pointInterpolation_t::VerbosityLevel::Detailed;
  this->interpolator->find(verboseLevel);
}

void neknek_t::findIntPoints()
{
  const dlong nsessions = this->nsessions_;
  const dlong sessionID = this->sessionID_;

  auto *mesh = nrs->meshV;

  // Setup findpts
  const dfloat tol = (sizeof(dfloat) == sizeof(double)) ? 5e-13 : 1e-6;
  constexpr dlong npt_max = 1;
  const dfloat bb_tol = 0.01;

  auto &device = platform->device.occaDevice();

  this->interpolator.reset();
  this->interpolator = std::make_shared<pointInterpolation_t>(nrs, bb_tol, tol, true, sessionID_, true);
  this->interpolator->setTimerLevel(TimerLevel::Basic);
  this->interpolator->setTimerName("neknek_t::");

  auto numPoints = computeNumInterpPoints(nrs);
  this->npt_ = numPoints;
  this->reserveAllocation();

  std::vector<dfloat> neknekX(numPoints, 0.0);
  std::vector<dfloat> neknekY(numPoints, 0.0);
  std::vector<dfloat> neknekZ(numPoints, 0.0);
  std::vector<dlong> session(numPoints, 0.0);

  auto fields = neknekSolveFields();

  std::vector<dlong> pointMap(mesh->Nlocal, -1);

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
            session[ip] = sessionID;

            pointMap[idM] = ip;
            ++ip;
          }
        }
      }
    }
  }

  this->o_pointMap_.copyFrom(pointMap.data());

  this->interpolator->setPoints(numPoints, neknekX.data(), neknekY.data(), neknekZ.data(), session.data());

  const auto verboseLevel = pointInterpolation_t::VerbosityLevel::Detailed;
  this->interpolator->find(verboseLevel);

  this->o_x_ = platform->device.malloc<dfloat>(this->npt_, neknekX.data());
  this->o_y_ = platform->device.malloc<dfloat>(this->npt_, neknekY.data());
  this->o_z_ = platform->device.malloc<dfloat>(this->npt_, neknekZ.data());
  this->o_session_ = platform->device.malloc<dlong>(this->npt_, session.data());
}

void neknek_t::setup()
{
  dlong globalRank;
  MPI_Comm_rank(platform->comm.mpiCommParent, &globalRank);

  const int nsessions = this->nsessions_;
  if(platform->comm.mpiRank == 0) {
    printf("configuring neknek with %d sessions\n", nsessions);
    std::fflush(stdout);
  }

  nrsCheck(static_cast<bool>(platform->options.compareArgs("CONSTANT FLOW RATE", "TRUE")), 
           platform->comm.mpiComm, EXIT_FAILURE,
           "%s\n", "constant flow rate support not supported");

  if(nrs->pSolver) {
    nrsCheck(static_cast<bool>(nrs->pSolver->allNeumann && platform->options.compareArgs("LOWMACH", "TRUE")), 
             platform->comm.mpiComm, EXIT_FAILURE, 
             "%s\n", "variable p0th is not supported!");
  }

  std::vector<int> Nscalars(nsessions, 0);
  Nscalars[this->sessionID_] = nrs->Nscalar;

  MPI_Allreduce(MPI_IN_PLACE, Nscalars.data(), nsessions, MPI_INT, MPI_MAX, platform->comm.mpiCommParent);

  auto minNscalar = *std::min_element(Nscalars.begin(),Nscalars.end());

  bool allSame = std::all_of(Nscalars.begin(), Nscalars.end(), [minNscalar] (auto v) {return v == minNscalar;});

  if(platform->comm.mpiRank == 0 && !allSame){
    std::cout << "WARNING: Nscalar is not the same across all sessions -> using the minimum value: " << minNscalar << "\n";
  }

  this->Nscalar_ = minNscalar;

  int movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
  MPI_Allreduce(MPI_IN_PLACE, &movingMesh, 1, MPI_INT, MPI_MAX, platform->comm.mpiCommParent);
  this->globalMovingMesh = movingMesh;

  this->findIntPoints();
}

neknek_t::neknek_t(nrs_t *nrs_, dlong nsessions, dlong sessionID)
    : nsessions_(nsessions), sessionID_(sessionID), nrs(nrs_)
{
  nrs->neknek = this;
  if (nrs->cds) {
    nrs->cds->neknek = this;
  }

  this->nEXT_ = 1;
  if(!platform->options.getArgs("NEKNEK BOUNDARY EXT ORDER").empty())
    platform->options.getArgs("NEKNEK BOUNDARY EXT ORDER", this->nEXT_);

  // set boundary ext order to report to user, if not specified
  platform->options.setArgs("NEKNEK BOUNDARY EXT ORDER", std::to_string(this->nEXT_));

  this->coeffEXT.resize(this->nEXT_);
  this->o_coeffEXT = platform->device.malloc<dfloat>(this->nEXT_);

  this->setup();

  this->copyNekNekPointsKernel = platform->kernels.get("copyNekNekPoints");
  this->computeFluxKernel = platform->kernels.get("computeFlux");
  this->fixSurfaceFluxKernel = platform->kernels.get("fixSurfaceFlux");
}

void neknek_t::updateBoundary(int tstep, int stage)
{
  // do not invoke barrier -- this is performed later
  platform->timer.tic("neknek update boundary", 0);

  // do not invoke barrier in timer_t::tic
  platform->timer.tic("neknek sync", 0);
  MPI_Barrier(platform->comm.mpiCommParent);
  platform->timer.toc("neknek sync");
  this->tSync_ = platform->timer.query("neknek sync", "HOST:MAX");

  if (this->globalMovingMesh) {
    platform->timer.tic("neknek updateInterpPoints", 1);
    this->updateInterpPoints();
    platform->timer.toc("neknek updateInterpPoints");

    this->recomputePartition = true;
  }

  platform->timer.tic("neknek exchange", 1);

  this->interpolator->eval(nrs->NVfields, nrs->fieldOffset, nrs->o_U, this->fieldOffset_, this->o_U_);

  if (this->Nscalar_) {
    this->interpolator->eval(this->Nscalar_, nrs->fieldOffset, nrs->cds->o_S, this->fieldOffset_, this->o_S_);
  }

  // lag state, update timestepper coefficients and compute extrapolated state
  if (stage == 1) {
    auto *mesh = nrs->meshV;
    int extOrder = std::min(tstep, this->nEXT_);
    int bdfOrder = std::min(tstep, nrs->nBDF);
    nek::extCoeff(this->coeffEXT.data(), nrs->dt, extOrder, bdfOrder);

    for (int i = this->nEXT_; i > extOrder; i--) {
      this->coeffEXT[i - 1] = 0.0;
    }

    this->o_coeffEXT.copyFrom(this->coeffEXT.data(), this->nEXT_);

    for (int s = this->nEXT_ + 1; s > 1; s--) {
      auto N = nrs->NVfields * this->fieldOffset_;
      this->o_U_.copyFrom(this->o_U_, N, (s - 1) * N, (s - 2) * N);

      N = this->Nscalar_ * this->fieldOffset_;
      this->o_S_.copyFrom(this->o_S_, N, (s - 1) * N, (s - 2) * N);
    }

    auto o_Uold = this->o_U_ + this->fieldOffset_ * nrs->NVfields;
    auto o_Sold = this->o_S_ + this->fieldOffset_ * this->Nscalar_;

    if (this->npt_) {
      nrs->extrapolateKernel(this->npt_,
                             nrs->NVfields,
                             this->nEXT_,
                             this->fieldOffset_,
                             this->o_coeffEXT,
                             o_Uold,
                             this->o_U_);
    }

    if (this->Nscalar_ && this->npt_) {
      nrs->extrapolateKernel(this->npt_,
                             this->Nscalar_,
                             this->nEXT_,
                             this->fieldOffset_,
                             this->o_coeffEXT,
                             o_Sold,
                             this->o_S_);
    }
  }

  platform->timer.toc("neknek exchange");

  this->tExch_ = platform->timer.query("neknek exchange", "DEVICE:MAX");
  this->ratio_ = this->tSync_ / this->tExch_;

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

  const dfloat tol = (sizeof(dfloat) == sizeof(double)) ? 5e-13 : 1e-6;
  constexpr dlong npt_max = 1;
  const dfloat bb_tol = 0.01;
  auto pointInterp = pointInterpolation_t(nrs, bb_tol, tol, true, sessionID_, true);
  auto o_dist = pointInterp.distance();

  auto mesh = nrs->meshV;

  auto o_sess = platform->o_memPool.reserve<dlong>(nrs->fieldOffset);
  auto o_sumDist = platform->o_memPool.reserve<dfloat>(nrs->fieldOffset);
  auto o_found = platform->o_memPool.reserve<dfloat>(nrs->fieldOffset);
  auto o_interpDist = platform->o_memPool.reserve<dfloat>(nrs->fieldOffset);
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

    pointInterp.setPoints(mesh->Nlocal, mesh->o_x, mesh->o_y, mesh->o_z, o_sess);
    pointInterp.find(pointInterpolation_t::VerbosityLevel::None, true);

    auto &data = pointInterp.data();
    for (int n = 0; n < mesh->Nlocal; ++n) {
      found[n] = (data.code[n] == findpts::CODE_NOT_FOUND) ? 0.0 : 1.0;
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
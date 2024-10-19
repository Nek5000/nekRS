#include "platform.hpp"
#include "linAlg.hpp"
#include "nrs.hpp"
#include "udf.hpp"
#include "alignment.hpp"
#include "bcMap.hpp"
#include "bdry.hpp"

namespace
{

dfloat constantFlowScale = 0.0;

inline dfloat distance(dfloat x1, dfloat x2, dfloat y1, dfloat y2, dfloat z1, dfloat z2)
{
  const dfloat dist_x = x1 - x2;
  const dfloat dist_y = y1 - y2;
  const dfloat dist_z = z1 - z2;
  return std::sqrt(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);
}

void computeDirection(dfloat x1, dfloat x2, dfloat y1, dfloat y2, dfloat z1, dfloat z2, dfloat *direction)
{
  direction[0] = x1 - x2;
  direction[1] = y1 - y2;
  direction[2] = z1 - z2;

  const dfloat invMagnitude = 1 / distance(x1, x1, y1, y2, z1, z2);

  direction[0] *= invMagnitude;
  direction[1] *= invMagnitude;
  direction[2] *= invMagnitude;
}

dfloat lengthScale;
dfloat baseFlowRate;
dfloat currentFlowRate;
dfloat postCorrectionFlowRate;
dfloat flowRate;

int fromBID;
int toBID;
dfloat flowDirection[3];

void compute(nrs_t *nrs, double time)
{
  auto mesh = nrs->mesh;

  double flops = 0.0;

  platform->timer.tic("pressure rhs", 1);
  auto o_Prhs = platform->deviceMemoryPool.reserve<dfloat>(nrs->fieldOffset);
  auto o_lambda0 = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nlocal);
  {
    occa::memory o_gradPCoeff = platform->deviceMemoryPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);

    platform->linAlg->adyz(mesh->Nlocal, 1.0, nrs->o_rho, o_lambda0);

    nrs->wgradientVolumeKernel(mesh->Nelements,
                               mesh->o_vgeo,
                               mesh->o_D,
                               nrs->fieldOffset,
                               o_lambda0,
                               o_gradPCoeff);

    double flopsGrad = 6 * mesh->Np * mesh->Nq + 18 * mesh->Np;
    flopsGrad *= static_cast<double>(mesh->Nelements);
    flops += flopsGrad;

    nrs->computeFieldDotNormalKernel(mesh->Nlocal,
                                     nrs->fieldOffset,
                                     flowDirection[0],
                                     flowDirection[1],
                                     flowDirection[2],
                                     o_gradPCoeff,
                                     o_Prhs);

    flops += 5 * mesh->Nlocal;
  }
  platform->timer.toc("pressure rhs");

  platform->timer.tic("pressureSolve", 1);
  nrs->pSolver->solve(o_lambda0, o_NULL, o_Prhs, nrs->o_Pc);
  platform->timer.toc("pressureSolve");
  o_Prhs.free();

  // solve homogenous Stokes problem
  platform->timer.tic("velocity rhs", 1);
  occa::memory o_RhsVel = platform->deviceMemoryPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);
  {
    nrs->gradientVolumeKernel(mesh->Nelements,
                              mesh->o_vgeo,
                              mesh->o_D,
                              nrs->fieldOffset,
                              nrs->o_Pc,
                              o_RhsVel);

    double flopsGrad = 6 * mesh->Np * mesh->Nq + 18 * mesh->Np;
    flopsGrad *= static_cast<double>(mesh->Nelements);
    flops += flopsGrad;

    platform->linAlg->scaleMany(mesh->Nlocal, nrs->NVfields, nrs->fieldOffset, -1.0, o_RhsVel);

    occa::memory o_JwF = platform->deviceMemoryPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);
    o_JwF.copyFrom(mesh->o_LMM, mesh->Nlocal, 0 * nrs->fieldOffset, 0);
    o_JwF.copyFrom(mesh->o_LMM, mesh->Nlocal, 1 * nrs->fieldOffset, 0);
    o_JwF.copyFrom(mesh->o_LMM, mesh->Nlocal, 2 * nrs->fieldOffset, 0);

    for (int dim = 0; dim < nrs->NVfields; ++dim) {
      const dlong offset = dim * nrs->fieldOffset;
      const dfloat n_dim = flowDirection[dim];
      platform->linAlg->axpby(mesh->Nlocal, n_dim, o_JwF, 1.0, o_RhsVel, offset, offset);
    }
  }
  platform->timer.toc("velocity rhs");

  platform->timer.tic("velocitySolve", 1);

  o_lambda0.free();
  o_lambda0 = nrs->o_mue;
  auto o_lambda1 = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nlocal);
  platform->linAlg->axpby(mesh->Nlocal, nrs->g0 / nrs->dt[0], nrs->o_rho, 0.0, o_lambda1);

  if (nrs->uvwSolver) {
    nrs->uvwSolver->solve(o_lambda0, o_lambda1, o_RhsVel, nrs->o_Uc);
  } else {
    occa::memory o_Ucx = nrs->o_Uc + 0 * nrs->fieldOffset;
    occa::memory o_Ucy = nrs->o_Uc + 1 * nrs->fieldOffset;
    occa::memory o_Ucz = nrs->o_Uc + 2 * nrs->fieldOffset;
    nrs->uSolver->solve(o_lambda0, o_lambda1, o_RhsVel.slice(0 * nrs->fieldOffset), o_Ucx);
    nrs->vSolver->solve(o_lambda0, o_lambda1, o_RhsVel.slice(1 * nrs->fieldOffset), o_Ucy);
    nrs->wSolver->solve(o_lambda0, o_lambda1, o_RhsVel.slice(2 * nrs->fieldOffset), o_Ucz);
  }
  platform->timer.toc("velocitySolve");

  platform->flopCounter->add("ConstantFlowRate::compute", flops);
}

bool checkIfRecompute(nrs_t *nrs, int tstep)
{

  mesh_t *mesh = nrs->mesh;

  constexpr int nPropertyFields = 2;
  const dfloat TOL = 1e-6;

  bool adjustFlowRate = false;

  // did the properties change?
  occa::memory o_propDelta = platform->deviceMemoryPool.reserve<dfloat>(nPropertyFields * nrs->fieldOffset);
  platform->linAlg->axpbyzMany(mesh->Nlocal,
                               nPropertyFields,
                               nrs->fieldOffset,
                               1.0,
                               nrs->o_prop,
                               -1.0,
                               nrs->o_prevProp,
                               o_propDelta);

  const dfloat delta = platform->linAlg->norm2Many(mesh->Nlocal,
                                                   nPropertyFields,
                                                   nrs->fieldOffset,
                                                   o_propDelta,
                                                   platform->comm.mpiComm);

  if (delta > TOL) {
    adjustFlowRate = true;
    nrs->o_prevProp.copyFrom(nrs->o_prop, nPropertyFields * nrs->fieldOffset);
  }

  adjustFlowRate |= platform->options.compareArgs("MOVING MESH", "TRUE");
  adjustFlowRate |= tstep <= std::max(nrs->nEXT, nrs->nBDF);
  adjustFlowRate |= abs(nrs->dt[0] - nrs->dt[1]) > TOL;

  static dfloat prevFlowRate = 0;
  if (std::abs(flowRate - prevFlowRate) > TOL) {
    adjustFlowRate |= true;
    prevFlowRate = flowRate;
  }

  return adjustFlowRate;
}

bool checkIfRecomputeDirection(nrs_t *nrs, int tstep)
{
  return platform->options.compareArgs("MOVING MESH", "TRUE") || tstep < 2;
}

} // namespace

void nrs_t::flowRatePrintInfo(bool verboseInfo)
{
  auto mesh = this->mesh;

  if (platform->comm.mpiRank != 0) {
    return;
  }

  std::string flowRateType = "flowRate";

  dfloat currentRate = currentFlowRate;
  dfloat finalFlowRate = postCorrectionFlowRate;
  dfloat userSpecifiedFlowRate = flowRate * mesh->volume / lengthScale;

  dfloat err = std::abs(userSpecifiedFlowRate - finalFlowRate);

  dfloat scale = constantFlowScale; // rho * meanGradP

  if (!platform->options.compareArgs("CONSTANT FLOW RATE TYPE", "VOLUMETRIC")) {
    flowRateType = "uBulk";

    // put in bulk terms, instead of volumetric
    currentRate *= lengthScale / mesh->volume;
    finalFlowRate *= lengthScale / mesh->volume;
    userSpecifiedFlowRate = flowRate;
    err = std::abs(userSpecifiedFlowRate - finalFlowRate);
  }
  if (verboseInfo) {
    printf("flowRate  : %s0 %.2e  %s %.2e  err %.2e  scale %.5e\n",
           flowRateType.c_str(),
           currentRate,
           flowRateType.c_str(),
           finalFlowRate,
           err,
           scale);
  }
}

bool nrs_t::adjustFlowRate(int tstep, double time)
{
  double flops = 0.0;

  mesh_t *mesh = this->mesh;
  platform->options.getArgs("FLOW RATE", flowRate);

  const bool movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

  const bool X_aligned = platform->options.compareArgs("CONSTANT FLOW DIRECTION", "X");
  const bool Y_aligned = platform->options.compareArgs("CONSTANT FLOW DIRECTION", "Y");
  const bool Z_aligned = platform->options.compareArgs("CONSTANT FLOW DIRECTION", "Z");
  const bool directionAligned = X_aligned || Y_aligned || Z_aligned;

  nekrsCheck(!directionAligned,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "Flow direction is not aligned in (X,Y,Z)");

  const bool recomputeBaseFlowRate = checkIfRecompute(this, tstep);
  const bool recomputeDirection = checkIfRecomputeDirection(this, tstep);

  if (recomputeDirection) {
    if (directionAligned) {
      occa::memory o_coord;
      if (X_aligned) {
        o_coord = mesh->o_x;
        flowDirection[0] = 1.0;
        flowDirection[1] = 0.0;
        flowDirection[2] = 0.0;
      }
      if (Y_aligned) {
        o_coord = mesh->o_y;
        flowDirection[0] = 0.0;
        flowDirection[1] = 1.0;
        flowDirection[2] = 0.0;
      }
      if (Z_aligned) {
        o_coord = mesh->o_z;
        flowDirection[0] = 0.0;
        flowDirection[1] = 0.0;
        flowDirection[2] = 1.0;
      }

      const dfloat maxCoord = platform->linAlg->max(mesh->Nlocal, o_coord, platform->comm.mpiComm);
      const dfloat minCoord = platform->linAlg->min(mesh->Nlocal, o_coord, platform->comm.mpiComm);
      lengthScale = maxCoord - minCoord;
    } else {

      platform->options.getArgs("CONSTANT FLOW FROM BID", fromBID);
      platform->options.getArgs("CONSTANT FLOW TO BID", toBID);

      occa::memory o_centroid =
          platform->deviceMemoryPool.reserve<dfloat>(this->NVfields * mesh->Nelements * mesh->Nfaces);
      platform->linAlg->fill(mesh->Nelements * mesh->Nfaces * 3, 0.0, o_centroid);

      occa::memory o_counts = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nelements * mesh->Nfaces);
      platform->linAlg->fill(mesh->Nelements * mesh->Nfaces, 0.0, o_counts);

      this->computeFaceCentroidKernel(mesh->Nelements,
                                      fromBID,
                                      mesh->o_EToB,
                                      mesh->o_vmapM,
                                      mesh->o_x,
                                      mesh->o_y,
                                      mesh->o_z,
                                      o_centroid,
                                      o_counts);
      flops += 3 * mesh->Nlocal;

      dfloat NfacesContrib =
          platform->linAlg->sum(mesh->Nelements * mesh->Nfaces, o_counts, platform->comm.mpiComm);
      dfloat sumFaceAverages_x = platform->linAlg->sum(mesh->Nelements * mesh->Nfaces,
                                                       o_centroid,
                                                       platform->comm.mpiComm,
                                                       0 * mesh->Nelements * mesh->Nfaces);
      dfloat sumFaceAverages_y = platform->linAlg->sum(mesh->Nelements * mesh->Nfaces,
                                                       o_centroid,
                                                       platform->comm.mpiComm,
                                                       1 * mesh->Nelements * mesh->Nfaces);
      dfloat sumFaceAverages_z = platform->linAlg->sum(mesh->Nelements * mesh->Nfaces,
                                                       o_centroid,
                                                       platform->comm.mpiComm,
                                                       2 * mesh->Nelements * mesh->Nfaces);

      const dfloat centroidFrom_x = sumFaceAverages_x / NfacesContrib;
      const dfloat centroidFrom_y = sumFaceAverages_y / NfacesContrib;
      const dfloat centroidFrom_z = sumFaceAverages_z / NfacesContrib;

      platform->linAlg->fill(mesh->Nelements * mesh->Nfaces * 3, 0.0, o_centroid);
      platform->linAlg->fill(mesh->Nelements * mesh->Nfaces, 0.0, o_counts);
      this->computeFaceCentroidKernel(mesh->Nelements,
                                      toBID,
                                      mesh->o_EToB,
                                      mesh->o_vmapM,
                                      mesh->o_x,
                                      mesh->o_y,
                                      mesh->o_z,
                                      o_centroid,
                                      o_counts);

      flops += 3 * mesh->Nlocal;

      NfacesContrib = platform->linAlg->sum(mesh->Nelements * mesh->Nfaces, o_counts, platform->comm.mpiComm);
      sumFaceAverages_x = platform->linAlg->sum(mesh->Nelements * mesh->Nfaces,
                                                o_centroid,
                                                platform->comm.mpiComm,
                                                0 * mesh->Nelements * mesh->Nfaces);
      sumFaceAverages_y = platform->linAlg->sum(mesh->Nelements * mesh->Nfaces,
                                                o_centroid,
                                                platform->comm.mpiComm,
                                                1 * mesh->Nelements * mesh->Nfaces);
      sumFaceAverages_z = platform->linAlg->sum(mesh->Nelements * mesh->Nfaces,
                                                o_centroid,
                                                platform->comm.mpiComm,
                                                2 * mesh->Nelements * mesh->Nfaces);

      const dfloat centroidTo_x = sumFaceAverages_x / NfacesContrib;
      const dfloat centroidTo_y = sumFaceAverages_y / NfacesContrib;
      const dfloat centroidTo_z = sumFaceAverages_z / NfacesContrib;

      lengthScale =
          distance(centroidFrom_x, centroidTo_x, centroidFrom_y, centroidTo_y, centroidFrom_z, centroidTo_z);

      computeDirection(centroidFrom_x,
                       centroidTo_x,
                       centroidFrom_y,
                       centroidTo_y,
                       centroidFrom_z,
                       centroidTo_z,
                       flowDirection);
    }
  }

  if (recomputeBaseFlowRate) {
    if (platform->verbose && platform->comm.mpiComm == 0) {
      std::cout << "recomputing base flow rate\n";
    }

    auto getSolverData = [](elliptic *solver) {
      if (solver) {
        std::tuple<int, dfloat, dfloat, dfloat> val(solver->Niter(),
                                                    solver->initialResidual(),
                                                    solver->initialGuessResidual(),
                                                    solver->finalResidual());
        return val;
      } else {
        std::tuple<int, dfloat, dfloat, dfloat> val(0, 0, 0, 0);
        return val;
      }
    };

    const auto [NiterUVW, res00NormUVW, res0NormUVW, resNormUVW] = getSolverData(this->uvwSolver);
    const auto [NiterU, res00NormU, res0NormU, resNormU] = getSolverData(this->uSolver);
    const auto [NiterV, res00NormV, res0NormV, resNormV] = getSolverData(this->vSolver);
    const auto [NiterW, res00NormW, res0NormW, resNormW] = getSolverData(this->wSolver);
    const auto [NiterP, res00NormP, res0NormP, resNormP] = getSolverData(this->pSolver);

    compute(this, time);

    // restore norms + update iteration count
    auto setSolverData = [](elliptic *solver, int Niter, dfloat res00Norm, dfloat res0Norm, dfloat resNorm) {
      solver->Niter(solver->Niter() + Niter);
      solver->initialResidual(res00Norm);
      solver->initialGuessResidual(res0Norm);
      solver->finalResidual(resNorm);
    };

    if (this->uvwSolver) {
      setSolverData(this->uvwSolver, NiterUVW, res00NormUVW, res0NormUVW, resNormUVW);
    } else {
      setSolverData(this->uSolver, NiterU, res00NormU, res0NormU, resNormU);
      setSolverData(this->vSolver, NiterV, res00NormV, res0NormV, resNormV);
      setSolverData(this->wSolver, NiterW, res00NormW, res0NormW, resNormW);
    }
    setSolverData(this->pSolver, NiterP, res00NormP, res0NormP, resNormP);
  }

  occa::memory o_currentFlowRate = platform->deviceMemoryPool.reserve<dfloat>(this->fieldOffset);

  this->computeFieldDotNormalKernel(mesh->Nlocal,
                                    this->fieldOffset,
                                    flowDirection[0],
                                    flowDirection[1],
                                    flowDirection[2],
                                    this->o_U,
                                    o_currentFlowRate);

  flops += 5 * mesh->Nlocal;

  // scale by mass matrix
  platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_LMM, o_currentFlowRate);

  currentFlowRate =
      platform->linAlg->sum(mesh->Nlocal, o_currentFlowRate, platform->comm.mpiComm) / lengthScale;

  if (recomputeBaseFlowRate) {
    occa::memory o_baseFlowRate = platform->deviceMemoryPool.reserve<dfloat>(this->fieldOffset);
    this->computeFieldDotNormalKernel(mesh->Nlocal,
                                      this->fieldOffset,
                                      flowDirection[0],
                                      flowDirection[1],
                                      flowDirection[2],
                                      this->o_Uc,
                                      o_baseFlowRate);
    flops += 5 * mesh->Nlocal;

    // scale by mass matrix
    platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_LMM, o_baseFlowRate);
    baseFlowRate = platform->linAlg->sum(mesh->Nlocal, o_baseFlowRate, platform->comm.mpiComm) / lengthScale;
  }

  // user specifies a mean velocity, not volumetric flow rate
  dfloat volumetricFlowRate = flowRate * mesh->volume / lengthScale;
  if (platform->options.compareArgs("CONSTANT FLOW RATE TYPE", "VOLUMETRIC")) {
    volumetricFlowRate = flowRate;
  }

  const dfloat deltaFlowRate = volumetricFlowRate - currentFlowRate;

  constantFlowScale = deltaFlowRate / baseFlowRate;

  // add corrections
  platform->linAlg->axpbyMany(mesh->Nlocal,
                              this->NVfields,
                              this->fieldOffset,
                              constantFlowScale,
                              this->o_Uc,
                              1.0,
                              this->o_U);

  platform->linAlg->axpby(mesh->Nlocal, constantFlowScale, this->o_Pc, 1.0, this->o_P);

  // compute flow rate after correction as diagnostic
  this->computeFieldDotNormalKernel(mesh->Nlocal,
                                    this->fieldOffset,
                                    flowDirection[0],
                                    flowDirection[1],
                                    flowDirection[2],
                                    this->o_U,
                                    o_currentFlowRate);

  flops += 5 * mesh->Nlocal;

  // scale by mass matrix
  platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_LMM, o_currentFlowRate);

  postCorrectionFlowRate =
      platform->linAlg->sum(mesh->Nlocal, o_currentFlowRate, platform->comm.mpiComm) / lengthScale;

  platform->flopCounter->add("ConstantFlowRate::adjust", flops);

  return recomputeBaseFlowRate;
}

dfloat nrs_t::flowRatescaleFactor()
{
  return constantFlowScale;
}

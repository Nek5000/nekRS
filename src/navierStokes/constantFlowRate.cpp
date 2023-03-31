#include "constantFlowRate.hpp"
#include "linAlg.hpp"
#include "nrs.hpp"
#include "udf.hpp"
#include <limits>
#include "alignment.hpp"
#include "bcMap.hpp"
#include "bdry.hpp"

namespace {

dfloat constantFlowScale = 0.0;

inline dfloat distance(
    dfloat x1, dfloat x2, dfloat y1, dfloat y2, dfloat z1, dfloat z2) {
  const dfloat dist_x = x1 - x2;
  const dfloat dist_y = y1 - y2;
  const dfloat dist_z = z1 - z2;
  return std::sqrt(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);
}

void computeDirection(dfloat x1,
    dfloat x2,
    dfloat y1,
    dfloat y2,
    dfloat z1,
    dfloat z2,
    dfloat *direction) {
  direction[0] = x1 - x2;
  direction[1] = y1 - y2;
  direction[2] = z1 - z2;

  const dfloat magnitude = distance(x1, x1, y1, y2, z1, z2);

  direction[0] /= magnitude;
  direction[1] /= magnitude;
  direction[2] /= magnitude;
}

dfloat lengthScale;
dfloat baseFlowRate;
dfloat currentFlowRate;
dfloat postCorrectionFlowRate;
dfloat flowRate;

int fromBID;
int toBID;
dfloat flowDirection[3];

void compute(nrs_t *nrs, dfloat time) {

  constexpr int ndim = 3;
  mesh_t *mesh = nrs->meshV;

  double flops = 0.0;

  platform->timer.tic("pressure rhs", 1);
  occa::memory &o_gradPCoeff = platform->o_mempool.slice0;
  occa::memory &o_Prhs = platform->o_mempool.slice3;

  nrs->setEllipticCoeffPressureKernel(
      mesh->Nlocal, nrs->fieldOffset, nrs->o_rho, nrs->o_ellipticCoeff);

  nrs->wgradientVolumeKernel(mesh->Nelements,
      mesh->o_vgeo,
      mesh->o_D,
      nrs->fieldOffset,
      nrs->o_ellipticCoeff,
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
  platform->timer.toc("pressure rhs");

  platform->timer.tic("pressureSolve", 1);
  ellipticSolve(nrs->pSolver, o_Prhs, nrs->o_Pc);
  platform->timer.toc("pressureSolve");

  // solve homogenous Stokes problem
  platform->timer.tic("velocity rhs", 1);
  occa::memory &o_RhsVel = platform->o_mempool.slice0;
  nrs->gradientVolumeKernel(mesh->Nelements,
      mesh->o_vgeo,
      mesh->o_D,
      nrs->fieldOffset,
      nrs->o_Pc,
      o_RhsVel); 

  flopsGrad = 6 * mesh->Np * mesh->Nq + 18 * mesh->Np;
  flopsGrad *= static_cast<double>(mesh->Nelements);
  flops += flopsGrad;

  platform->linAlg->scaleMany(
      mesh->Nlocal, nrs->NVfields, nrs->fieldOffset, -1.0, o_RhsVel);

  occa::memory &o_BF = platform->o_mempool.slice3;
  o_BF.copyFrom(mesh->o_LMM,
      mesh->Nlocal * sizeof(dfloat),
      0 * nrs->fieldOffset * sizeof(dfloat),
      0);
  o_BF.copyFrom(mesh->o_LMM,
      mesh->Nlocal * sizeof(dfloat),
      1 * nrs->fieldOffset * sizeof(dfloat),
      0);
  o_BF.copyFrom(mesh->o_LMM,
      mesh->Nlocal * sizeof(dfloat),
      2 * nrs->fieldOffset * sizeof(dfloat),
      0);

  for (int dim = 0; dim < ndim; ++dim) {
    const dlong offset = dim * nrs->fieldOffset;
    const dfloat n_dim = flowDirection[dim];
    platform->linAlg->axpby(
        mesh->Nlocal, n_dim, o_BF, 1.0, o_RhsVel, offset, offset);
  }

  platform->timer.toc("velocity rhs");

  platform->timer.tic("velocitySolve", 1);
  nrs->setEllipticCoeffKernel(mesh->Nlocal,
      nrs->g0 * nrs->idt,
      0 * nrs->fieldOffset,
      nrs->fieldOffset,
      0,
      nrs->o_mue,
      nrs->o_rho,
      o_NULL,
      nrs->o_ellipticCoeff);

  if (nrs->uvwSolver) {
    ellipticSolve(nrs->uvwSolver, o_RhsVel, nrs->o_Uc);
  } else {
    occa::memory o_Ucx = nrs->o_Uc + (0 * sizeof(dfloat)) * nrs->fieldOffset;
    occa::memory o_Ucy = nrs->o_Uc + (1 * sizeof(dfloat)) * nrs->fieldOffset;
    occa::memory o_Ucz = nrs->o_Uc + (2 * sizeof(dfloat)) * nrs->fieldOffset;
    ellipticSolve(nrs->uSolver, platform->o_mempool.slice0, o_Ucx);
    ellipticSolve(nrs->vSolver, platform->o_mempool.slice1, o_Ucy);
    ellipticSolve(nrs->wSolver, platform->o_mempool.slice2, o_Ucz);
  }
  platform->timer.toc("velocitySolve");

  platform->flopCounter->add("ConstantFlowRate::compute", flops);
}

bool checkIfRecompute(nrs_t *nrs, int tstep) {

  mesh_t *mesh = nrs->meshV;

  constexpr int nPropertyFields = 2;
  const dfloat TOL = 1e-8;

  bool adjustFlowRate = false;

  platform->linAlg->axpbyzMany(mesh->Nlocal,
      nPropertyFields,
      nrs->fieldOffset,
      1.0,
      nrs->o_prop,
      -1.0,
      nrs->o_prevProp,
      platform->o_mempool.slice0);

  const dfloat delta = platform->linAlg->norm2Many(mesh->Nlocal,
      nPropertyFields,
      nrs->fieldOffset,
      platform->o_mempool.slice0,
      platform->comm.mpiComm);

  if (delta > TOL) {
    adjustFlowRate = true;
    nrs->o_prevProp.copyFrom(
        nrs->o_prop, nPropertyFields * nrs->fieldOffset * sizeof(dfloat));
  }

  adjustFlowRate |= platform->options.compareArgs("MOVING MESH", "TRUE");
  adjustFlowRate |= tstep <= std::max(nrs->nEXT, nrs->nBDF);
  adjustFlowRate |= abs(nrs->dt[0] - nrs->dt[1]) > TOL;
  return adjustFlowRate;
}

bool checkIfRecomputeDirection(nrs_t *nrs, int tstep) {
  return platform->options.compareArgs("MOVING MESH", "TRUE") || tstep < 2;
}

} // namespace


namespace ConstantFlowRate {

bool adjust(nrs_t *nrs, int tstep, dfloat time) {

  double flops = 0.0;

  constexpr int ndim = 3;
  mesh_t *mesh = nrs->meshV;
  platform->options.getArgs("FLOW RATE", flowRate);

  const bool movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

  const bool X_aligned =
      platform->options.compareArgs("CONSTANT FLOW DIRECTION", "X");
  const bool Y_aligned =
      platform->options.compareArgs("CONSTANT FLOW DIRECTION", "Y");
  const bool Z_aligned =
      platform->options.compareArgs("CONSTANT FLOW DIRECTION", "Z");
  const bool directionAligned = X_aligned || Y_aligned || Z_aligned;

  nrsCheck(!directionAligned, platform->comm.mpiComm, EXIT_FAILURE,
           "Flow direction is not aligned in (X,Y,Z)\n", "");

  const bool recomputeBaseFlowRate =
      checkIfRecompute(nrs, tstep);
  const bool recomputeDirection =
      checkIfRecomputeDirection(nrs, tstep);

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

      const dfloat maxCoord =
          platform->linAlg->max(mesh->Nlocal, o_coord, platform->comm.mpiComm);
      const dfloat minCoord =
          platform->linAlg->min(mesh->Nlocal, o_coord, platform->comm.mpiComm);
      lengthScale = maxCoord - minCoord;
    } else {

      platform->options.getArgs("CONSTANT FLOW FROM BID", fromBID);
      platform->options.getArgs("CONSTANT FLOW TO BID", toBID);

      occa::memory o_centroid = platform->o_mempool.slice0;
      occa::memory o_counts = platform->o_mempool.slice3;
      platform->linAlg->fill(
          mesh->Nelements * mesh->Nfaces * 3, 0.0, o_centroid);
      platform->linAlg->fill(mesh->Nelements * mesh->Nfaces, 0.0, o_counts);
      nrs->computeFaceCentroidKernel(mesh->Nelements,
          fromBID,
          mesh->o_EToB,
          mesh->o_vmapM,
          mesh->o_x,
          mesh->o_y,
          mesh->o_z,
          o_centroid,
          o_counts);
      flops += 3 * mesh->Nlocal;

      dfloat NfacesContrib = platform->linAlg->sum(
          mesh->Nelements * mesh->Nfaces, o_counts, platform->comm.mpiComm);
      dfloat sumFaceAverages_x =
          platform->linAlg->sum(mesh->Nelements * mesh->Nfaces,
              o_centroid,
              platform->comm.mpiComm,
              0 * mesh->Nelements * mesh->Nfaces);
      dfloat sumFaceAverages_y =
          platform->linAlg->sum(mesh->Nelements * mesh->Nfaces,
              o_centroid,
              platform->comm.mpiComm,
              1 * mesh->Nelements * mesh->Nfaces);
      dfloat sumFaceAverages_z =
          platform->linAlg->sum(mesh->Nelements * mesh->Nfaces,
              o_centroid,
              platform->comm.mpiComm,
              2 * mesh->Nelements * mesh->Nfaces);

      const dfloat centroidFrom_x = sumFaceAverages_x / NfacesContrib;
      const dfloat centroidFrom_y = sumFaceAverages_y / NfacesContrib;
      const dfloat centroidFrom_z = sumFaceAverages_z / NfacesContrib;

      platform->linAlg->fill(
          mesh->Nelements * mesh->Nfaces * 3, 0.0, o_centroid);
      platform->linAlg->fill(mesh->Nelements * mesh->Nfaces, 0.0, o_counts);
      nrs->computeFaceCentroidKernel(mesh->Nelements,
          toBID,
          mesh->o_EToB,
          mesh->o_vmapM,
          mesh->o_x,
          mesh->o_y,
          mesh->o_z,
          o_centroid,
          o_counts);

      flops += 3 * mesh->Nlocal;

      NfacesContrib = platform->linAlg->sum(
          mesh->Nelements * mesh->Nfaces, o_counts, platform->comm.mpiComm);
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

      lengthScale = distance(centroidFrom_x,
          centroidTo_x,
          centroidFrom_y,
          centroidTo_y,
          centroidFrom_z,
          centroidTo_z);

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
    int NiterU, NiterV, NiterW, NiterUVW, NiterP;
    double res00NormU, res00NormV, res00NormW, res00NormUVW, res00NormP;
    double res0NormU, res0NormV, res0NormW, res0NormUVW, res0NormP;
    double resNormU, resNormV, resNormW, resNormUVW, resNormP;
    if(nrs->uvwSolver){
      NiterUVW = nrs->uvwSolver->Niter;
      res00NormUVW = nrs->uvwSolver->res00Norm;
      res0NormUVW = nrs->uvwSolver->res0Norm;
      resNormUVW = nrs->uvwSolver->resNorm;
    }
    else {
      NiterU = nrs->uSolver->Niter;
      res00NormU = nrs->uSolver->res00Norm;
      res0NormU = nrs->uSolver->res0Norm;
      resNormU = nrs->uSolver->resNorm;

      NiterV = nrs->vSolver->Niter;
      res00NormV = nrs->vSolver->res00Norm;
      res0NormV = nrs->vSolver->res0Norm;
      resNormV = nrs->vSolver->resNorm;

      NiterW = nrs->wSolver->Niter;
      res00NormW = nrs->wSolver->res00Norm;
      res0NormW = nrs->wSolver->res0Norm;
      resNormW = nrs->wSolver->resNorm;

    }

    NiterP = nrs->pSolver->Niter;
    res00NormP = nrs->pSolver->res00Norm;
    res0NormP = nrs->pSolver->res0Norm;
    resNormP = nrs->pSolver->resNorm;

    compute(nrs, time);

    if(nrs->uvwSolver){
      nrs->uvwSolver->Niter += NiterUVW;
      nrs->uvwSolver->res00Norm = res00NormUVW;
      nrs->uvwSolver->res0Norm = res0NormUVW;
      nrs->uvwSolver->resNorm = resNormUVW;
    }
    else {
      nrs->uSolver->Niter += NiterU;
      nrs->uSolver->res00Norm = res00NormU;
      nrs->uSolver->res0Norm = res0NormU;
      nrs->uSolver->resNorm = resNormU;

      nrs->vSolver->Niter += NiterV;
      nrs->vSolver->res00Norm = res00NormV;
      nrs->vSolver->res0Norm = res0NormV;
      nrs->vSolver->resNorm = resNormV;

      nrs->wSolver->Niter += NiterW;
      nrs->wSolver->res00Norm = res00NormW;
      nrs->wSolver->res0Norm = res0NormW;
      nrs->wSolver->resNorm = resNormW;
    }
    nrs->pSolver->Niter += NiterP;
    nrs->pSolver->res00Norm = res00NormP;
    nrs->pSolver->res0Norm = res0NormP;
    nrs->pSolver->resNorm = resNormP;
  }

  occa::memory &o_currentFlowRate = platform->o_mempool.slice0;
  occa::memory &o_baseFlowRate = platform->o_mempool.slice1;

  nrs->computeFieldDotNormalKernel(mesh->Nlocal,
      nrs->fieldOffset,
      flowDirection[0],
      flowDirection[1],
      flowDirection[2],
      nrs->o_U,
      o_currentFlowRate);

  flops += 5 * mesh->Nlocal;

  // scale by mass matrix
  platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_LMM, o_currentFlowRate);

  currentFlowRate =
      platform->linAlg->sum(
          mesh->Nlocal, o_currentFlowRate, platform->comm.mpiComm) /
      lengthScale;

  if (recomputeBaseFlowRate) {
    nrs->computeFieldDotNormalKernel(mesh->Nlocal,
        nrs->fieldOffset,
        flowDirection[0],
        flowDirection[1],
        flowDirection[2],
        nrs->o_Uc,
        o_baseFlowRate);
    flops += 5 * mesh->Nlocal;

    // scale by mass matrix
    platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_LMM, o_baseFlowRate);
    baseFlowRate = platform->linAlg->sum(
                       mesh->Nlocal, o_baseFlowRate, platform->comm.mpiComm) /
                   lengthScale;
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
      nrs->NVfields,
      nrs->fieldOffset,
      constantFlowScale,
      nrs->o_Uc,
      1.0,
      nrs->o_U);

  platform->linAlg->axpby(mesh->Nlocal, constantFlowScale, nrs->o_Pc, 1.0, nrs->o_P);

  // compute flow rate after correction as diagnostic
  nrs->computeFieldDotNormalKernel(mesh->Nlocal,
      nrs->fieldOffset,
      flowDirection[0],
      flowDirection[1],
      flowDirection[2],
      nrs->o_U,
      o_currentFlowRate);

  flops += 5 * mesh->Nlocal;

  // scale by mass matrix
  platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_LMM, o_currentFlowRate);

  postCorrectionFlowRate =
      platform->linAlg->sum(
          mesh->Nlocal, o_currentFlowRate, platform->comm.mpiComm) /
      lengthScale;

  platform->flopCounter->add("ConstantFlowRate::adjust", flops);

  return recomputeBaseFlowRate;
}

dfloat scaleFactor(){
  return constantFlowScale;
}


void printInfo(mesh_t* mesh, bool verboseInfo)
{
  if(platform->comm.mpiRank != 0) return;

  std::string flowRateType = "flowRate";

  dfloat currentRate = currentFlowRate;
  dfloat finalFlowRate = postCorrectionFlowRate;
  dfloat userSpecifiedFlowRate = flowRate * mesh->volume / lengthScale;

  dfloat err = std::abs(userSpecifiedFlowRate - finalFlowRate);

  // scale is invariant to uBulk/volumetric flow rate, since it's a unitless ratio
  dfloat scale = constantFlowScale;

  if(!platform->options.compareArgs("CONSTANT FLOW RATE TYPE", "VOLUMETRIC")){
    flowRateType = "uBulk";

    // put in bulk terms, instead of volumetric
    currentRate *= lengthScale / mesh->volume;
    finalFlowRate *= lengthScale / mesh->volume;
    userSpecifiedFlowRate = flowRate;
    err = std::abs(userSpecifiedFlowRate - finalFlowRate);
  }
  if(verboseInfo) 
    printf("flowRate : %s0 %.2e  %s %.2e  err %.2e  scale %.5e\n",
      flowRateType.c_str(), currentRate, flowRateType.c_str(), finalFlowRate, err, scale);
}

} // namespace ConstantFlowRate

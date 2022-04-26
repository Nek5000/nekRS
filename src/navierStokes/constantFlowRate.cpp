#include "constantFlowRate.hpp"
#include "linAlg.hpp"
#include "nrs.hpp"
#include "udf.hpp"
#include <limits>
#include "alignment.hpp"
#include "bcMap.hpp"
#include "bdry.hpp"

namespace {
static dfloat constantFlowScale = 0.0;
inline dfloat distance(
    dfloat x1, dfloat x2, dfloat y1, dfloat y2, dfloat z1, dfloat z2) {
  const dfloat dist_x = x1 - x2;
  const dfloat dist_y = y1 - y2;
  const dfloat dist_z = z1 - z2;
  return std::sqrt(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);
}
inline void computeDirection(dfloat x1,
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

static dfloat lengthScale;
static dfloat baseFlowRate;
static dfloat currentFlowRate;
static dfloat postCorrectionFlowRate;
static dfloat flowRate;

static int fromBID;
static int toBID;
static dfloat flowDirection[3];

} // namespace

namespace ConstantFlowRate {

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

bool apply(nrs_t *nrs, int tstep, dfloat time) {

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

  if (!directionAligned) {
    if (platform->comm.mpiRank == 0)
      printf("Flow direction is not aligned in (X,Y,Z).\n"
             "Currently, only (X,Y,Z) aligned flow directions are supported in "
             "the "
             "constant flow rate driver.\n");
    ABORT(1);
  }

  const bool recomputeBaseFlowRate =
      ConstantFlowRate::checkIfRecompute(nrs, tstep);
  const bool recomputeDirection =
      ConstantFlowRate::checkIfRecomputeDirection(nrs, tstep);

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

    ConstantFlowRate::compute(nrs, lengthScale, time);

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

  // vx += scale * vxc
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

  platform->flopCounter->add("ConstantFlowRate::apply", flops);

  return recomputeBaseFlowRate;
}

dfloat scaleFactor(){
  return constantFlowScale;
}

void compute(nrs_t *nrs, double lengthScale, dfloat time) {

  constexpr int ndim = 3;
  mesh_t *mesh = nrs->meshV;

  double flops = 0.0;

  platform->timer.tic("pressureSolve", 1);
  {
    platform->timer.tic("pressure rhs", 1);
    occa::memory &o_gradPCoeff = platform->o_mempool.slice0;
    occa::memory &o_Prhs = platform->o_mempool.slice3;

    nrs->setEllipticCoeffPressureKernel(
        mesh->Nlocal, nrs->fieldOffset, nrs->o_rho, nrs->o_ellipticCoeff);

    nrs->gradientVolumeKernel(mesh->Nelements,
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

    // enforce Dirichlet BCs
    platform->linAlg->fill(nrs->fieldOffset,
        -1.0*std::numeric_limits<dfloat>::max(),
        platform->o_mempool.slice6);
    for (int sweep = 0; sweep < 2; sweep++) {
      nrs->pressureDirichletBCKernel(mesh->Nelements,
          time,
          nrs->fieldOffset,
          mesh->o_sgeo,
          mesh->o_x,
          mesh->o_y,
          mesh->o_z,
          mesh->o_vmapM,
          mesh->o_EToB,
          nrs->o_EToB,
          nrs->o_usrwrk,
          nrs->o_U,
          platform->o_mempool.slice6);

      // take care of Neumann-Dirichlet shared edges across elements
      if (sweep == 0)
        oogs::startFinish(platform->o_mempool.slice6,
            1,
            nrs->fieldOffset,
            ogsDfloat,
            ogsMax,
            nrs->gsh);
      if (sweep == 1)
        oogs::startFinish(platform->o_mempool.slice6,
            1,
            nrs->fieldOffset,
            ogsDfloat,
            ogsMin,
            nrs->gsh);
    }

    if (nrs->pSolver->Nmasked)
      nrs->maskCopyKernel(nrs->pSolver->Nmasked,
          0,
          nrs->pSolver->o_maskIds,
          platform->o_mempool.slice6,
          nrs->o_Pc);

    platform->timer.toc("pressure rhs");
    ellipticSolve(nrs->pSolver, o_Prhs, nrs->o_Pc);
  }
  platform->timer.toc("pressureSolve");

  platform->timer.tic("velocitySolve", 1);
  {
    platform->timer.tic("velocity rhs", 1);
    nrs->setEllipticCoeffKernel(mesh->Nlocal,
        nrs->g0 * nrs->idt,
        0 * nrs->fieldOffset,
        nrs->fieldOffset,
        nrs->o_mue,
        nrs->o_rho,
        nrs->o_ellipticCoeff);
    occa::memory &o_RhsVel = platform->o_mempool.slice0;
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
    nrs->gradientVolumeKernel(mesh->Nelements,
        mesh->o_vgeo,
        mesh->o_D,
        nrs->fieldOffset,
        nrs->o_Pc,
        o_RhsVel // <- rhs = -\grad{P_c}
    );

    double flopsGrad = 6 * mesh->Np * mesh->Nq + 18 * mesh->Np;
    flopsGrad *= static_cast<double>(mesh->Nelements);
    flops += flopsGrad;

    // rhs = -\grad{P_c} + BF n_i
    platform->linAlg->scaleMany(
        mesh->Nlocal, nrs->NVfields, nrs->fieldOffset, -1.0, o_RhsVel);

    for (int dim = 0; dim < ndim; ++dim) {
      const dlong offset = dim * nrs->fieldOffset;
      const dfloat n_dim = flowDirection[dim];
      platform->linAlg->axpby(
          mesh->Nlocal, n_dim, o_BF, 1.0, o_RhsVel, offset, offset);
    }

    if (nrs->uvwSolver) {
      if(nrs->uvwSolver->Nmasked)
        nrs->maskKernel(
            nrs->uvwSolver->Nmasked, nrs->uvwSolver->o_maskIds, o_RhsVel);
    } else {
      if (nrs->uSolver->Nmasked)
        nrs->maskKernel(nrs->uSolver->Nmasked,
            nrs->uSolver->o_maskIds,
            o_RhsVel + 0 * nrs->fieldOffset * sizeof(dfloat));
      if (nrs->vSolver->Nmasked)
        nrs->maskKernel(nrs->vSolver->Nmasked,
            nrs->vSolver->o_maskIds,
            o_RhsVel + 1 * nrs->fieldOffset * sizeof(dfloat));
      if (nrs->wSolver->Nmasked)
        nrs->maskKernel(nrs->wSolver->Nmasked,
            nrs->wSolver->o_maskIds,
            o_RhsVel + 2 * nrs->fieldOffset * sizeof(dfloat));
    }

    platform->linAlg->fill(nrs->NVfields * nrs->fieldOffset,
        -1.0*std::numeric_limits<dfloat>::max(),
        platform->o_mempool.slice3);
    for (int sweep = 0; sweep < 2; sweep++) {

      nrs->velocityDirichletBCKernel(mesh->Nelements,
          nrs->fieldOffset,
          time,
          mesh->o_sgeo,
          nrs->o_zeroNormalMaskVelocity,
          mesh->o_x,
          mesh->o_y,
          mesh->o_z,
          mesh->o_vmapM,
          mesh->o_EToB,
          nrs->o_EToB,
          nrs->o_usrwrk,
          nrs->o_Uc,
          platform->o_mempool.slice3);

      nrs->velocityNeumannBCKernel(mesh->Nelements,
          nrs->fieldOffset,
          mesh->o_sgeo,
          mesh->o_vmapM,
          mesh->o_EToB,
          nrs->o_EToB,
          time,
          mesh->o_x,
          mesh->o_y,
          mesh->o_z,
          nrs->o_usrwrk,
          nrs->o_U,
          platform->o_mempool.slice3);

      // take care of Neumann-Dirichlet shared edges across elements
      if (sweep == 0)
        oogs::startFinish(platform->o_mempool.slice3,
            nrs->NVfields,
            nrs->fieldOffset,
            ogsDfloat,
            ogsMax,
            nrs->gsh);
      if (sweep == 1)
        oogs::startFinish(platform->o_mempool.slice3,
            nrs->NVfields,
            nrs->fieldOffset,
            ogsDfloat,
            ogsMin,
            nrs->gsh);
    }
    if (nrs->uvwSolver) {

      if (nrs->uvwSolver->Nmasked)
        nrs->maskCopyKernel(nrs->uvwSolver->Nmasked,
            0 * nrs->fieldOffset,
            nrs->uvwSolver->o_maskIds,
            platform->o_mempool.slice3,
            nrs->o_Uc);
      if (bcMap::unalignedBoundary(mesh->cht, "velocity")) {
        applyZeroNormalMask(nrs, nrs->uvwSolver->o_EToB, nrs->o_zeroNormalMaskVelocity, nrs->o_Uc);
      }
    } else {
      if (nrs->uSolver->Nmasked)
        nrs->maskCopyKernel(nrs->uSolver->Nmasked,
            0 * nrs->fieldOffset,
            nrs->uSolver->o_maskIds,
            platform->o_mempool.slice3,
            nrs->o_Uc);
      if (nrs->vSolver->Nmasked)
        nrs->maskCopyKernel(nrs->vSolver->Nmasked,
            1 * nrs->fieldOffset,
            nrs->vSolver->o_maskIds,
            platform->o_mempool.slice3,
            nrs->o_Uc);
      if (nrs->wSolver->Nmasked)
        nrs->maskCopyKernel(nrs->wSolver->Nmasked,
            2 * nrs->fieldOffset,
            nrs->wSolver->o_maskIds,
            platform->o_mempool.slice3,
            nrs->o_Uc);
    }
    platform->timer.toc("velocity rhs");

    if (nrs->uvwSolver) {
      ellipticSolve(nrs->uvwSolver, o_RhsVel, nrs->o_Uc);
    } else {
      occa::memory o_Ucx = nrs->o_Uc + 0 * nrs->fieldOffset * sizeof(dfloat);
      occa::memory o_Ucy = nrs->o_Uc + 1 * nrs->fieldOffset * sizeof(dfloat);
      occa::memory o_Ucz = nrs->o_Uc + 2 * nrs->fieldOffset * sizeof(dfloat);
      ellipticSolve(nrs->uSolver, platform->o_mempool.slice0, o_Ucx);
      ellipticSolve(nrs->vSolver, platform->o_mempool.slice1, o_Ucy);
      ellipticSolve(nrs->wSolver, platform->o_mempool.slice2, o_Ucz);
    }
  }
  platform->timer.toc("velocitySolve");

  platform->flopCounter->add("ConstantFlowRate::compute", flops);

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
    printf("  flowRate : %s0 %.2e  %s %.2e  err %.2e  scale %.2e\n",
      flowRateType.c_str(), currentRate, flowRateType.c_str(), finalFlowRate, err, scale);
}

} // namespace ConstantFlowRate

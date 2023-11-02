#include "nrs.hpp"
#include "neknek.hpp"
#include "avm.hpp"
#include "cfl.hpp"
#include "constantFlowRate.hpp"
#include "nekInterfaceAdapter.hpp"
#include "timeStepper.hpp"
#include "tombo.hpp"
#include "subCycling.hpp"
#include "udf.hpp"
#include "bcMap.hpp"
#include "applyDirichlet.hpp"
#include "bdry.hpp"
#include "Urst.hpp"

static void lagFields(nrs_t *nrs)
{
  // lag velocity
  for (int s = std::max(nrs->nBDF, nrs->nEXT); s > 1; s--) {
    const auto N = nrs->NVfields * nrs->fieldOffset;
    nrs->o_U.copyFrom(nrs->o_U, N, (s - 1) * N, (s - 2) * N);
  }

  // lag scalars
  if (nrs->Nscalar) {
    auto cds = nrs->cds;
    if (cds->anyEllipticSolver) {
      for (int s = std::max(cds->nBDF, cds->nEXT); s > 1; s--) {
        const auto N = cds->fieldOffsetSum;
        cds->o_S.copyFrom(cds->o_S, N, (s - 1) * N, (s - 2) * N);
      }
    }
  }

  // lag mesh velocity
  const bool movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
  if (movingMesh) {
    auto mesh = nrs->_mesh;
    for (int s = std::max(nrs->nEXT, mesh->nAB); s > 1; s--) {
      const auto N = nrs->NVfields * nrs->fieldOffset;
      mesh->o_U.copyFrom(mesh->o_U, N, (s - 1) * N, (s - 2) * N);
    }
  }
}

static void extrapolate(nrs_t *nrs)
{
  mesh_t *mesh = nrs->meshV;

  {
    nrs->extrapolateKernel(mesh->Nlocal,
                           nrs->NVfields,
                           nrs->nEXT,
                           nrs->fieldOffset,
                           nrs->o_coeffEXT,
                           nrs->o_U,
                           nrs->o_Ue);

    if (nrs->flow && platform->options.compareArgs("LOWMACH", "TRUE")) {
      if (nrs->pSolver->allNeumann) {
        nrs->p0the = 0.0;
        for (int ext = 0; ext < nrs->nEXT; ++ext) {
          nrs->p0the += nrs->coeffEXT[ext] * nrs->p0th[ext];
        }
      }
    }
  }

  if (platform->options.compareArgs("MOVING MESH", "TRUE")) {
    if (nrs->cht) {
      mesh = nrs->_mesh;
    }
    nrs->extrapolateKernel(mesh->Nlocal,
                           nrs->NVfields,
                           nrs->nEXT,
                           nrs->fieldOffset,
                           nrs->o_coeffEXT,
                           mesh->o_U,
                           mesh->o_Ue);
  }

  if (nrs->Nscalar > 0) {
    auto cds = nrs->cds;
    if (!cds->o_Se.isInitialized()) {
      return;
    }
    nrs->extrapolateKernel(cds->mesh[0]->Nlocal,
                           cds->NSfields,
                           cds->nEXT,
                           cds->fieldOffset[0],
                           cds->o_coeffEXT,
                           cds->o_S,
                           cds->o_Se);
  }
}

static void computeDivUErr(nrs_t *nrs, dfloat &divUErrVolAvg, dfloat &divUErrL2)
{
  mesh_t *mesh = nrs->meshV;

  auto o_divErr = platform->o_memPool.reserve<dfloat>(nrs->fieldOffset);

  nrs->divergenceVolumeKernel(mesh->Nelements, mesh->o_vgeo, mesh->o_D, nrs->fieldOffset, nrs->o_U, o_divErr);

  double flops = 18 * (mesh->Np * mesh->Nq + mesh->Np);
  flops *= static_cast<double>(mesh->Nelements);

  platform->flopCounter->add("divergenceVolumeKernel", flops);

  oogs::startFinish(o_divErr, 1, nrs->fieldOffset, ogsDfloat, ogsAdd, nrs->gsh);
  platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_invLMM, o_divErr);

  platform->linAlg->axpby(mesh->Nlocal, 1.0, nrs->o_div, -1.0, o_divErr);
  divUErrL2 = platform->linAlg->weightedNorm2(mesh->Nlocal, mesh->o_LMM, o_divErr, platform->comm.mpiComm) /
              sqrt(mesh->volume);

  divUErrVolAvg =
      platform->linAlg->innerProd(mesh->Nlocal, mesh->o_LMM, o_divErr, platform->comm.mpiComm) / mesh->volume;

  divUErrVolAvg = std::abs(divUErrVolAvg);
}

namespace timeStepper
{

void advectionFlops(mesh_t *mesh, int Nfields)
{
  const auto cubNq = mesh->cubNq;
  const auto cubNp = mesh->cubNp;
  const auto Nq = mesh->Nq;
  const auto Np = mesh->Np;
  const auto Nelements = mesh->Nelements;
  double flopCount = 0.0; // per elem basis
  if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
    flopCount += 4. * Nq * (cubNp + cubNq * cubNq * Nq + cubNq * Nq * Nq); // interpolation
    flopCount += 6. * cubNp * cubNq;                                       // apply Dcub
    flopCount += 5 * cubNp; // compute advection term on cubature mesh
    flopCount += mesh->Np;  // weight by inv. mass matrix
  } else {
    flopCount += 8 * (Np * Nq + Np);
  }

  flopCount *= Nelements;
  flopCount *= Nfields;

  platform->flopCounter->add("advection", flopCount);
}

void adjustDt(nrs_t *nrs, int tstep)
{
  const double TOLToZero = (sizeof(dfloat) == sizeof(double)) ? 1e-8 : 1e-5;
  auto initialTimeStepProvided = true;
  if (nrs->dt[0] < TOLToZero && tstep == 1) {
    nrs->dt[0] = 1.0; // startup without any initial timestep guess
    initialTimeStepProvided = false;
  }

  dfloat targetCFL;
  platform->options.getArgs("TARGET CFL", targetCFL);

  const auto CFLmax = 1.2 * targetCFL;
  const auto CFLmin = 0.8 * targetCFL;

  const auto CFL = computeCFL(nrs);

  if (!initialTimeStepProvided) {
    if (CFL > TOLToZero) {
      nrs->dt[0] = targetCFL / CFL * nrs->dt[0];
      nrs->unitTimeCFL = CFL / nrs->dt[0];
    } else {
      // estimate from userf
      if (udf.uEqnSource) {
        platform->linAlg->fill(nrs->fieldOffset * nrs->NVfields, 0.0, nrs->o_FU);
        platform->timer.tic("udfUEqnSource", 1);
        double startTime;
        platform->options.getArgs("START TIME", startTime);
        udf.uEqnSource(nrs, startTime, nrs->o_U, nrs->o_FU);
        platform->timer.toc("udfUEqnSource");

        occa::memory o_FUx = nrs->o_FU + 0 * nrs->fieldOffset;
        occa::memory o_FUy = nrs->o_FU + 1 * nrs->fieldOffset;
        occa::memory o_FUz = nrs->o_FU + 2 * nrs->fieldOffset;

        platform->linAlg->abs(nrs->NVfields * nrs->fieldOffset, nrs->o_FU);

        const auto maxFUx = platform->linAlg->max(nrs->meshV->Nlocal, o_FUx, platform->comm.mpiComm);
        const auto maxFUy = platform->linAlg->max(nrs->meshV->Nlocal, o_FUy, platform->comm.mpiComm);
        const auto maxFUz = platform->linAlg->max(nrs->meshV->Nlocal, o_FUz, platform->comm.mpiComm);
        const auto maxFU = std::max({maxFUx, maxFUy, maxFUz});

        const auto minRho = platform->linAlg->min(nrs->meshV->Nlocal, nrs->o_rho, platform->comm.mpiComm);
        const auto maxU = maxFU / minRho;
        const auto x = nrs->meshV->x;
        const auto y = nrs->meshV->y;
        const auto z = nrs->meshV->z;

        auto minLengthScale = 10 * std::numeric_limits<double>::max();
        for (int i = 0; nrs->meshV->Nlocal; i++) {
          const double lengthScale = sqrt((x[0] - x[1]) * (x[0] - x[1]) + (y[0] - y[1]) * (y[0] - y[1]) +
                                          (z[0] - z[1]) * (z[0] - z[1]));
          minLengthScale = std::min(lengthScale, minLengthScale);
        }

        MPI_Allreduce(MPI_IN_PLACE, &minLengthScale, 1, MPI_DOUBLE, MPI_MIN, platform->comm.mpiComm);
        if (maxU > TOLToZero) {
          nrs->dt[0] = sqrt(targetCFL * minLengthScale / maxU);
        } else {
          nrsAbort(platform->comm.mpiComm,
                   EXIT_FAILURE,
                   "%s\n",
                   "Zero velocity and body force! Please specify an initial timestep!");
        }
      }
    }
    nrs->CFL = CFL;
    return;
  }

  const auto unitTimeCFLold = (tstep == 1) ? CFL / nrs->dt[0] : nrs->unitTimeCFL;
  const auto CFLold = (tstep == 1) ? CFL : nrs->CFL;

  nrs->CFL = CFL;
  nrs->unitTimeCFL = CFL / nrs->dt[0];

  const auto CFLpred = 2.0 * nrs->CFL - CFLold;

  const dfloat TOL = 1e-3;

  if (nrs->CFL > CFLmax || CFLpred > CFLmax || nrs->CFL < CFLmin) {
    const double A = (nrs->unitTimeCFL - unitTimeCFLold) / nrs->dt[0];
    const double B = nrs->unitTimeCFL;
    const double C = -targetCFL;
    const double descriminant = B * B - 4 * A * C;
    nrs->dt[1] = nrs->dt[0];
    if (descriminant <= 0.0) {
      nrs->dt[0] = nrs->dt[0] * (targetCFL / nrs->CFL);
    } else if (std::abs((nrs->unitTimeCFL - unitTimeCFLold) / nrs->unitTimeCFL) < TOL) {
      nrs->dt[0] = nrs->dt[0] * (targetCFL / nrs->CFL);
    } else {
      const double dtLow = (-B + sqrt(descriminant)) / (2.0 * A);
      const double dtHigh = (-B - sqrt(descriminant)) / (2.0 * A);
      if (dtHigh > 0.0 && dtLow > 0.0) {
        nrs->dt[0] = std::min(dtLow, dtHigh);
      } else if (dtHigh <= 0.0 && dtLow <= 0.0) {
        nrs->dt[0] = nrs->dt[0] * targetCFL / nrs->CFL;
      } else {
        nrs->dt[0] = std::max(dtHigh, dtLow);
      }
    }

    // limit dt change
    dfloat maxAdjustDtRatio = 1;
    dfloat minAdjustDtRatio = 1;
    platform->options.getArgs("MAX ADJUST DT RATIO", maxAdjustDtRatio);
    platform->options.getArgs("MIN ADJUST DT RATIO", minAdjustDtRatio);
    if (tstep > 1) {
      nrs->dt[0] = std::max(nrs->dt[0], minAdjustDtRatio * nrs->dt[1]);
    }
    if (tstep > 1) {
      nrs->dt[0] = std::min(nrs->dt[0], maxAdjustDtRatio * nrs->dt[1]);
    }
  }

  nrs->dt[0] = setPrecision(nrs->dt[0], 4); // to avoid accumulation of roundoff errors
}

void initStep(nrs_t *nrs, double time, dfloat dt, int tstep)
{
  nrs->timePrevious = time;

  nrs->tstep = tstep;

  cds_t *cds = nrs->cds;

  const bool movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

  setDt(nrs, dt, tstep);

  extrapolate(nrs);

  if (nrs->Nsubsteps) {
    mesh_t *mesh = nrs->meshV;
    if (nrs->cht) {
      mesh = nrs->cds->mesh[0];
    }
    const auto NCubature = nrs->NVfields * nrs->cubatureOffset;
    for (int s = nrs->nEXT; s > 1; s--) {
      const auto N = nrs->fieldOffset;
      if (movingMesh) {
        mesh->o_divU.copyFrom(mesh->o_divU, N, (s - 1) * N, (s - 2) * N);
        nrs->o_relUrst.copyFrom(nrs->o_relUrst, NCubature, (s - 1) * NCubature, (s - 2) * NCubature);
      } else {
        nrs->o_Urst.copyFrom(nrs->o_Urst, NCubature, (s - 1) * NCubature, (s - 2) * NCubature);
      }
    }
    if (movingMesh) {
      double flops = 18 * (mesh->Np * mesh->Nq + mesh->Np);
      flops *= static_cast<double>(mesh->Nelements);
      nrs->divergenceVolumeKernel(mesh->Nelements,
                                  mesh->o_vgeo,
                                  mesh->o_D,
                                  nrs->fieldOffset,
                                  mesh->o_U,
                                  mesh->o_divU);
      platform->flopCounter->add("divergenceVolumeKernel", flops);
    }
  }

  computeUrst(nrs, movingMesh && nrs->Nsubsteps, platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"));

  if (nrs->Nscalar) {
    if (cds->anyEllipticSolver) {
      platform->timer.tic("makeq", 1);
      platform->linAlg->fill(cds->fieldOffsetSum, 0.0, cds->o_FS);
      makeq(nrs, time, tstep, cds->o_FS, cds->o_BF);
      platform->timer.toc("makeq");
    }
  }

  if (nrs->flow) {
    platform->timer.tic("makef", 1);
    platform->linAlg->fill(nrs->fieldOffset * nrs->NVfields, 0.0, nrs->o_FU);
    makef(nrs, time, tstep, nrs->o_FU, nrs->o_BF);
    platform->timer.toc("makef");
  }

  if (movingMesh) {
    mesh_t *mesh = nrs->_mesh;
    for (int s = std::max(nrs->nBDF, nrs->nEXT); s > 1; s--) {
      const auto N = nrs->fieldOffset;
      mesh->o_LMM.copyFrom(mesh->o_LMM, N, (s - 1) * N, (s - 2) * N);
      mesh->o_invLMM.copyFrom(mesh->o_invLMM, N, (s - 1) * N, (s - 2) * N);
    }

    mesh->move();

    if (nrs->flow) {
      if (bcMap::unalignedMixedBoundary("velocity")) {
        createZeroNormalMask(nrs,
                             nrs->meshV,
                             nrs->uvwSolver->o_EToB,
                             nrs->o_EToBVVelocity,
                             nrs->o_zeroNormalMaskVelocity);
      }
    }

    if (!platform->options.compareArgs("MESH SOLVER", "NONE")) {
      if (bcMap::unalignedMixedBoundary("mesh")) {
        createZeroNormalMask(nrs,
                             mesh,
                             nrs->meshSolver->o_EToB,
                             nrs->o_EToBVMeshVelocity,
                             nrs->o_zeroNormalMaskMeshVelocity);
      }
    }

    if (nrs->cht) {
      nrs->meshV->computeInvLMM();
    }
  }
}

void finishStep(nrs_t *nrs)
{
  nrs->dt[2] = nrs->dt[1];
  nrs->dt[1] = nrs->dt[0];

  if (nrs->flow) {
    if (nrs->pSolver->allNeumann && platform->options.compareArgs("LOWMACH", "TRUE")) {
      nrs->p0the = nrs->p0th[0];
    }
  }

  platform->device.finish();
}

bool runStep(nrs_t *nrs, std::function<bool(int)> convergenceCheck, int stage)
{
  mesh_t *mesh = nrs->meshV;
  cds_t *cds = nrs->cds;

  const auto tstep = nrs->tstep;
  const double timeNew =
      nrs->timePrevious + setPrecision(nrs->dt[0], std::numeric_limits<decltype(nrs->dt[0])>::digits10 + 1);

  const int isOutputStep = nrs->isOutputStep;

  if (nrs->neknek) {
    nrs->neknek->updateBoundary(tstep, stage);
  }

  if (nrs->cvode) {
    scalarSolveCvode(nrs, nrs->timePrevious, timeNew, cds->o_S, stage, tstep);
  }

  if (stage == 1) {
    lagFields(nrs);
  }

  applyDirichlet(nrs, timeNew); // actual + extrapolated solution

  if (nrs->Nscalar) {
    scalarSolve(nrs, timeNew, cds->o_S, stage);
  }

  if (udf.postScalar) {
    udf.postScalar(nrs, timeNew, tstep);
  }

  evaluateProperties(nrs, timeNew);

  if (udf.div) {
    platform->timer.tic("udfDiv", 1);
    platform->linAlg->fill(mesh->Nlocal, 0.0, nrs->o_div);
    udf.div(nrs, timeNew, nrs->o_div);
    platform->timer.toc("udfDiv");
  }

  if (udf.preFluid) {
    udf.preFluid(nrs, timeNew, tstep);
  }

  if (nrs->flow) {
    fluidSolve(nrs, timeNew, nrs->o_P, nrs->o_U, stage, tstep);
  }

  if (!platform->options.compareArgs("MESH SOLVER", "NONE")) {
    occa::memory o_Unew = meshSolve(nrs, timeNew, stage);
    nrs->meshV->o_U.copyFrom(o_Unew, nrs->NVfields * nrs->fieldOffset);
  }

  nrs->timeStepConverged = convergenceCheck(stage);

  platform->timer.tic("udfExecuteStep", 1);
  nek::ifoutfld(0);
  nrs->isOutputStep = 0;
  if (isOutputStep && nrs->timeStepConverged) {
    nek::ifoutfld(1);
    nrs->isOutputStep = 1;
  }
  if (udf.executeStep) {
    udf.executeStep(nrs, timeNew, tstep);
  }
  platform->timer.toc("udfExecuteStep");

  if (!nrs->timeStepConverged) {
    printInfo(nrs, timeNew, tstep, false, true);
  }

  return nrs->timeStepConverged;
}

void setDt(nrs_t *nrs, dfloat dt, int tstep)
{
  nrs->dt[0] = dt;
  nrsCheck(std::isnan(nrs->dt[0]) || std::isinf(nrs->dt[0]),
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "%s",
           "Unreasonable dt! Dying ...\n");

  nrs->idt = 1 / nrs->dt[0];

  const int bdfOrder = std::min(tstep, nrs->nBDF);
  const int extOrder = std::min(tstep, nrs->nEXT);

  nek::bdfCoeff(&nrs->g0, nrs->coeffBDF, nrs->dt, bdfOrder);
  nek::extCoeff(nrs->coeffEXT, nrs->dt, extOrder, bdfOrder);

  for (int i = nrs->nEXT; i > extOrder; i--) {
    nrs->coeffEXT[i - 1] = 0.0;
  }
  for (int i = nrs->nBDF; i > bdfOrder; i--) {
    nrs->coeffBDF[i - 1] = 0.0;
  }

  if (platform->options.compareArgs("MOVING MESH", "TRUE")) {
    mesh_t *mesh = nrs->meshV;
    if (nrs->cht) {
      mesh = nrs->cds->mesh[0];
    }
    const int meshOrder = std::min(tstep, mesh->nAB);
    nek::coeffAB(mesh->coeffAB, nrs->dt, meshOrder);
    for (int i = 0; i < meshOrder; ++i) {
      mesh->coeffAB[i] *= nrs->dt[0];
    }
    for (int i = mesh->nAB; i > meshOrder; i--) {
      mesh->coeffAB[i - 1] = 0.0;
    }
    mesh->o_coeffAB.copyFrom(mesh->coeffAB, mesh->nAB);
  }

  nrs->ig0 = 1.0 / nrs->g0;
  nrs->o_coeffBDF.copyFrom(nrs->coeffBDF);
  nrs->o_coeffEXT.copyFrom(nrs->coeffEXT);

  if (nrs->Nscalar) {
    nrs->cds->idt = 1 / nrs->cds->dt[0];
    nrs->cds->g0 = nrs->g0;
    nrs->cds->ig0 = nrs->ig0;
  }
}

void makeq(nrs_t *nrs, double time, int tstep, occa::memory o_FS, occa::memory o_BF)
{
  cds_t *cds = nrs->cds;

  if (udf.sEqnSource) {
    platform->timer.tic("udfSEqnSource", 1);
    udf.sEqnSource(nrs, time, cds->o_S, o_FS);
    platform->timer.toc("udfSEqnSource");
  }

  for (int is = 0; is < cds->NSfields; is++) {
    if (!cds->compute[is] || cds->cvodeSolve[is]) {
      continue;
    }

    std::string sid = scalarDigitStr(is);

    mesh_t *mesh;
    (is) ? mesh = cds->meshV : mesh = cds->mesh[0];
    const dlong isOffset = cds->fieldOffsetScan[is];

    if (platform->options.compareArgs("SCALAR" + sid + " REGULARIZATION METHOD", "HPFRT")) {
      cds->filterRTKernel(cds->meshV->Nelements,
                          is,
                          1,
                          nrs->fieldOffset,
                          cds->o_applyFilterRT,
                          cds->o_filterRT,
                          cds->o_filterS,
                          cds->o_rho,
                          cds->o_S,
                          o_FS);

      double flops = 6 * mesh->Np * mesh->Nq + 4 * mesh->Np;
      flops *= static_cast<double>(mesh->Nelements);
      platform->flopCounter->add("scalarFilterRT", flops);
    }

    const int movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
    if (movingMesh && !cds->Nsubsteps) {
      cds->advectMeshVelocityKernel(cds->meshV->Nelements,
                                    mesh->o_vgeo,
                                    mesh->o_D,
                                    isOffset,
                                    cds->vFieldOffset,
                                    cds->o_rho,
                                    mesh->o_U,
                                    cds->o_S,
                                    o_FS);
      double flops = 18 * mesh->Np * mesh->Nq + 21 * mesh->Np;
      flops *= static_cast<double>(mesh->Nelements);
      platform->flopCounter->add("scalar advectMeshVelocity", flops);
    }

    occa::memory o_Usubcycling;
    if (platform->options.compareArgs("ADVECTION", "TRUE")) {
      if (cds->Nsubsteps) {
        o_Usubcycling = scalarSubCycle(nrs, std::min(tstep, cds->nEXT), time, is);
      } else {
        if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
          cds->strongAdvectionCubatureVolumeKernel(cds->meshV->Nelements,
                                                   1,
                                                   0, /* weighted */
                                                   0, /* sharedRho */
                                                   mesh->o_vgeo,
                                                   mesh->o_cubDiffInterpT,
                                                   mesh->o_cubInterpT,
                                                   mesh->o_cubProjectT,
                                                   cds->o_compute + is,
                                                   cds->o_fieldOffsetScan + is,
                                                   cds->vFieldOffset,
                                                   cds->vCubatureOffset,
                                                   cds->o_S,
                                                   cds->o_Urst,
                                                   cds->o_rho,
                                                   o_FS);
        } else {
          cds->strongAdvectionVolumeKernel(cds->meshV->Nelements,
                                           1,
                                           0, /* weighted */
                                           mesh->o_vgeo,
                                           mesh->o_D,
                                           cds->o_compute + is,
                                           cds->o_fieldOffsetScan + is,
                                           cds->vFieldOffset,
                                           cds->o_S,
                                           cds->o_Urst,
                                           cds->o_rho,
                                           o_FS);
        }
        advectionFlops(cds->mesh[0], 1);
      }
    }

    cds->sumMakefKernel(mesh->Nlocal,
                        mesh->o_LMM,
                        cds->idt,
                        cds->o_coeffEXT,
                        cds->o_coeffBDF,
                        cds->fieldOffsetSum,
                        cds->fieldOffset[is],
                        isOffset,
                        cds->o_S,
                        o_Usubcycling,
                        o_FS,
                        cds->o_rho,
                        o_BF);

    if (platform->verbose) {
      const dfloat debugNorm = platform->linAlg->weightedNorm2Many(mesh->Nlocal,
                                                                   nrs->NVfields,
                                                                   nrs->fieldOffset,
                                                                   mesh->ogs->o_invDegree,
                                                                   o_BF,
                                                                   platform->comm.mpiComm);
      if (platform->comm.mpiRank == 0) {
        printf("BF scalar=%d norm: %.15e\n", is, debugNorm);
      }
    }

    dfloat scalarSumMakef = (3 * cds->nEXT + 3) * static_cast<double>(mesh->Nlocal);
    scalarSumMakef += (cds->Nsubsteps) ? mesh->Nlocal : 3 * cds->nBDF * static_cast<double>(mesh->Nlocal);
    platform->flopCounter->add("scalarSumMakef", scalarSumMakef);
  }

  for (int s = std::max(cds->nBDF, cds->nEXT); s > 1; s--) {
    const auto N = cds->fieldOffsetSum;
    o_FS.copyFrom(o_FS, N, (s - 1) * N, (s - 2) * N);
  }
}

void scalarSolve(nrs_t *nrs, double time, occa::memory o_S, int stage)
{
  cds_t *cds = nrs->cds;

  // avoid invoking the scalarSolve tic in the case there are no elliptic scalars
  if (!cds->anyEllipticSolver) {
    return;
  }

  platform->timer.tic("scalarSolve", 1);
  for (int is = 0; is < cds->NSfields; is++) {
    if (!cds->compute[is] || cds->cvodeSolve[is]) {
      continue;
    }

    mesh_t *mesh;
    (is) ? mesh = cds->meshV : mesh = cds->mesh[0];

    cds->setEllipticCoeffKernel(mesh->Nlocal,
                                cds->g0 * cds->idt,
                                cds->fieldOffsetScan[is],
                                nrs->fieldOffset,
                                static_cast<int>(cds->o_BFDiag.isInitialized()),
                                cds->o_diff,
                                cds->o_rho,
                                cds->o_BFDiag,
                                cds->o_ellipticCoeff);

    occa::memory o_Snew = cdsSolve(is, cds, time, stage);
    o_Snew.copyTo(o_S, cds->fieldOffset[is], cds->fieldOffsetScan[is]);
  }
  platform->timer.toc("scalarSolve");
}

void scalarSolveCvode(nrs_t *nrs, double tn, double time, occa::memory o_S, int stage, int tstep)
{
  // CVODE is not supported with sub-stepping
  if (!nrs->cvode || stage > 1) {
    return;
  }

  nrs->cvode->solve(tn, time, tstep);
}

void makef(nrs_t *nrs, double time, int tstep, occa::memory o_FU, occa::memory o_BF)
{
  mesh_t *mesh = nrs->meshV;
  const int verbose = platform->options.compareArgs("VERBOSE", "TRUE");
  const int movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

  if (udf.uEqnSource) {
    platform->timer.tic("udfUEqnSource", 1);
    udf.uEqnSource(nrs, time, nrs->o_U, o_FU);
    platform->timer.toc("udfUEqnSource");
  }

  if (platform->options.compareArgs("VELOCITY REGULARIZATION METHOD", "HPFRT")) {
    nrs->filterRTKernel(mesh->Nelements, nrs->o_filterRT, nrs->filterS, nrs->fieldOffset, nrs->o_U, o_FU);
    double flops = 24 * mesh->Np * mesh->Nq + 3 * mesh->Np;
    flops *= static_cast<double>(mesh->Nelements);
    platform->flopCounter->add("velocityFilterRT", flops);
  }


  if (movingMesh && !nrs->Nsubsteps) {
    nrs->advectMeshVelocityKernel(mesh->Nelements,
                                  mesh->o_vgeo,
                                  mesh->o_D,
                                  nrs->fieldOffset,
                                  mesh->o_U,
                                  nrs->o_U,
                                  o_FU);
    double flops = 54 * mesh->Np * mesh->Nq + 63 * mesh->Np;
    flops *= static_cast<double>(mesh->Nelements);
    platform->flopCounter->add("velocity advectMeshVelocity", flops);
  }

  occa::memory o_Usubcycling;
  if (platform->options.compareArgs("ADVECTION", "TRUE")) {
    if (nrs->Nsubsteps) {
      o_Usubcycling = velocitySubCycle(nrs, std::min(tstep, nrs->nEXT), time);
    } else {
      auto o_adv = platform->o_memPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);

      if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
        nrs->strongAdvectionCubatureVolumeKernel(mesh->Nelements,
                                                 mesh->o_vgeo,
                                                 mesh->o_cubDiffInterpT,
                                                 mesh->o_cubInterpT,
                                                 mesh->o_cubProjectT,
                                                 nrs->fieldOffset,
                                                 nrs->cubatureOffset,
                                                 nrs->o_U,
                                                 nrs->o_Urst,
                                                 o_adv);
      } else {
        nrs->strongAdvectionVolumeKernel(mesh->Nelements,
                                         mesh->o_vgeo,
                                         mesh->o_D,
                                         nrs->fieldOffset,
                                         nrs->o_U,
                                         nrs->o_Urst,
                                         o_adv);
      }

      platform->linAlg->axpby(nrs->NVfields * nrs->fieldOffset, -1.0, o_adv, 1.0, o_FU);

      advectionFlops(nrs->meshV, nrs->NVfields);
    }
  }

  nrs->sumMakefKernel(mesh->Nlocal,
                      mesh->o_LMM,
                      nrs->idt,
                      nrs->o_coeffEXT,
                      nrs->o_coeffBDF,
                      nrs->fieldOffset,
                      nrs->o_U,
                      o_Usubcycling,
                      o_FU,
                      o_BF);


  dfloat sumMakefFlops = 0.0;
  if (nrs->Nsubsteps) {
    sumMakefFlops += (6 + 6 * nrs->nEXT) * static_cast<double>(mesh->Nlocal);
  } else {
    sumMakefFlops += (6 * nrs->nEXT + 12 * nrs->nBDF) * static_cast<double>(mesh->Nlocal);
  }

  platform->flopCounter->add("sumMakef", sumMakefFlops);

  if (verbose) {
    const dfloat debugNorm = platform->linAlg->weightedNorm2Many(mesh->Nlocal,
                                                                 nrs->NVfields,
                                                                 nrs->fieldOffset,
                                                                 mesh->ogs->o_invDegree,
                                                                 o_BF,
                                                                 platform->comm.mpiComm);
    if (platform->comm.mpiRank == 0) {
      printf("BF norm: %.15e\n", debugNorm);
    }
  }

  for (int s = std::max(nrs->nBDF, nrs->nEXT); s > 1; s--) {
    const auto N = nrs->NVfields * nrs->fieldOffset;
    o_FU.copyFrom(o_FU, N, (s - 1) * N, (s - 2) * N);
  }
}

void fluidSolve(nrs_t *nrs, double time, occa::memory o_P, occa::memory o_U, int stage, int tstep)
{
  mesh_t *mesh = nrs->meshV;

  platform->timer.tic("pressureSolve", 1);
  {
    platform->timer.tic("pressureSolve", 1);
    nrs->setEllipticCoeffPressureKernel(mesh->Nlocal, nrs->fieldOffset, nrs->o_rho, nrs->o_ellipticCoeff);
    occa::memory o_Pnew = tombo::pressureSolve(nrs, time, stage);
    o_P.copyFrom(o_Pnew);
  }
  platform->timer.toc("pressureSolve");

  platform->timer.tic("velocitySolve", 1);
  {
    nrs->setEllipticCoeffKernel(mesh->Nlocal,
                                nrs->g0 * nrs->idt,
                                0 * nrs->fieldOffset,
                                nrs->fieldOffset,
                                static_cast<int>(nrs->o_BFDiag.isInitialized()),
                                nrs->o_mue,
                                nrs->o_rho,
                                nrs->o_BFDiag,
                                nrs->o_ellipticCoeff);

    occa::memory o_Unew = tombo::velocitySolve(nrs, time, stage);
    o_U.copyFrom(o_Unew, nrs->NVfields * nrs->fieldOffset);
  }
  platform->timer.toc("velocitySolve");

  if (platform->options.compareArgs("CONSTANT FLOW RATE", "TRUE")) {
    ConstantFlowRate::adjust(nrs, tstep, time);
  }
}

void printInfo(nrs_t *nrs, double time, int tstep, bool printStepInfo, bool printVerboseInfo)
{
  cds_t *cds = nrs->cds;

  const double elapsedStep = platform->timer.query("elapsedStep", "DEVICE:MAX");
  const double elapsedStepSum = platform->timer.query("elapsedStepSum", "DEVICE:MAX");
  bool verboseInfo = platform->options.compareArgs("VERBOSE SOLVER INFO", "TRUE");
  const auto cfl = computeCFL(nrs);
  dfloat divUErrVolAvg, divUErrL2;

  if (verboseInfo) {
    computeDivUErr(nrs, divUErrVolAvg, divUErrL2);
  }
  if (platform->comm.mpiRank == 0) {
    if (verboseInfo && printVerboseInfo) {
      bool cvodePrinted = false;
      for (int is = 0; is < nrs->Nscalar; is++) {
        if (cds->compute[is] && !cds->cvodeSolve[is]) {
          elliptic_t *solver = cds->solver[is];
          if (solver->solutionProjection) {
            const int prevVecs = solver->solutionProjection->getPrevNumVecsProjection();
            if (prevVecs > 0) {
              printf("projS%02d  : resNorm0 %.2e  resNorm %.2e  ratio = %.3e  %d/%d\n",
                     is,
                     solver->res00Norm,
                     solver->res0Norm,
                     solver->res00Norm / solver->res0Norm,
                     prevVecs,
                     solver->solutionProjection->getMaxNumVecsProjection());
            }
          }
          printf("S%02d      : iter %03d  resNorm0 %.2e  "
                 "resNorm %.2e\n",
                 is,
                 solver->Niter,
                 solver->res0Norm,
                 solver->resNorm);
        } else if (cds->cvodeSolve[is] && !cvodePrinted) {
          nrs->cvode->printInfo(true);
          cvodePrinted = true;
        }
      }

      if (nrs->neknek) {
        printf("neknek   : sync %.2e  exchange %.2e\n", nrs->neknek->tSync(), nrs->neknek->tExch());
      }

      if (nrs->flow) {
        elliptic_t *solver = nrs->pSolver;
        if (solver->solutionProjection) {
          const int prevVecs = solver->solutionProjection->getPrevNumVecsProjection();
          if (prevVecs > 0) {
            printf("projP    : resNorm0 %.2e  resNorm %.2e  ratio = %.3e  %d/%d\n",
                   solver->res00Norm,
                   solver->res0Norm,
                   solver->res00Norm / solver->res0Norm,
                   prevVecs,
                   solver->solutionProjection->getMaxNumVecsProjection());
          }
        }
        printf("P        : iter %03d  resNorm0 %.2e  resNorm %.2e\n",
               solver->Niter,
               solver->res0Norm,
               solver->resNorm);

        if (nrs->uvwSolver) {
          solver = nrs->uvwSolver;
          if (solver->solutionProjection) {
            const int prevVecs = solver->solutionProjection->getPrevNumVecsProjection();
            if (prevVecs > 0) {
              printf("projUVW  : resNorm0 %.2e  resNorm %.2e  ratio = %.3e  %d/%d\n",
                     solver->res00Norm,
                     solver->res0Norm,
                     solver->res00Norm / solver->res0Norm,
                     prevVecs,
                     solver->solutionProjection->getMaxNumVecsProjection());
            }
          }
          printf("UVW      : iter %03d  resNorm0 %.2e  "
                 "resNorm %.2e  divErrNorms %.2e %.2e\n",
                 solver->Niter,
                 solver->res0Norm,
                 solver->resNorm,
                 divUErrVolAvg,
                 divUErrL2);
        } else {
          solver = nrs->uSolver;
          if (solver->solutionProjection) {
            const int prevVecs = solver->solutionProjection->getPrevNumVecsProjection();
            if (prevVecs > 0) {
              printf("projU    : resNorm0 %.2e  resNorm %.2e  ratio = %.3e  %d/%d\n",
                     solver->res00Norm,
                     solver->res0Norm,
                     solver->res00Norm / solver->res0Norm,
                     prevVecs,
                     solver->solutionProjection->getMaxNumVecsProjection());
            }
          }
          printf("U        : iter %03d  resNorm0 %.2e  "
                 "resNorm %.2e  divErrNorms %.2e %.2e\n",
                 solver->Niter,
                 solver->res0Norm,
                 solver->resNorm,
                 divUErrVolAvg,
                 divUErrL2);
          solver = nrs->vSolver;
          if (solver->solutionProjection) {
            const int prevVecs = solver->solutionProjection->getPrevNumVecsProjection();
            if (prevVecs > 0) {
              printf("projV    : resNorm0 %.2e  resNorm %.2e  ratio = %.3e  %d/%d\n",
                     solver->res00Norm,
                     solver->res0Norm,
                     solver->res00Norm / solver->res0Norm,
                     prevVecs,
                     solver->solutionProjection->getMaxNumVecsProjection());
            }
          }
          printf("V        : iter %03d  resNorm0 %.2e  "
                 "resNorm %.2e\n",
                 solver->Niter,
                 solver->res0Norm,
                 solver->resNorm);
          solver = nrs->wSolver;
          if (solver->solutionProjection) {
            const int prevVecs = solver->solutionProjection->getPrevNumVecsProjection();
            if (prevVecs > 0) {
              printf("projW    : resNorm0 %.2e  resNorm %.2e  ratio = %.3e  %d/%d\n",
                     solver->res00Norm,
                     solver->res0Norm,
                     solver->res00Norm / solver->res0Norm,
                     prevVecs,
                     solver->solutionProjection->getMaxNumVecsProjection());
            }
          }
          printf("W        : iter %03d  resNorm0 %.2e  "
                 "resNorm %.2e\n",
                 solver->Niter,
                 solver->res0Norm,
                 solver->resNorm);
        }
      }

      if (nrs->meshSolver) {
        elliptic_t *solver = nrs->meshSolver;
        if (solver->solutionProjection) {
          const int prevVecs = solver->solutionProjection->getPrevNumVecsProjection();
          if (prevVecs > 0) {
            printf("projMSH  : resNorm0 %.2e  resNorm %.2e  ratio = %.3e  %d/%d\n",
                   solver->res00Norm,
                   solver->res0Norm,
                   solver->res00Norm / solver->res0Norm,
                   prevVecs,
                   solver->solutionProjection->getMaxNumVecsProjection());
          }
        }
        printf("MSH      : iter %03d  resNorm0 %.2e  resNorm %.2e\n",
               solver->Niter,
               solver->res0Norm,
               solver->resNorm);
      }
    }

    if (platform->options.compareArgs("CONSTANT FLOW RATE", "TRUE")) {
      ConstantFlowRate::printInfo(nrs->meshV, verboseInfo && printVerboseInfo);
    }

    if (printStepInfo) {
      printf("step= %d  t= %.8e  dt=%.1e  C= %.2f", tstep, time, nrs->dt[0], cfl);
    }

    if (!verboseInfo) {
      bool cvodePrinted = false;
      for (int is = 0; is < nrs->Nscalar; is++) {
        if (cds->compute[is] && !cds->cvodeSolve[is]) {
          printf("  S: %d", cds->solver[is]->Niter);
        } else if (cds->cvodeSolve[is] && !cvodePrinted) {
          nrs->cvode->printInfo(false);
          cvodePrinted = true;
        }
      }

      if (nrs->flow) {
        printf("  P: %d", nrs->pSolver->Niter);
        if (nrs->uvwSolver) {
          printf("  UVW: %d", nrs->uvwSolver->Niter);
        } else {
          printf("  U: %d  V: %d  W: %d", nrs->uSolver->Niter, nrs->vSolver->Niter, nrs->wSolver->Niter);
        }
      }
      if (nrs->meshSolver) {
        printf("  MSH: %d", nrs->meshSolver->Niter);
      }
    }

    if (nrs->timeStepConverged && printStepInfo) {
      printf("  elapsedStep= %.2es  elapsedStepSum= %.5es\n", elapsedStep, elapsedStepSum);
    }
  }

  if (nrs->cvode)
    nrs->cvode->resetCounters();

  bool largeCFLCheck = (cfl > 30) && numberActiveFields(nrs);

  nrsCheck(largeCFLCheck || std::isnan(cfl) || std::isinf(cfl),
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "%s\n",
           "Unreasonable CFL!");
}

} // namespace timeStepper

#include "nrs.hpp"
#include "avm.hpp"
#include "nekInterfaceAdapter.hpp"
#include "tombo.hpp"
#include "udf.hpp"
#include "bcMap.hpp"
#include "applyDirichlet.hpp"
#include "bdry.hpp"

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
    auto mesh = (nrs->cht) ? nrs->cds->mesh[0] : nrs->mesh;
    for (int s = std::max(nrs->nEXT, mesh->nAB); s > 1; s--) {
      const auto N = nrs->NVfields * nrs->fieldOffset;
      mesh->o_U.copyFrom(mesh->o_U, N, (s - 1) * N, (s - 2) * N);
    }
  }
}

static void extrapolate(nrs_t *nrs)
{
  auto mesh = nrs->mesh;

  nrs->extrapolateKernel(mesh->Nlocal,
                         nrs->NVfields,
                         nrs->nEXT,
                         nrs->fieldOffset,
                         nrs->o_coeffEXT,
                         nrs->o_U,
                         nrs->o_Ue);

  if (platform->options.compareArgs("MOVING MESH", "TRUE")) {
    auto mesh = (nrs->cht) ? nrs->cds->mesh[0] : nrs->mesh;
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

static void setEXTBDFCoeffs(nrs_t *nrs, int tstep)
{
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
    auto mesh = nrs->mesh;
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

  nrs->o_coeffBDF.copyFrom(nrs->coeffBDF);
  nrs->o_coeffEXT.copyFrom(nrs->coeffEXT);
}

static void computeUrst(nrs_t *nrs)
{
  const bool movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

  const bool relative = movingMesh && nrs->Nsubsteps;
  const bool cubature = platform->options.compareArgs("ADVECTION TYPE", "CUBATURE");

  occa::memory &o_Urst = relative ? nrs->o_relUrst : nrs->o_Urst;
  auto mesh = nrs->mesh;
  double flopCount = 0.0;

  if (cubature) {
    nrs->UrstCubatureKernel(mesh->Nelements,
                            static_cast<int>(relative),
                            mesh->o_cubvgeo,
                            mesh->o_cubInterpT,
                            nrs->fieldOffset,
                            nrs->cubatureOffset,
                            nrs->o_U,
                            mesh->o_U,
                            o_Urst);
    flopCount += 6 * mesh->Np * mesh->cubNq;
    flopCount += 6 * mesh->Nq * mesh->Nq * mesh->cubNq * mesh->cubNq;
    flopCount += 6 * mesh->Nq * mesh->cubNp;
    flopCount += 24 * mesh->cubNp;
    flopCount *= mesh->Nelements;
  } else {
    nrs->UrstKernel(mesh->Nelements,
                    static_cast<int>(relative),
                    mesh->o_vgeo,
                    nrs->fieldOffset,
                    nrs->o_U,
                    mesh->o_U,
                    o_Urst);
    flopCount += 24 * static_cast<double>(mesh->Nlocal);
  }
  platform->flopCounter->add("Urst", flopCount);
}

dfloat nrs_t::adjustDt(int tstep)
{
  static auto firstTime = true;
  static dfloat unitTimeCFL;
  static dfloat CFL;

  const double TOLToZero = 1e-6;

  double dt_ = -1;

  dfloat targetCFL;
  platform->options.getArgs("TARGET CFL", targetCFL);

  if (tstep == 1) {
    unitTimeCFL = this->computeCFL(1.0);
    CFL = unitTimeCFL;

    if (unitTimeCFL > TOLToZero) {
      dt_ = targetCFL / unitTimeCFL;
    } else {
      // estimate from userf
      if (this->userVelocitySource) {
        platform->linAlg->fill(this->fieldOffset * this->NVfields, 0.0, this->o_NLT);
        platform->timer.tic("udfUEqnSource", 1);
        double startTime;
        platform->options.getArgs("START TIME", startTime);
        this->userVelocitySource(startTime);
        platform->timer.toc("udfUEqnSource");

        occa::memory o_FUx = this->o_NLT + 0 * this->fieldOffset;
        occa::memory o_FUy = this->o_NLT + 1 * this->fieldOffset;
        occa::memory o_FUz = this->o_NLT + 2 * this->fieldOffset;

        platform->linAlg->abs(this->NVfields * this->fieldOffset, this->o_NLT);

        const auto maxFUx = platform->linAlg->max(this->mesh->Nlocal, o_FUx, platform->comm.mpiComm);
        const auto maxFUy = platform->linAlg->max(this->mesh->Nlocal, o_FUy, platform->comm.mpiComm);
        const auto maxFUz = platform->linAlg->max(this->mesh->Nlocal, o_FUz, platform->comm.mpiComm);
        const auto maxFU = std::max({maxFUx, maxFUy, maxFUz});

        const auto minRho = platform->linAlg->min(this->mesh->Nlocal, this->o_rho, platform->comm.mpiComm);
        const auto maxU = maxFU / minRho;

        std::vector<dfloat> Jw(mesh->Nlocal);
        mesh->o_Jw.copyTo(Jw.data(), Jw.size());

        auto minLengthScale = 10 * std::numeric_limits<dfloat>::max();
        for (int i = 0; i < this->mesh->Nlocal; i++) {
          minLengthScale = std::min(std::cbrt(Jw[i]), minLengthScale);
        }
        MPI_Allreduce(MPI_IN_PLACE, &minLengthScale, 1, MPI_DFLOAT, MPI_MIN, platform->comm.mpiComm);
        if (maxU > TOLToZero) {
          dt_ = sqrt(targetCFL * minLengthScale / maxU);
        } else {
          nekrsAbort(platform->comm.mpiComm,
                     EXIT_FAILURE,
                     "%s\n",
                     "Zero velocity and body force! Please specify an initial timestep!");
        }
      }
    }
    firstTime = false;
    return dt_;
  }

  dt_ = this->dt[0];
  const auto dtOld = this->dt[0];

  auto CFLold = CFL;
  CFL = this->computeCFL();
  if (firstTime) {
    CFLold = CFL;
  }

  auto unitTimeCFLold = unitTimeCFL;
  unitTimeCFL = CFL / dtOld;
  if (firstTime) {
    unitTimeCFLold = unitTimeCFL;
  }

  const auto CFLpred = 2.0 * CFL - CFLold;
  const dfloat TOL = 1e-3;
  const auto CFLmax = 1.2 * targetCFL;
  const auto CFLmin = 0.8 * targetCFL;

  if (CFL > CFLmax || CFLpred > CFLmax || CFL < CFLmin) {
    const double A = (unitTimeCFL - unitTimeCFLold) / dtOld;
    const double B = unitTimeCFL;
    const double C = -targetCFL;
    const double descriminant = B * B - 4 * A * C;

    if (descriminant <= 0.0) {
      dt_ = dtOld * (targetCFL / CFL);
    } else if (std::abs((unitTimeCFL - unitTimeCFLold) / unitTimeCFL) < TOL) {
      dt_ = dtOld * (targetCFL / CFL);
    } else {
      const double dtLow = (-B + sqrt(descriminant)) / (2.0 * A);
      const double dtHigh = (-B - sqrt(descriminant)) / (2.0 * A);
      if (dtHigh > 0.0 && dtLow > 0.0) {
        dt_ = std::min(dtLow, dtHigh);
      } else if (dtHigh <= 0.0 && dtLow <= 0.0) {
        dt_ = dtOld * targetCFL / CFL;
      } else {
        dt_ = std::max(dtHigh, dtLow);
      }
    }
  }
  firstTime = false;

  if (platform->verbose && platform->comm.mpiRank == 0) {
    printf("adjustDt: dt=%g CFL= %g CFLpred= %g CFLmax= %g CFLmin= %g\n", dt_, CFL, CFLpred, CFLmax, CFLmin);
  }

  return dt_;
}

static occa::memory meshSolve(nrs_t *nrs, double time, int iter)
{
  auto mesh = (nrs->cht) ? nrs->cds->mesh[0] : nrs->mesh;
  linAlg_t *linAlg = platform->linAlg;

  auto o_rhs = platform->deviceMemoryPool.reserve<dfloat>(mesh->dim * nrs->fieldOffset);
  platform->linAlg->fill(mesh->dim * nrs->fieldOffset, 0, o_rhs);

  platform->timer.tic("meshSolve", 1);

  auto o_lambda0 = nrs->o_meshMue;

  auto o_U = platform->deviceMemoryPool.reserve<dfloat>(mesh->dim * nrs->fieldOffset);
  if (platform->options.compareArgs("MESH INITIAL GUESS", "EXTRAPOLATION") && iter == 1) {
    o_U.copyFrom(mesh->o_Ue);
  } else {
    o_U.copyFrom(mesh->o_U);
  }

  nrs->meshSolver->solve(o_lambda0, o_NULL, o_rhs, o_U);

  platform->timer.toc("meshSolve");

  return o_U;
}

void nrs_t::initStep(double time, dfloat _dt, int tstep)
{
  if (platform->options.compareArgs("NEKNEK MULTIRATE TIMESTEPPER", "TRUE")) {
    initOuterStep(time, _dt, tstep);
  } else {
    initInnerStep(time, _dt, tstep);
  }
}

void nrs_t::initInnerStep(double time, dfloat _dt, int tstep)
{
  this->timePrevious = time;

  this->tstep = tstep;

  cds_t *cds = this->cds;

  const bool movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

  this->dt[0] = _dt;
  nekrsCheck(this->dt[0] <= 0 || std::isnan(this->dt[0]) || std::isinf(this->dt[0]),
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s",
             "Unreasonable dt!\n");

  setEXTBDFCoeffs(this, tstep);

  extrapolate(this);

  if (this->Nsubsteps) {
    auto mesh = (this->cht) ? this->cds->mesh[0] : this->mesh;
    const auto NCubature = this->NVfields * this->cubatureOffset;
    for (int s = this->nEXT; s > 1; s--) {
      const auto N = this->fieldOffset;
      if (movingMesh) {
        mesh->o_divU.copyFrom(mesh->o_divU, N, (s - 1) * N, (s - 2) * N);
        this->o_relUrst.copyFrom(this->o_relUrst, NCubature, (s - 1) * NCubature, (s - 2) * NCubature);
      } else {
        this->o_Urst.copyFrom(this->o_Urst, NCubature, (s - 1) * NCubature, (s - 2) * NCubature);
      }
    }
    if (movingMesh) {
      double flops = 18 * (mesh->Np * mesh->Nq + mesh->Np);
      flops *= static_cast<double>(mesh->Nelements);
      this->divergenceVolumeKernel(mesh->Nelements,
                                   mesh->o_vgeo,
                                   mesh->o_D,
                                   this->fieldOffset,
                                   mesh->o_U,
                                   mesh->o_divU);
      platform->flopCounter->add("divergenceVolumeKernel", flops);
    }
  }

  computeUrst(this);

  if (this->Nscalar) {
    if (cds->anyEllipticSolver) {
      platform->timer.tic("makeq", 1);

      platform->linAlg->fill(cds->fieldOffsetSum, 0.0, cds->o_NLT);

      if (this->userScalarSource) {
        platform->timer.tic("udfSEqnSource", 1);
        this->userScalarSource(time);
        platform->timer.toc("udfSEqnSource");
      }

      for (int is = 0; is < this->Nscalar; is++) {
        if (!cds->compute[is] || cds->cvodeSolve[is]) {
          continue;
        }

        const std::string sid = scalarDigitStr(is);

        auto mesh = cds->mesh[is];
        const dlong isOffset = cds->fieldOffsetScan[is];

        occa::memory o_Usubcycling;
        cds->makeNLT(is, time, tstep, o_Usubcycling);

        this->sumMakefKernel(mesh->Nlocal,
                             0,
                             mesh->o_LMM,
                             1 / this->dt[0],
                             this->o_coeffEXT,
                             this->o_coeffBDF,
                             isOffset,
                             cds->fieldOffsetSum,
                             this->fieldOffset,
                             cds->o_rho,
                             cds->o_S,
                             o_Usubcycling,
                             cds->o_NLT,
                             cds->o_JwF);

        dfloat scalarSumMakef = (3 * this->nEXT + 3);
        scalarSumMakef += (this->Nsubsteps) ? 1 : 3 * this->nBDF;
        platform->flopCounter->add("scalarSumMakef", scalarSumMakef * static_cast<double>(mesh->Nlocal));

        if (platform->verbose) {
          const dfloat debugNorm = platform->linAlg->weightedNorm2Many(mesh->Nlocal,
                                                                       1,
                                                                       isOffset,
                                                                       mesh->ogs->o_invDegree,
                                                                       cds->o_JwF + cds->fieldOffsetScan[is],
                                                                       platform->comm.mpiComm);
          if (platform->comm.mpiRank == 0) {
            printf("BF scalar=%d norm: %.15e\n", is, debugNorm);
          }
        }
      }

      for (int s = std::max(this->nBDF, this->nEXT); s > 1; s--) {
        const auto N = cds->fieldOffsetSum;
        cds->o_NLT.copyFrom(cds->o_NLT, N, (s - 1) * N, (s - 2) * N);
      }

      platform->timer.toc("makeq");
    }
  }

  if (this->flow) {
    platform->timer.tic("makef", 1);

    platform->linAlg->fill(this->fieldOffset * this->NVfields, 0.0, this->o_NLT);

    if (this->userVelocitySource) {
      platform->timer.tic("udfUEqnSource", 1);
      this->userVelocitySource(time);
      platform->timer.toc("udfUEqnSource");
    }

    occa::memory o_Usubcycling;
    this->makeNLT(time, tstep, o_Usubcycling);

    this->sumMakefKernel(mesh->Nlocal,
                         1,
                         mesh->o_LMM,
                         1 / this->dt[0],
                         this->o_coeffEXT,
                         this->o_coeffBDF,
                         this->fieldOffset,
                         this->NVfields * this->fieldOffset,
                         this->fieldOffset,
                         this->o_rho,
                         this->o_U,
                         o_Usubcycling,
                         this->o_NLT,
                         this->o_JwF);

    dfloat sumMakefFlops = (this->Nsubsteps) ? (6 + 6 * this->nEXT) : (6 * this->nEXT + 12 * this->nBDF);
    platform->flopCounter->add("sumMakef", sumMakefFlops * static_cast<double>(mesh->Nlocal));

    if (platform->verbose) {
      const dfloat debugNorm = platform->linAlg->weightedNorm2Many(mesh->Nlocal,
                                                                   this->NVfields,
                                                                   this->fieldOffset,
                                                                   mesh->ogs->o_invDegree,
                                                                   this->o_JwF,
                                                                   platform->comm.mpiComm);
      if (platform->comm.mpiRank == 0) {
        printf("BF velocity norm: %.15e\n", debugNorm);
      }
    }

    for (int s = std::max(this->nBDF, this->nEXT); s > 1; s--) {
      const auto N = this->NVfields * this->fieldOffset;
      this->o_NLT.copyFrom(this->o_NLT, N, (s - 1) * N, (s - 2) * N);
    }

    platform->timer.toc("makef");
  }

  if (movingMesh) {
    auto mesh = (cht) ? cds->mesh[0] : this->mesh;
    for (int s = std::max(this->nBDF, this->nEXT); s > 1; s--) {
      const auto N = this->fieldOffset;
      mesh->o_LMM.copyFrom(mesh->o_LMM, N, (s - 1) * N, (s - 2) * N);
      mesh->o_invLMM.copyFrom(mesh->o_invLMM, N, (s - 1) * N, (s - 2) * N);
    }

    mesh->move();
    if (cht) {
      mesh->computeInvLMM();
    }

    if (this->flow) {
      if (bcMap::unalignedMixedBoundary("velocity")) {
        createZeroNormalMask(this,
                             this->mesh,
                             this->uvwSolver->o_EToB(),
                             this->o_EToBVVelocity,
                             this->o_zeroNormalMaskVelocity);
      }
    }

    if (!platform->options.compareArgs("MESH SOLVER", "NONE")) {
      if (bcMap::unalignedMixedBoundary("mesh")) {
        createZeroNormalMask(this,
                             mesh,
                             this->meshSolver->o_EToB(),
                             this->o_EToBVMeshVelocity,
                             this->o_zeroNormalMaskMeshVelocity);
      }
    }
  }
}

void nrs_t::finishStep()
{
  if (platform->options.compareArgs("NEKNEK MULTIRATE TIMESTEPPER", "TRUE")) {
    finishOuterStep();
  } else {
    finishInnerStep();
  }
}

void nrs_t::finishInnerStep()
{
  this->dt[2] = this->dt[1];
  this->dt[1] = this->dt[0];

  platform->device.finish();
}

bool nrs_t::runStep(std::function<bool(int)> convergenceCheck, int iter)
{
  timeStepConverged = false;

  if (platform->options.compareArgs("NEKNEK MULTIRATE TIMESTEPPER", "TRUE")) {
    runOuterStep(convergenceCheck, iter);
  } else {
    runInnerStep(convergenceCheck, iter, true);
  }

  return timeStepConverged;
}

bool nrs_t::runInnerStep(std::function<bool(int)> convergenceCheck, int iter, bool outerConverged)
{
  this->outerCorrector = iter;

  const auto tstep = this->tstep;
  const double timeNew = this->timePrevious + setPrecision(this->dt[0], 5);

  const auto checkpointStep0 = checkpointStep;

  if (this->neknek) {
    this->neknek->updateBoundary(tstep, iter, timeNew);
  }

  if (cds) {
    if (cds->cvode && iter == 1) {
      cds->cvode->solve(this->timePrevious, timeNew, tstep);
    }
  }

  if (iter == 1) {
    lagFields(this);
  }

  {
    if (this->flow) {
      applyDirichletVelocity(this, timeNew, this->o_U, this->o_Ue, this->o_P);
    }

    if (this->Nscalar) {
      applyDirichletScalars(this, timeNew, this->cds->o_S, this->cds->o_Se);
    }

    if (!platform->options.compareArgs("MESH SOLVER", "NONE")) {
      auto mesh = (cht) ? cds->mesh[0] : this->mesh;
      applyDirichletMesh(this, timeNew, mesh->o_U, mesh->o_Ue, this->o_U);
    }

    if (this->neknek) {
      this->neknek->fixCoupledSurfaceFlux(this->o_U);
    }
  }

  if (this->Nscalar) {
    cds->solve(timeNew, iter);
  }

  if (this->postScalar) {
    this->postScalar(timeNew, tstep);
  }

  this->evaluateProperties(timeNew);

  this->evaluateDivergence(timeNew);

  if (this->preFluid) {
    this->preFluid(timeNew, tstep);
  }

  if (this->flow) {
    platform->timer.tic("pressureSolve", 1);
    {
      const auto &o_Pnew = tombo::pressureSolve(this, timeNew, iter);
      this->o_P.copyFrom(o_Pnew, mesh->Nlocal);
    }
    platform->timer.toc("pressureSolve");

    platform->timer.tic("velocitySolve", 1);
    {
      auto o_Unew = tombo::velocitySolve(this, timeNew, iter);
      this->o_U.copyFrom(o_Unew, this->NVfields * this->fieldOffset);
      o_Unew.free();
    }
    platform->timer.toc("velocitySolve");

    if (platform->options.compareArgs("CONSTANT FLOW RATE", "TRUE")) {
      this->adjustFlowRate(tstep, timeNew);
    }
  }

  if (!platform->options.compareArgs("MESH SOLVER", "NONE")) {
    auto o_Unew = meshSolve(this, timeNew, iter);
    mesh->o_U.copyFrom(o_Unew, this->NVfields * this->fieldOffset);
  }

  const auto converged = convergenceCheck(iter);

  timeStepConverged = outerConverged && converged;

  nek::ifoutfld(0);
  checkpointStep = 0;
  if (checkpointStep0 && timeStepConverged) {
    nek::ifoutfld(1);
    checkpointStep = 1;
  }

  platform->timer.tic("udfExecuteStep", 1);
  if (udf.executeStep) {
    udf.executeStep(timeNew, tstep);
  }
  platform->timer.toc("udfExecuteStep");

  return converged;
}

void nrs_t::saveSolutionState()
{
  const bool movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
  if (!o_Usave.isInitialized()) {
    o_Usave = platform->device.malloc<dfloat>(o_U.length());
    o_Psave = platform->device.malloc<dfloat>(o_P.length());
    o_NLTsave = platform->device.malloc<dfloat>(o_NLT.length());
    o_Upropsave = platform->device.malloc<dfloat>(o_prop.length());
    if (Nsubsteps) {
      auto o_Urst_ = movingMesh ? o_relUrst : o_Urst;
      o_Urstsave = platform->device.malloc<dfloat>(o_Urst_.length());
    }
    if (movingMesh) {
      auto mesh = (cht) ? cds->mesh[0] : this->mesh;
      o_LMMsave = platform->device.malloc<dfloat>(mesh->o_LMM.length());
      o_Umeshsave = platform->device.malloc<dfloat>(mesh->o_U.length());
      o_xsave = platform->device.malloc<dfloat>(mesh->o_x.length());
      o_ysave = platform->device.malloc<dfloat>(mesh->o_y.length());
      o_zsave = platform->device.malloc<dfloat>(mesh->o_z.length());
    }
  }

  o_Usave.copyFrom(o_U, o_U.length());
  o_Psave.copyFrom(o_P, o_P.length());
  o_NLTsave.copyFrom(o_NLT, o_NLT.length());
  o_Upropsave.copyFrom(o_prop, o_prop.length());
  if (Nsubsteps) {
    auto o_Urst_ = movingMesh ? o_relUrst : o_Urst;
    o_Urstsave.copyFrom(o_Urst, o_Urst_.length());
  }
  if (cds) {
    cds->saveSolutionState();
  }

  if (movingMesh) {
    auto mesh = (cht) ? cds->mesh[0] : this->mesh;
    o_LMMsave.copyFrom(mesh->o_LMM, mesh->o_LMM.length());
    o_Umeshsave.copyFrom(mesh->o_U, mesh->o_U.length());
    o_xsave.copyFrom(mesh->o_x, mesh->o_x.length());
    o_ysave.copyFrom(mesh->o_y, mesh->o_y.length());
    o_zsave.copyFrom(mesh->o_z, mesh->o_z.length());
  }
}

void nrs_t::restoreSolutionState()
{
  const bool movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
  o_Usave.copyTo(o_U, o_U.length());
  o_Psave.copyTo(o_P, o_P.length());
  o_NLTsave.copyTo(o_NLT, o_NLT.length());
  o_Upropsave.copyTo(o_prop, o_prop.length());
  if (cds) {
    cds->restoreSolutionState();
  }
  if (Nsubsteps) {
    auto o_Urst_ = movingMesh ? o_relUrst : o_Urst;
    o_Urstsave.copyTo(o_Urst, o_Urst_.length());
  }
  if (movingMesh) {
    auto mesh = (cht) ? cds->mesh[0] : this->mesh;
    o_LMMsave.copyTo(mesh->o_LMM, mesh->o_LMM.length());
    o_Umeshsave.copyTo(mesh->o_U, mesh->o_U.length());
    o_xsave.copyTo(mesh->o_x, mesh->o_x.length());
    o_ysave.copyTo(mesh->o_y, mesh->o_y.length());
    o_zsave.copyTo(mesh->o_z, mesh->o_z.length());
    mesh->update();
  }
}

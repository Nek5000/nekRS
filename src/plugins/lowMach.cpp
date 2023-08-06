#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"

#include "lowMach.hpp"
#include "linAlg.hpp"

namespace {

nrs_t *the_nrs = nullptr;
linAlg_t *the_linAlg = nullptr;

int qThermal = 0;
dfloat alpha0 = 1.0;

occa::memory o_beta;
occa::memory o_kappa;

occa::memory h_scratch;

occa::kernel qtlKernel;
occa::kernel p0thHelperKernel;
occa::kernel surfaceFluxKernel;
}

void lowMach::buildKernel(occa::properties kernelInfo)
{
  int rank = platform->comm.mpiRank;
  const std::string path = getenv("NEKRS_KERNEL_DIR") + std::string("/plugins/");
  std::string kernelName, fileName;
  const std::string extension = ".okl";
  {
    kernelName = "qtlHex3D";
    fileName = path + kernelName + extension;
    qtlKernel = platform->device.buildKernel(fileName, kernelInfo, true);

    kernelName = "p0thHelper";
    fileName = path + kernelName + extension;
    p0thHelperKernel = platform->device.buildKernel(fileName, kernelInfo, true);

    kernelName = "surfaceFlux";
    fileName = path + kernelName + extension;
    surfaceFluxKernel = platform->device.buildKernel(fileName, kernelInfo, true);
  }

  platform->options.setArgs("LOWMACH", "TRUE");
}

void lowMach::setup(nrs_t *nrs, dfloat alpha_, occa::memory& o_beta_, occa::memory& o_kappa_)
{
  the_nrs = nrs;

  alpha0 = alpha_;
  nrs->alpha0Ref = alpha0;
  o_beta = o_beta_;
  o_kappa = o_kappa_;

  the_linAlg = platform->linAlg;
  mesh_t *mesh = nrs->meshV;
  int err = 1;
  if (platform->options.compareArgs("SCALAR00 IS TEMPERATURE", "TRUE"))
    err = 0;

  nrsCheck(err, platform->comm.mpiComm, EXIT_FAILURE,
           "%s\n", "requires solving for temperature!");

  h_scratch = platform->device.mallocHost(mesh->Nelements * sizeof(dfloat));
}

void lowMach::qThermalSingleComponent(dfloat time, occa::memory& o_div)
{
  qThermal = 1;
  nrs_t *nrs = the_nrs;
  cds_t *cds = nrs->cds;
  mesh_t *mesh = nrs->meshV;
  linAlg_t *linAlg = platform->linAlg;

  bool rhsCVODE = false;
  std::string scope = "udfDiv::";
  if(cds->cvode){
    rhsCVODE = cds->cvode->isRhsEvaluation();
    if(rhsCVODE){
      scope = cds->cvode->scope() + "::";
    }
  }

  nrs->gradientVolumeKernel(mesh->Nelements,
                            mesh->o_vgeo,
                            mesh->o_D,
                            nrs->fieldOffset,
                            cds->o_S,
                            platform->o_mempool.slice0);

  double flopsGrad = 6 * mesh->Np * mesh->Nq + 18 * mesh->Np;
  flopsGrad *= static_cast<double>(mesh->Nelements);

  oogs::startFinish(platform->o_mempool.slice0, nrs->NVfields, nrs->fieldOffset, ogsDfloat, ogsAdd, nrs->gsh);

  platform->linAlg
      ->axmyVector(mesh->Nlocal, nrs->fieldOffset, 0, 1.0, nrs->meshV->o_invLMM, platform->o_mempool.slice0);

  platform->linAlg->fill(mesh->Nelements * mesh->Np, 0.0, platform->o_mempool.slice3);
  if (udf.sEqnSource) {
    platform->timer.tic(scope + "udfSEqnSource", 1);
    udf.sEqnSource(nrs, time, cds->o_S, platform->o_mempool.slice3);
    platform->timer.toc(scope + "udfSEqnSource");
  }

  qtlKernel(mesh->Nelements,
            mesh->o_vgeo,
            mesh->o_D,
            nrs->fieldOffset,
            platform->o_mempool.slice0,
            o_beta,
            cds->o_diff,
            cds->o_rho,
            platform->o_mempool.slice3,
            o_div);

  double flopsQTL = 18 * mesh->Np * mesh->Nq + 23 * mesh->Np;
  flopsQTL *= static_cast<double>(mesh->Nelements);

  oogs::startFinish(o_div, 1, nrs->fieldOffset, ogsDfloat, ogsAdd, nrs->gsh);

  platform->linAlg->axmy(mesh->Nlocal, 1.0, nrs->meshV->o_invLMM, o_div);

  double surfaceFlops = 0.0;

  if (nrs->pSolver) {
    const bool closedVolume = nrs->pSolver->allNeumann;
    if(!closedVolume)
     return;

    const dlong Nlocal = mesh->Nlocal;

    linAlg->axmyz(Nlocal, 1.0, mesh->o_LMM, o_div, platform->o_mempool.slice0);
    const dfloat termQ = linAlg->sum(Nlocal, platform->o_mempool.slice0, platform->comm.mpiComm);

    surfaceFluxKernel(mesh->Nelements,
                      mesh->o_sgeo,
                      mesh->o_vmapM,
                      nrs->o_EToB,
                      nrs->fieldOffset,
                      rhsCVODE ? nrs->o_U : nrs->o_Ue,
                      platform->o_mempool.slice0);

    double surfaceFluxFlops = 13 * mesh->Nq * mesh->Nq;
    surfaceFluxFlops *= static_cast<double>(mesh->Nelements);

    platform->o_mempool.slice0.copyTo(h_scratch.ptr(), mesh->Nelements * sizeof(dfloat));
    auto scratch = (dfloat *) h_scratch.ptr();

    dfloat termV = 0.0;
    for (int i = 0; i < mesh->Nelements; ++i) {
      termV += scratch[i];
    }
    MPI_Allreduce(MPI_IN_PLACE, &termV, 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);

    p0thHelperKernel(Nlocal,
                     alpha0,
                     nrs->p0th[0],
                     o_beta,
                     o_kappa,
                     cds->o_rho,
                     nrs->meshV->o_LMM,
                     platform->o_mempool.slice0,
                     platform->o_mempool.slice1);

    double p0thHelperFlops = 4 * mesh->Nlocal;

    const dfloat prhs =
        (termQ - termV) / linAlg->sum(Nlocal, platform->o_mempool.slice0, platform->comm.mpiComm);
    linAlg->axpby(Nlocal, -prhs, platform->o_mempool.slice1, 1.0, o_div);

    const auto *coeff = rhsCVODE ? nrs->cvode->coeffBDF() : nrs->coeffBDF;
    dfloat Saqpq = 0.0;
    for (int i = 0; i < nrs->nBDF; ++i) {
      Saqpq += coeff[i] * nrs->p0th[i];
    }

    const auto g0 = rhsCVODE ? nrs->cvode->g0() : nrs->g0;
    const auto dt = rhsCVODE ? nrs->cvode->dt() : nrs->dt[0];

    const auto pcoef = (g0 - dt * prhs);
    const auto p0thn = Saqpq / pcoef;

    // only update p0th when not inside a CVODE evaluation
    if(!rhsCVODE){
      nrs->p0th[2] = nrs->p0th[1];
      nrs->p0th[1] = nrs->p0th[0];
      nrs->p0th[0] = p0thn;
    }

    nrs->dp0thdt = prhs * p0thn;

    surfaceFlops += surfaceFluxFlops + p0thHelperFlops;
  }

  qThermal = 0;

  double flops = surfaceFlops + flopsGrad + flopsQTL;
  platform->flopCounter->add("lowMach::qThermalRealGasSingleComponent", flops);
}

void lowMach::dpdt(occa::memory& o_FU)
{
  nrs_t *nrs = the_nrs;
  mesh_t *mesh = nrs->meshV;

  if(nrs->cds->cvodeSolve[0]) return; // contribution is not applied here

  if (!qThermal) {
    platform->linAlg->add(mesh->Nlocal, nrs->dp0thdt * alpha0, o_FU);
  }
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"

#include "lowMach.hpp"
#include "linAlg.hpp"

namespace {

nrs_t *_nrs = nullptr;
linAlg_t *the_linAlg = nullptr;

int qThermal = 0;
dfloat alpha0 = 1.0;

occa::memory o_beta;
occa::memory o_kappa;

occa::memory o_bID;

occa::kernel qtlKernel;
occa::kernel p0thHelperKernel;

static bool buildKernelCalled = false;
static bool setupCalled = false;

}

void lowMach::buildKernel(occa::properties kernelInfo)
{
  static bool isInitialized = false;
  if (isInitialized) return;
  isInitialized = true;

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
  }

  platform->options.setArgs("LOWMACH", "TRUE");
  buildKernelCalled = true;
}

void lowMach::setup(nrs_t *nrs, dfloat alpha_, occa::memory& o_beta_, occa::memory& o_kappa_)
{
  static bool isInitialized = false;
  if (isInitialized) return;
  isInitialized = true;

  _nrs = nrs;

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

  std::vector<int> bID;
  for (auto& [key, bcID] :  bcMap::map()) {
    const auto field = key.first;
    if (field == "velocity") {
      if (bcID == bcMap::bcTypeV || bcID == bcMap::bcTypeINT) {
        bID.push_back(key.second + 1);
      }
    }
  }
  o_bID = platform->device.malloc<int>(bID.size());
  o_bID.copyFrom(bID.data());

  setupCalled = true; 
}

void lowMach::qThermalSingleComponent(double time, occa::memory& o_div)
{
  nrsCheck(!setupCalled || !buildKernelCalled, MPI_COMM_SELF, EXIT_FAILURE,
           "%s\n", "called prior to tavg::setup()!");

  qThermal = 1;
  nrs_t *nrs = _nrs;
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

  auto o_gradT = platform->o_memPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);
  nrs->gradientVolumeKernel(mesh->Nelements,
                            mesh->o_vgeo,
                            mesh->o_D,
                            nrs->fieldOffset,
                            cds->o_S,
                            o_gradT);

  double flopsGrad = 6 * mesh->Np * mesh->Nq + 18 * mesh->Np;
  flopsGrad *= static_cast<double>(mesh->Nelements);

  oogs::startFinish(o_gradT, nrs->NVfields, nrs->fieldOffset, ogsDfloat, ogsAdd, nrs->gsh);

  platform->linAlg
      ->axmyVector(mesh->Nlocal, nrs->fieldOffset, 0, 1.0, nrs->meshV->o_invLMM, o_gradT);


  auto o_src = platform->o_memPool.reserve<dfloat>(nrs->fieldOffset);
  platform->linAlg->fill(mesh->Nlocal, 0.0, o_src);
  if (udf.sEqnSource) {
    platform->timer.tic(scope + "udfSEqnSource", 1);
    udf.sEqnSource(nrs, time, cds->o_S, o_src);
    platform->timer.toc(scope + "udfSEqnSource");
  }

  qtlKernel(mesh->Nelements,
            mesh->o_vgeo,
            mesh->o_D,
            nrs->fieldOffset,
            o_gradT,
            o_beta,
            cds->o_diff,
            cds->o_rho,
            o_src,
            o_div);

  o_gradT.free();
  o_src.free();

  double flopsQTL = 18 * mesh->Np * mesh->Nq + 23 * mesh->Np;
  flopsQTL *= static_cast<double>(mesh->Nelements);

  oogs::startFinish(o_div, 1, nrs->fieldOffset, ogsDfloat, ogsAdd, nrs->gsh);

  platform->linAlg->axmy(mesh->Nlocal, 1.0, nrs->meshV->o_invLMM, o_div);

  double surfaceFlops = 0.0;

  if (nrs->pSolver) {
    if(!nrs->pSolver->allNeumann) return;

    const auto termQ = [&]() 
    {
      auto o_tmp = platform->o_memPool.reserve<dfloat>(nrs->fieldOffset);
      linAlg->axmyz(mesh->Nlocal, 1.0, mesh->o_LMM, o_div, o_tmp);
      return linAlg->sum(mesh->Nlocal, o_tmp, platform->comm.mpiComm);
    }();

    auto o_tmp1 = platform->o_memPool.reserve<dfloat>(nrs->fieldOffset);
    auto o_tmp2 = platform->o_memPool.reserve<dfloat>(nrs->fieldOffset);
    p0thHelperKernel(mesh->Nlocal,
                     alpha0,
                     nrs->p0th[0],
                     o_beta,
                     o_kappa,
                     cds->o_rho,
                     nrs->meshV->o_LMM,
                     o_tmp1,
                     o_tmp2);

    double p0thHelperFlops = 4 * mesh->Nlocal;

    const auto flux = mesh->surfaceIntegralVector(nrs->fieldOffset, 
                                                   o_bID.length(), 
                                                   o_bID, 
                                                   rhsCVODE ? nrs->o_U : nrs->o_Ue);
    const auto termV = std::accumulate(flux.begin(), flux.end(), 0.0); 

    double surfaceFluxFlops = 13 * mesh->Nq * mesh->Nq;
    surfaceFluxFlops *= static_cast<double>(mesh->Nelements);

    const auto prhs = (termQ - termV) / linAlg->sum(mesh->Nlocal, o_tmp1, platform->comm.mpiComm);
    linAlg->axpby(mesh->Nlocal, -prhs, o_tmp2, 1.0, o_div);
    o_tmp1.free();
    o_tmp2.free();

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
  nrsCheck(!setupCalled || !buildKernelCalled, MPI_COMM_SELF, EXIT_FAILURE,
           "%s\n", "called prior to tavg::setup()!");

  nrs_t *nrs = _nrs;
  mesh_t *mesh = nrs->meshV;

  if(nrs->cds->cvodeSolve[0]) return; // contribution is not applied here

  if (!qThermal) {
    platform->linAlg->add(mesh->Nlocal, nrs->dp0thdt * alpha0, o_FU);
  }
}

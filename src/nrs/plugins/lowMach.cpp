#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"

#include "lowMach.hpp"
#include "linAlg.hpp"

namespace
{

nrs_t *_nrs = nullptr;

int qThermal = 0;
dfloat alpha0 = 1.0;

occa::memory o_beta;
occa::memory o_kappa;

occa::memory o_bID;

occa::kernel qtlKernel;
occa::kernel p0thHelperKernel;

static bool buildKernelCalled = false;
static bool setupCalled = false;

} // namespace

void lowMach::buildKernel(occa::properties kernelInfo)
{
  auto buildKernel = [&kernelInfo](const std::string &kernelName) {
    const auto path = getenv("NEKRS_KERNEL_DIR") + std::string("/nrs/plugins/");
    const auto fileName = path + "lowMach.okl";
    const auto reqName = "lowMach::";
    if (platform->options.compareArgs("REGISTER ONLY", "TRUE")) {
      platform->kernelRequests.add(reqName, fileName, kernelInfo);
      return occa::kernel();
    } else {
      buildKernelCalled = 1;
      return platform->kernelRequests.load(reqName, kernelName);
    }
  };

  qtlKernel = buildKernel("qtlHex3D");
  p0thHelperKernel = buildKernel("p0thHelper");

  platform->options.setArgs("LOWMACH", "TRUE");
}

void lowMach::setup(dfloat alpha_, const occa::memory &o_beta_, const occa::memory &o_kappa_)
{
  static bool isInitialized = false;
  if (isInitialized) {
    return;
  }
  isInitialized = true;

  _nrs = dynamic_cast<nrs_t *>(platform->solver);
  ;

  alpha0 = alpha_;
  _nrs->alpha0Ref = alpha0;
  o_beta = o_beta_;
  o_kappa = o_kappa_;

  int err = 1;
  if (platform->options.compareArgs("SCALAR00 IS TEMPERATURE", "TRUE")) {
    err = 0;
  }

  nekrsCheck(err, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "requires solving for temperature!");

  std::vector<int> bID;
  for (auto &[key, bcID] : bcMap::map()) {
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

void lowMach::qThermalSingleComponent(double time)
{
  auto o_div = _nrs->o_div;
  nekrsCheck(!setupCalled || !buildKernelCalled,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s\n",
             "called prior to tavg::setup()!");

  qThermal = 1;
  auto nrs = _nrs;
  auto cds = nrs->cds;
  auto mesh = nrs->mesh;
  linAlg_t *linAlg = platform->linAlg;

  std::string scope = "udfDiv::";

  bool rhsCVODE = false;
  if (cds->cvode) {
    rhsCVODE = cds->cvode->isRhsEvaluation();
  }

  auto o_gradT = platform->deviceMemoryPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);
  nrs->gradientVolumeKernel(mesh->Nelements, mesh->o_vgeo, mesh->o_D, nrs->fieldOffset, cds->o_S, o_gradT);

  double flopsGrad = 6 * mesh->Np * mesh->Nq + 18 * mesh->Np;
  flopsGrad *= static_cast<double>(mesh->Nelements);

  oogs::startFinish(o_gradT, nrs->NVfields, nrs->fieldOffset, ogsDfloat, ogsAdd, nrs->gsh);

  platform->linAlg->axmyVector(mesh->Nlocal, nrs->fieldOffset, 0, 1.0, nrs->mesh->o_invLMM, o_gradT);

  auto o_src = platform->deviceMemoryPool.reserve<dfloat>(nrs->fieldOffset);
  platform->linAlg->fill(mesh->Nlocal, 0.0, o_src);
  if (cds->userSource) {
    platform->timer.tic(scope + "udfSEqnSource", 1);
    auto o_saveFS = cds->o_NLT;
    cds->o_NLT = o_src;
    cds->userSource(time);
    cds->o_NLT = o_saveFS;
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

  platform->linAlg->axmy(mesh->Nlocal, 1.0, nrs->mesh->o_invLMM, o_div);

  double surfaceFlops = 0.0;

  if (nrs->pSolver) {
    if (!nrs->pSolver->nullSpace()) {
      return;
    }

    nekrsCheck(rhsCVODE,
               MPI_COMM_SELF,
               EXIT_FAILURE,
               "%s\n",
               "computing p0th and dp0thdt using CVODE is not supported!");

    const auto termQ = [&]() {
      auto o_tmp = platform->deviceMemoryPool.reserve<dfloat>(nrs->fieldOffset);
      linAlg->axmyz(mesh->Nlocal, 1.0, mesh->o_LMM, o_div, o_tmp);
      return linAlg->sum(mesh->Nlocal, o_tmp, platform->comm.mpiComm);
    }();

    auto o_tmp1 = platform->deviceMemoryPool.reserve<dfloat>(nrs->fieldOffset);
    auto o_tmp2 = platform->deviceMemoryPool.reserve<dfloat>(nrs->fieldOffset);
    p0thHelperKernel(mesh->Nlocal,
                     alpha0,
                     nrs->p0th[0],
                     o_beta,
                     o_kappa,
                     cds->o_rho,
                     nrs->mesh->o_LMM,
                     o_tmp1,
                     o_tmp2);

    double p0thHelperFlops = 4 * mesh->Nlocal;

    const auto flux =
        mesh->surfaceAreaNormalMultiplyIntegrate(nrs->fieldOffset, o_bID.length(), o_bID, nrs->o_Ue);
    const auto termV = std::accumulate(flux.begin(), flux.end(), 0.0);

    double surfaceFluxFlops = 13 * mesh->Nq * mesh->Nq;
    surfaceFluxFlops *= static_cast<double>(mesh->Nelements);

    const auto prhs = (termQ - termV) / linAlg->sum(mesh->Nlocal, o_tmp1, platform->comm.mpiComm);
    linAlg->axpby(mesh->Nlocal, -prhs, o_tmp2, 1.0, o_div);
    o_tmp1.free();
    o_tmp2.free();

    const auto *coeff = nrs->coeffBDF;
    dfloat Saqpq = 0.0;
    for (int i = 0; i < nrs->nBDF; ++i) {
      Saqpq += coeff[i] * nrs->p0th[i];
    }

    const auto g0 = nrs->g0;
    const auto dt = nrs->dt[0];

    const auto pcoef = (g0 - dt * prhs);
    const auto p0thn = Saqpq / pcoef;

    nrs->p0th[2] = nrs->p0th[1];
    nrs->p0th[1] = nrs->p0th[0];
    nrs->p0th[0] = p0thn;

    nrs->dp0thdt = prhs * p0thn;

    surfaceFlops += surfaceFluxFlops + p0thHelperFlops;
  }

  qThermal = 0;

  double flops = surfaceFlops + flopsGrad + flopsQTL;
  platform->flopCounter->add("lowMach::qThermalRealGasSingleComponent", flops);
}

void lowMach::dpdt(occa::memory &o_FU)
{
  nekrsCheck(!setupCalled || !buildKernelCalled,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s\n",
             "called prior to tavg::setup()!");

  auto nrs = _nrs;
  auto mesh = nrs->mesh;

  if (nrs->cds->cvodeSolve[0]) {
    return; // contribution is not applied here
  }

  if (!qThermal) {
    platform->linAlg->add(mesh->Nlocal, nrs->dp0thdt * alpha0, o_FU);
  }
}

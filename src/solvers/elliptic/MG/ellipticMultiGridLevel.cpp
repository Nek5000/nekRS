/*

   The MIT License (MIT)

   Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

 */
#include <type_traits>
#include "elliptic.h"
#include "ellipticPrecon.h"
#include "ellipticMultiGrid.h"
#include "linAlg.hpp"
#include <iostream>
void pMGLevel::Ax(occa::memory o_x, occa::memory o_Ax)
{
  ellipticOperator(elliptic,o_x,o_Ax, pfloatString);
}

void pMGLevel::residual(occa::memory o_rhs, occa::memory o_x, occa::memory o_res)
{
  ellipticOperator(elliptic, o_x, o_res, pfloatString);
  platform->linAlg->paxpbyMany(Nrows, elliptic->Nfields, elliptic->fieldOffset, 1.0, o_rhs, -1.0, o_res);
}

void pMGLevel::coarsen(occa::memory o_x, occa::memory o_Rx)
{
  double flopCounter = 0.0;
  if (options.compareArgs("DISCRETIZATION", "CONTINUOUS")) {
    platform->linAlg->paxmy(mesh->Nelements * NpF, 1.0, o_invDegreeFine, o_x);
    flopCounter += static_cast<double>(mesh->Nelements) * NpF;
  }

  const auto NqC = elliptic->mesh->Nq;
  const auto NqF = std::cbrt(NpF);

  precon_t *precon = elliptic->precon;

  precon->coarsenKernel(mesh->Nelements, o_R, o_x, o_Rx);
  const auto workPerElem = 2 * (NqF * NqF * NqF * NqC + NqF * NqF * NqC * NqC + NqF * NqC * NqC * NqC);
  flopCounter += static_cast<double>(mesh->Nelements) * workPerElem;

  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    oogs::startFinish(o_Rx, elliptic->Nfields, elliptic->fieldOffset, ogsPfloat, ogsAdd, elliptic->oogs);
    ellipticApplyMask(elliptic, o_Rx, pfloatString); // apply mask again because coarsenKernel do not preserve it
  }

  const double factor = std::is_same<pfloat, float>::value ? 0.5 : 1.0;
  platform->flopCounter->add("pMGLevel::coarsen, N=" + std::to_string(mesh->N), factor * flopCounter);

}

void pMGLevel::prolongate(occa::memory o_x, occa::memory o_Px)
{
  precon_t *precon = elliptic->precon;

  precon->prolongateKernel(mesh->Nelements, o_R, o_x, o_Px);
  const auto NqC = elliptic->mesh->Nq;
  const auto NqF = std::cbrt(NpF);
  double flopCounter = 2 * (NqF * NqF * NqF * NqC + NqF * NqF * NqC * NqC + NqF * NqC * NqC * NqC);
  flopCounter += NqF * NqF * NqF;
  flopCounter *= static_cast<double>(mesh->Nelements);

  const double factor = std::is_same<pfloat, float>::value ? 0.5 : 1.0;
  platform->flopCounter->add("pMGLevel::prolongate, N=" + std::to_string(mesh->N), factor * flopCounter);
}

void pMGLevel::smooth(occa::memory o_rhs, occa::memory o_x, bool x_is_zero)
{
  platform->timer.tic(elliptic->name + " preconditioner smoother N=" + std::to_string(mesh->N), 1);

  if(!x_is_zero && smootherType == SmootherType::ASM) return;
  if(!x_is_zero && smootherType == SmootherType::RAS) return;

  if (smootherType == SmootherType::CHEBYSHEV)
    this->smoothChebyshev(o_rhs, o_x, x_is_zero);
  else if (smootherType == SmootherType::OPT_FOURTH_CHEBYSHEV || smootherType == SmootherType::FOURTH_CHEBYSHEV)
    this->smoothFourthKindChebyshev(o_rhs, o_x, x_is_zero);
  else if (smootherType == SmootherType::ASM)
    this->smoothSchwarz(o_rhs, o_x, x_is_zero);
  else if (smootherType == SmootherType::RAS)
    this->smoothSchwarz(o_rhs, o_x, x_is_zero);
  else if (smootherType == SmootherType::JACOBI)
    this->smoothJacobi(o_rhs, o_x, x_is_zero);

  platform->timer.toc(elliptic->name + " preconditioner smoother N=" + std::to_string(mesh->N));
}

void pMGLevel::smoother(occa::memory o_x, occa::memory o_Sx, bool x_is_zero)
{
  if (chebySmootherType == ChebyshevSmootherType::JACOBI){
    this->smootherJacobi(o_x, o_Sx);
  } else {
    this->smoothSchwarz(o_x, o_Sx, true);
  }
}

void pMGLevel::smoothJacobi (occa::memory &o_r, occa::memory &o_x, bool xIsZero)
{
  occa::memory o_res = o_smootherResidual;
  occa::memory o_Ad  = o_smootherResidual2;
  occa::memory o_d   = o_smootherUpdate;

  const pfloat one = 1.0;
  const pfloat mone = -1.0;
  const pfloat zero = 0.0;

  double flopCount = 0.0;

  if(xIsZero) { //skip the Ax if x is zero
    //res = Sr
    platform->linAlg->paxmyz(Nrows,one,o_invDiagA,o_r,o_x);
    flopCount += Nrows;
  } else {
    //res = S(r-Ax)
    this->Ax(o_x,o_res);
    platform->linAlg->paxpby(Nrows, one, o_r, mone, o_res);
    platform->linAlg->paxmyz(Nrows, one, o_invDiagA, o_res, o_d);
    platform->linAlg->paxpby(Nrows, one, o_d, one, o_x);
    // two saxpy's + collocation
    flopCount += 7 * Nrows;
  }
  auto mesh = elliptic->mesh;
  const double factor = std::is_same<pfloat, float>::value ? 0.5 : 1.0;
  platform->flopCounter->add("pMGLevel::smoothJacobi, N=" + std::to_string(mesh->N), factor * flopCount);
}

void pMGLevel::smoothChebyshev (occa::memory &o_r, occa::memory &o_x, bool xIsZero)
{
  const auto ChebyshevDegree = xIsZero ? DownLegChebyshevDegree : UpLegChebyshevDegree;

  // p_0(0) = I -> no-op smoothing
  if (ChebyshevDegree == 0)
    return;

  const pfloat theta = 0.5 * (lambda1 + lambda0);
  const pfloat delta = 0.5 * (lambda1 - lambda0);
  const pfloat invTheta = 1.0 / theta;
  const pfloat sigma = theta / delta;
  pfloat rho_n = 1. / sigma;
  pfloat rho_np1;

  pfloat one = 1., mone = -1., zero = 0.0;

  occa::memory o_res = o_smootherResidual;
  occa::memory o_Ad  = o_smootherResidual2;
  occa::memory o_d   = o_smootherUpdate;

  double flopCount = 0.0;

  if (xIsZero) {
    platform->linAlg->pfill(Nrows, zero, o_x);
  }

  // res = S(r-Ax)
  if (!xIsZero) {
    this->Ax(o_x,o_res);
    platform->linAlg->paxpby(Nrows, one, o_r, mone, o_res);
    flopCount += 2 * Nrows;
  } else {
    o_res.copyFrom(o_r, Nrows*sizeof(pfloat));
  }
  this->smoother(o_res, o_res, xIsZero);

  // d = invTheta*res
  platform->linAlg->paxpby(Nrows, invTheta, o_res, zero, o_d);
  flopCount += Nrows;

  for (int k = 1; k < ChebyshevDegree; k++) {

    // SAd_k
    this->Ax(o_d,o_Ad);
    this->smoother(o_Ad, o_Ad, xIsZero);

    // x_k+1 = x_k + d_k
    // r_k+1 = r_k - SAd_k
    // d_k+1 = (rho_k+1*rho_k)*d_k  + (2*rho_k+1/delta)*r_k+1

    const pfloat rhoSave = rho_n;
    rho_n = 1.0 / (2.0 * sigma - rho_n);

    const pfloat rCoeff = 2.0 * rho_n / delta;
    const pfloat dCoeff = rho_n * rhoSave;

    elliptic->updateChebyshevKernel(Nrows, dCoeff, rCoeff, o_Ad, o_d, o_res, o_x);

    flopCount += 5 * Nrows;
  }
  //x_k+1 = x_k + d_k
  platform->linAlg->paxpby(Nrows, one, o_d, one, o_x);
  flopCount += Nrows;
  ellipticApplyMask(elliptic, o_x, pfloatString);
  const double factor = std::is_same<pfloat, float>::value ? 0.5 : 1.0;
  platform->flopCounter->add("pMGLevel::smoothChebyshev, N=" + std::to_string(mesh->N), factor * flopCount);
}

void pMGLevel::smoothFourthKindChebyshev (occa::memory &o_r, occa::memory &o_x, bool xIsZero)
{
  const auto ChebyshevDegree = xIsZero ? DownLegChebyshevDegree : UpLegChebyshevDegree;
  auto &betas = xIsZero ? DownLegBetas : UpLegBetas;

  // p_0(0) = I -> no-op smoothing
  if (ChebyshevDegree == 0)
    return;

  pfloat one = 1., mone = -1., zero = 0.0;

  occa::memory o_res = o_smootherResidual;
  occa::memory o_Ad = o_smootherResidual2;
  occa::memory o_d = o_smootherUpdate;

  const auto rho = this->lambda1;

  double flopCount = 0.0;

  // r = b - Ax
  if (xIsZero) {
    platform->linAlg->pfill(Nrows, zero, o_x);
    o_res.copyFrom(o_r, Nrows * sizeof(pfloat));
  }
  else {
    this->Ax(o_x,o_res);
    platform->linAlg->paxpby(Nrows, one, o_r, mone, o_res);
    flopCount += Nrows;
  }

  // d = \dfrac{4}{3} \dfrac{1}{\rho(SA)} Sr
  this->smoother(o_res, o_Ad, xIsZero);
  const pfloat coeff = 4.0 / (3.0 * rho);
  platform->linAlg->paxpby(Nrows, coeff, o_Ad, zero, o_d);

  for (int k = 1; k < ChebyshevDegree; k++) {

    // Ad_k
    this->Ax(o_d, o_Ad);

    // x_k+1 = x_k + \beta_k d_k
    // r_k+1 = r_k - Ad_k
    elliptic->updateFourthKindChebyshevKernel(Nrows, betas[k - 1], o_Ad, o_d, o_res, o_x);

    this->smoother(o_res, o_Ad, xIsZero);

    // d_k+1 = \dfrac{2k-1}{2k+3} d_k + \dfrac{8k+4}{2k+3} \dfrac{1}{\rho(SA)} S r_k+1
    const pfloat dCoeff = (2.0 * k - 1.0) / (2.0 * k + 3.0);
    const pfloat rCoeff = (8.0 * k + 4.0) / ((2.0 * k + 3.0) * rho);
    platform->linAlg->paxpby(Nrows, rCoeff, o_Ad, dCoeff, o_d);
  }

  //x_k+1 = x_k + \beta_k d_k
  platform->linAlg->paxpby(Nrows, betas.back(), o_d, one, o_x);
  flopCount += 2 * Nrows;
  ellipticApplyMask(elliptic, o_x, pfloatString);
  const double factor = std::is_same<pfloat, float>::value ? 0.5 : 1.0;
  platform->flopCounter->add("pMGLevel::smoothOptChebyshev, N=" + std::to_string(mesh->N), factor * flopCount);
}

void pMGLevel::smootherJacobi(occa::memory &o_r, occa::memory &o_Sr)
{
  platform->linAlg->paxmyz(Nrows, 1.0f, o_invDiagA, o_r, o_Sr);
  const double factor = std::is_same<pfloat, float>::value ? 0.5 : 1.0;
  platform->flopCounter->add("pMGLevel::smootherJacobi, N=" + std::to_string(mesh->N), factor * Nrows);
}

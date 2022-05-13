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
#include "linAlg.hpp"
#include <iostream>
void MGLevel::Ax(occa::memory o_x, occa::memory o_Ax)
{
  ellipticOperator(elliptic,o_x,o_Ax, pfloatString);
}

void MGLevel::residual(occa::memory o_rhs, occa::memory o_x, occa::memory o_res)
{
  if(stype != SmootherType::SCHWARZ) {
    ellipticOperator(elliptic,o_x,o_res, dfloatString);
    // subtract r = b - A*x
    platform->linAlg->axpbyMany(Nrows, elliptic->Nfields, elliptic->Ntotal, 1.0, o_rhs, -1.0, o_res);
  } else {
    o_res.copyFrom(o_rhs, Nrows*sizeof(dfloat));
  }
}

void MGLevel::coarsen(occa::memory o_x, occa::memory o_Rx)
{
  double flopCounter = 0.0;
  if (options.compareArgs("DISCRETIZATION", "CONTINUOUS")) {
    platform->linAlg->axmy(mesh->Nelements * NpF, 1.0, o_invDegree, o_x);
    flopCounter += static_cast<double>(mesh->Nelements) * NpF;
  }

  const auto NqC = elliptic->mesh->Nq;
  const auto NqF = std::cbrt(NpF);

  elliptic->precon->coarsenKernel(mesh->Nelements, o_R, o_x, o_Rx);
  const auto workPerElem = 2 * (NqF * NqF * NqF * NqC + NqF * NqF * NqC * NqC + NqF * NqC * NqC * NqC);
  flopCounter += static_cast<double>(mesh->Nelements) * workPerElem;

  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    oogs::startFinish(o_Rx, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, ogsAdd, elliptic->oogs);
    // ellipticApplyMask(elliptic, o_Rx, dfloatString);
  }

  platform->flopCounter->add("MGLevel::coarsen, N=" + std::to_string(mesh->N), flopCounter);
}

void MGLevel::prolongate(occa::memory o_x, occa::memory o_Px)
{
  elliptic->precon->prolongateKernel(mesh->Nelements, o_R, o_x, o_Px);
  const auto NqC = elliptic->mesh->Nq;
  const auto NqF = std::cbrt(NpF);
  double flopCounter = 2 * (NqF * NqF * NqF * NqC + NqF * NqF * NqC * NqC + NqF * NqC * NqC * NqC);
  flopCounter += NqF * NqF * NqF;
  flopCounter *= static_cast<double>(mesh->Nelements);
  platform->flopCounter->add("MGLevel::prolongate, N=" + std::to_string(mesh->N), flopCounter);
}

void MGLevel::smooth(occa::memory o_rhs, occa::memory o_x, bool x_is_zero)
{
  platform->timer.tic(elliptic->name + " preconditioner smoother", 1);
  if(!x_is_zero && stype == SmootherType::SCHWARZ) return;
  if(!strstr(pfloatString,dfloatString)) {
    elliptic->fusedCopyDfloatToPfloatKernel(Nrows, o_x, o_rhs, o_xPfloat, o_rhsPfloat);
    if (stype == SmootherType::CHEBYSHEV)
      this->smoothChebyshev(o_rhsPfloat, o_xPfloat, x_is_zero);
    else if (stype == SmootherType::SCHWARZ)
      this->smoothSchwarz(o_rhsPfloat, o_xPfloat, x_is_zero);
    else if (stype == SmootherType::JACOBI)
      this->smoothJacobi(o_rhsPfloat, o_xPfloat, x_is_zero);
      
    elliptic->copyPfloatToDPfloatKernel(Nrows, o_xPfloat, o_x);
  } else {
    if (stype == SmootherType::CHEBYSHEV)
      this->smoothChebyshev(o_rhs, o_x, x_is_zero);
    else if (stype == SmootherType::SCHWARZ)
      this->smoothSchwarz(o_rhs, o_x, x_is_zero);
    else if (stype == SmootherType::JACOBI)
      this->smoothJacobi(o_rhs, o_x, x_is_zero);
  }
  platform->timer.toc(elliptic->name + " preconditioner smoother");
}

void MGLevel::smoother(occa::memory o_x, occa::memory o_Sx, bool x_is_zero)
{
  // x_is_zero = true <-> downward leg
  if(x_is_zero) {
    if (smtypeDown == SecondarySmootherType::JACOBI)
      this->smootherJacobi(o_x, o_Sx);
    else
      this->smoothSchwarz(o_x, o_Sx, true); // no-op if false
  } else {
    if (smtypeUp == SecondarySmootherType::JACOBI)
      this->smootherJacobi(o_x, o_Sx);
    else
      this->smoothSchwarz(o_x, o_Sx, true); // no-op if false
  }
}

void MGLevel::smoothJacobi (occa::memory &o_r, occa::memory &o_x, bool xIsZero)
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
    elliptic->dotMultiplyPfloatKernel(Nrows,o_invDiagA,o_r,o_x);
    flopCount += Nrows;
  } else {
    //res = S(r-Ax)
    this->Ax(o_x,o_res);
    elliptic->scaledAddPfloatKernel(Nrows, one, o_r, mone, o_res);
    elliptic->dotMultiplyPfloatKernel(Nrows, o_invDiagA, o_res, o_d);
    elliptic->scaledAddPfloatKernel(Nrows, one, o_d, one, o_x);
    // two saxpy's + collocation
    flopCount += 7 * Nrows;
  }
  auto mesh = elliptic->mesh;
  const double factor = std::is_same<pfloat, float>::value ? 0.5 : 1.0;
  platform->flopCounter->add("MGLevel::smoothJacobi, N=" + std::to_string(mesh->N), factor * flopCount);
}
void MGLevel::smoothChebyshevOneIteration (occa::memory &o_r, occa::memory &o_x, bool xIsZero)
{
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

  if(xIsZero) { //skip the Ax if x is zero
    //res = Sr
    this->smoother(o_r, o_res, xIsZero);
    elliptic->updateSmoothedSolutionVecKernel(Nrows, invTheta, o_res, one, o_d, zero, o_x);
    flopCount += 4 * Nrows;
  } else {
    //res = S(r-Ax)
    this->Ax(o_x,o_res);
    elliptic->scaledAddPfloatKernel(Nrows, one, o_r, mone, o_res);
    this->smoother(o_res, o_res, xIsZero);
    elliptic->updateSmoothedSolutionVecKernel(Nrows, invTheta, o_res, one, o_d, one, o_x);

    flopCount += 7 * Nrows;
  }

  //r_k+1 = r_k - SAd_k
  this->Ax(o_d,o_Ad);
  this->smoother(o_Ad, o_Ad, xIsZero);
  rho_np1 = 1.0 / (2. * sigma - rho_n);
  pfloat rhoDivDelta = 2.0 * rho_np1 / delta;
  elliptic->updateChebyshevSolutionVecKernel(Nrows, rhoDivDelta, rho_np1, rho_n, o_Ad, o_res, o_d, o_x);
  ellipticApplyMask(elliptic, o_x, pfloatString);

  flopCount += 6 * Nrows;

  auto mesh = elliptic->mesh;
  const double factor = std::is_same<pfloat, float>::value ? 0.5 : 1.0;
  platform->flopCounter->add("MGLevel::smoothChebyshevOneIteration, N=" + std::to_string(mesh->N),
                             factor * flopCount);
}
void MGLevel::smoothChebyshevTwoIteration (occa::memory &o_r, occa::memory &o_x, bool xIsZero)
{
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

  if(xIsZero) { //skip the Ax if x is zero
    //res = Sr
    this->smoother(o_r, o_res, xIsZero);

    elliptic->updateSmoothedSolutionVecKernel(Nrows, invTheta, o_res, one, o_d, zero, o_x);
    flopCount += 4 * Nrows;
  } else {
    //res = S(r-Ax)
    this->Ax(o_x,o_res);
    elliptic->scaledAddPfloatKernel(Nrows, one, o_r, mone, o_res);
    this->smoother(o_res, o_res, xIsZero);

    elliptic->updateSmoothedSolutionVecKernel(Nrows, invTheta, o_res, one, o_d, one, o_x);
    flopCount += 7 * Nrows;
  }


  //r_k+1 = r_k - SAd_k
  this->Ax(o_d,o_Ad);
  this->smoother(o_Ad, o_Ad, xIsZero);
  rho_np1 = 1.0 / (2. * sigma - rho_n);
  pfloat rhoDivDelta = 2.0 * rho_np1 / delta;

  elliptic->updateIntermediateSolutionVecKernel(Nrows, rhoDivDelta, rho_n, rho_np1, o_Ad, o_res, o_d, o_x);
  flopCount += 6 * Nrows;

  rho_n = rho_np1;
  //r_k+1 = r_k - SAd_k
  this->Ax(o_d,o_Ad);
  this->smoother(o_Ad, o_Ad, xIsZero);
  rho_np1 = 1.0 / (2. * sigma - rho_n);
  rhoDivDelta = 2.0 * rho_np1 / delta;

  elliptic->updateIntermediateSolutionVecKernel(Nrows, rhoDivDelta, rho_n, rho_np1, o_Ad, o_res, o_d, o_x);
  flopCount += 6 * Nrows;

  ellipticApplyMask(elliptic, o_x, pfloatString);
  const double factor = std::is_same<pfloat, float>::value ? 0.5 : 1.0;
  platform->flopCounter->add("MGLevel::smoothChebyshevTwoIteration, N=" + std::to_string(mesh->N),
                             factor * flopCount);
}

void MGLevel::smoothChebyshev (occa::memory &o_r, occa::memory &o_x, bool xIsZero)
{
  if(ChebyshevIterations == 1) {
    smoothChebyshevOneIteration(o_r,o_x,xIsZero);
    return;
  } else if (ChebyshevIterations == 2) {
    smoothChebyshevTwoIteration(o_r,o_x,xIsZero);
    return;
  }
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

  if(xIsZero) { //skip the Ax if x is zero
    //res = Sr
    this->smoother(o_r, o_res, xIsZero);

    //d = invTheta*res
    elliptic->scaledAddPfloatKernel(Nrows, invTheta, o_res, zero, o_d);
    flopCount += Nrows;
  } else {
    //res = S(r-Ax)
    this->Ax(o_x,o_res);
    elliptic->scaledAddPfloatKernel(Nrows, one, o_r, mone, o_res);
    this->smoother(o_res, o_res, xIsZero);
    flopCount += 2 * Nrows;

    //d = invTheta*res
    elliptic->scaledAddPfloatKernel(Nrows, invTheta, o_res, zero, o_d);
    flopCount += Nrows;
  }

  for (int k = 0; k < ChebyshevIterations; k++) {
    //x_k+1 = x_k + d_k
    if (xIsZero && (k == 0)) {
      elliptic->scaledAddPfloatKernel(Nrows, one, o_d, zero, o_x);
    }
    else {
      elliptic->scaledAddPfloatKernel(Nrows, one, o_d, one, o_x);
      flopCount += 1 * Nrows;
    }

    //r_k+1 = r_k - SAd_k
    this->Ax(o_d,o_Ad);
    this->smoother(o_Ad, o_Ad, xIsZero);
    elliptic->scaledAddPfloatKernel(Nrows, mone, o_Ad, one, o_res);
    flopCount += Nrows;

    rho_np1 = 1.0 / (2. * sigma - rho_n);
    pfloat rhoDivDelta = 2.0 * rho_np1 / delta;

    //d_k+1 = rho_k+1*rho_k*d_k  + 2*rho_k+1*r_k+1/delta
    elliptic->scaledAddPfloatKernel(Nrows, rhoDivDelta, o_res, rho_np1 * rho_n, o_d);
    flopCount += 4 * Nrows;

    rho_n = rho_np1;
  }
  //x_k+1 = x_k + d_k
  elliptic->scaledAddPfloatKernel(Nrows, one, o_d, one, o_x);
  flopCount += Nrows;
  ellipticApplyMask(elliptic, o_x, pfloatString);
  const double factor = std::is_same<pfloat, float>::value ? 0.5 : 1.0;
  platform->flopCounter->add("MGLevel::smoothChebyshev, N=" + std::to_string(mesh->N), factor * flopCount);
}

void MGLevel::smootherJacobi(occa::memory &o_r, occa::memory &o_Sr)
{
  elliptic->dotMultiplyPfloatKernel(Nrows, o_invDiagA, o_r, o_Sr);
  const double factor = std::is_same<pfloat, float>::value ? 0.5 : 1.0;
  platform->flopCounter->add("MGLevel::smootherJacobi, N=" + std::to_string(mesh->N), factor * Nrows);
}

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
    const dlong Nlocal = mesh->Np * mesh->Nelements;
    platform->linAlg->axpbyMany(
      Nlocal,
      elliptic->Nfields,
      elliptic->Ntotal,
      1.0,
      o_rhs,
      -1.0,
      o_res
    );
  } else {
    o_res.copyFrom(o_rhs, Nrows*sizeof(dfloat));
  }
}

void MGLevel::coarsen(occa::memory o_x, occa::memory o_Rx)
{
  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) 
    platform->linAlg->axmy(mesh->Nelements * NpF, 1.0, o_invDegree, o_x);

  elliptic->precon->coarsenKernel(mesh->Nelements, o_R, o_x, o_Rx);

  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    oogs::startFinish(o_Rx, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, ogsAdd, elliptic->oogs);
    //if (elliptic->Nmasked) mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, o_Rx);
  }
}

void MGLevel::prolongate(occa::memory o_x, occa::memory o_Px)
{
  elliptic->precon->prolongateKernel(mesh->Nelements, o_R, o_x, o_Px);
}

void MGLevel::smooth(occa::memory o_rhs, occa::memory o_x, bool x_is_zero)
{
  platform->timer.tic(elliptic->name + " preconditioner smoother", 1);
  if(!x_is_zero && stype == SmootherType::SCHWARZ) return;
  if(!strstr(pfloatString,dfloatString)) {
    elliptic->fusedCopyDfloatToPfloatKernel(Nrows, o_x, o_rhs, o_xPfloat, o_rhsPfloat);
    if (stype == SmootherType::CHEBYSHEV)
      this->smoothChebyshev(o_rhsPfloat, o_xPfloat, x_is_zero);
    else
      this->smoothSchwarz(o_rhsPfloat, o_xPfloat, x_is_zero);
    elliptic->copyPfloatToDPfloatKernel(Nrows, o_xPfloat, o_x);
  } else {
    if (stype == SmootherType::CHEBYSHEV)
      this->smoothChebyshev(o_rhs, o_x, x_is_zero);
    else
      this->smoothSchwarz(o_rhs, o_x, x_is_zero);
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

  if(xIsZero) { //skip the Ax if x is zero
    //res = Sr
    this->smoother(o_r, o_res, xIsZero);
    elliptic->updateSmoothedSolutionVecKernel(Nrows, invTheta, o_res, one, o_d, zero, o_x);
  } else {
    //res = S(r-Ax)
    this->Ax(o_x,o_res);
    elliptic->scaledAddPfloatKernel(Nrows, one, o_r, mone, o_res);
    this->smoother(o_res, o_res, xIsZero);
    elliptic->updateSmoothedSolutionVecKernel(Nrows, invTheta, o_res, one, o_d, one, o_x);
  }

  //r_k+1 = r_k - SAd_k
  this->Ax(o_d,o_Ad);
  this->smoother(o_Ad, o_Ad, xIsZero);
  rho_np1 = 1.0 / (2. * sigma - rho_n);
  pfloat rhoDivDelta = 2.0 * rho_np1 / delta;
  elliptic->updateChebyshevSolutionVecKernel(Nrows, rhoDivDelta, rho_np1, rho_n, o_Ad, o_res, o_d, o_x);
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

  if(xIsZero) { //skip the Ax if x is zero
    //res = Sr
    this->smoother(o_r, o_res, xIsZero);

    elliptic->updateSmoothedSolutionVecKernel(Nrows, invTheta, o_res, one, o_d, zero, o_x);
  } else {
    //res = S(r-Ax)
    this->Ax(o_x,o_res);
    elliptic->scaledAddPfloatKernel(Nrows, one, o_r, mone, o_res);
    this->smoother(o_res, o_res, xIsZero);

    elliptic->updateSmoothedSolutionVecKernel(Nrows, invTheta, o_res, one, o_d, one, o_x);
  }


  //r_k+1 = r_k - SAd_k
  this->Ax(o_d,o_Ad);
  this->smoother(o_Ad, o_Ad, xIsZero);
  rho_np1 = 1.0 / (2. * sigma - rho_n);
  pfloat rhoDivDelta = 2.0 * rho_np1 / delta;

  elliptic->updateIntermediateSolutionVecKernel(Nrows, rhoDivDelta, rho_n, rho_np1, o_Ad, o_res, o_d, o_x);

  rho_n = rho_np1;
  //r_k+1 = r_k - SAd_k
  this->Ax(o_d,o_Ad);
  this->smoother(o_Ad, o_Ad, xIsZero);
  rho_np1 = 1.0 / (2. * sigma - rho_n);
  rhoDivDelta = 2.0 * rho_np1 / delta;

  elliptic->updateIntermediateSolutionVecKernel(Nrows, rhoDivDelta, rho_n, rho_np1, o_Ad, o_res, o_d, o_x);

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

  if(xIsZero) { //skip the Ax if x is zero
    //res = Sr
    this->smoother(o_r, o_res, xIsZero);

    //d = invTheta*res
    elliptic->scaledAddPfloatKernel(Nrows, invTheta, o_res, zero, o_d);
  } else {
    //res = S(r-Ax)
    this->Ax(o_x,o_res);
    elliptic->scaledAddPfloatKernel(Nrows, one, o_r, mone, o_res);
    this->smoother(o_res, o_res, xIsZero);

    //d = invTheta*res
    elliptic->scaledAddPfloatKernel(Nrows, invTheta, o_res, zero, o_d);
  }

  for (int k = 0; k < ChebyshevIterations; k++) {
    //x_k+1 = x_k + d_k
    if (xIsZero && (k == 0))
      elliptic->scaledAddPfloatKernel(Nrows, one, o_d, zero, o_x);
    else
      elliptic->scaledAddPfloatKernel(Nrows, one, o_d, one, o_x);

    //r_k+1 = r_k - SAd_k
    this->Ax(o_d,o_Ad);
    this->smoother(o_Ad, o_Ad, xIsZero);
    elliptic->scaledAddPfloatKernel(Nrows, mone, o_Ad, one, o_res);

    rho_np1 = 1.0 / (2. * sigma - rho_n);
    pfloat rhoDivDelta = 2.0 * rho_np1 / delta;

    //d_k+1 = rho_k+1*rho_k*d_k  + 2*rho_k+1*r_k+1/delta
    elliptic->scaledAddPfloatKernel(Nrows, rhoDivDelta, o_res, rho_np1 * rho_n, o_d);

    rho_n = rho_np1;
  }
  //x_k+1 = x_k + d_k
  elliptic->scaledAddPfloatKernel(Nrows, one, o_d, one, o_x);
}

void MGLevel::smootherJacobi(occa::memory &o_r, occa::memory &o_Sr)
{
  elliptic->dotMultiplyPfloatKernel(mesh->Np * mesh->Nelements,o_invDiagA,o_r,o_Sr);
}

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
#include <iostream>
void MGLevel::Ax(occa::memory o_x, occa::memory o_Ax)
{
  ellipticOperator(elliptic,o_x,o_Ax, pfloatString);
}

void MGLevel::residual(occa::memory o_rhs, occa::memory o_x, occa::memory o_res)
{
  if(stype != SCHWARZ) {
    ellipticOperator(elliptic,o_x,o_res, dfloatString);
    // subtract r = b - A*x
    ellipticScaledAdd(elliptic, 1.f, o_rhs, -1.f, o_res);
  } else {
    //o_res.copyFrom(o_rhs, Nrows*sizeof(dfloat));
    dfloat zero = 0.0;
    dfloat one = 1.0;
    elliptic->scaledAddKernel(Nrows, one, o_rhs, zero, o_res);
  }
}

void MGLevel::coarsen(occa::memory o_x, occa::memory o_Rx)
{
  if (options.compareArgs("DISCRETIZATION","CONTINUOUS"))
    elliptic->dotMultiplyKernel(mesh->Nelements * NpF, o_invDegree, o_x, o_x);

  elliptic->precon->coarsenKernel(mesh->Nelements, o_R, o_x, o_Rx);

  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    oogs::startFinish(o_Rx, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, ogsAdd, elliptic->oogs);
    if (elliptic->Nmasked) mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, o_Rx);
  }
}

void MGLevel::prolongate(occa::memory o_x, occa::memory o_Px)
{
  elliptic->precon->prolongateKernel(mesh->Nelements, o_R, o_x, o_Px);
}

void MGLevel::smooth(occa::memory o_rhs, occa::memory o_x, bool x_is_zero)
{
  if(!strstr(pfloatString,dfloatString)){
    elliptic->copyDfloatToPfloatKernel(Nrows, o_xPfloat, o_x);
    elliptic->copyDfloatToPfloatKernel(Nrows, o_rhsPfloat, o_rhs);
    if (stype == RICHARDSON)
      this->smoothRichardson(o_rhsPfloat, o_xPfloat, x_is_zero);
    else if (stype == CHEBYSHEV)
      this->smoothChebyshev(o_rhsPfloat, o_xPfloat, x_is_zero);
    else if (stype == SCHWARZ)
      this->smoothSchwarz(o_rhsPfloat, o_xPfloat, x_is_zero);
    elliptic->copyPfloatToDPfloatKernel(Nrows, o_xPfloat, o_x);
    elliptic->copyPfloatToDPfloatKernel(Nrows, o_rhsPfloat, o_rhs);
  } else {
    if (stype == RICHARDSON)
      this->smoothRichardson(o_rhs, o_x, x_is_zero);
    else if (stype == CHEBYSHEV)
      this->smoothChebyshev(o_rhs, o_x, x_is_zero);
    else if (stype == SCHWARZ)
      this->smoothSchwarz(o_rhs, o_x, x_is_zero);
  }
}

void MGLevel::smoother(occa::memory o_x, occa::memory o_Sx, bool x_is_zero)
{
  // x_is_zero = true <-> downward leg
  if(x_is_zero) {
    if (smtypeDown == JACOBI)
      this->smootherJacobi(o_x, o_Sx);
    else if (smtypeDown == SCHWARZ_SMOOTH)
      //this->smootherSchwarz(o_x, o_Sx);
      this->smoothSchwarz(o_x, o_Sx, true); // no-op if false
  } else {
    if (smtypeUp == JACOBI)
      this->smootherJacobi(o_x, o_Sx);
    else if (smtypeUp == SCHWARZ_SMOOTH)
      //this->smootherSchwarz(o_x, o_Sx);
      this->smoothSchwarz(o_x, o_Sx, true); // no-op if false
  }
}

void MGLevel::smoothRichardson(occa::memory &o_r, occa::memory &o_x, bool xIsZero)
{
  occa::memory o_res = o_smootherResidual;

  if (xIsZero) {
    this->smoother(o_r, o_x, xIsZero);
    return;
  }

  pfloat one = 1.;
  pfloat mone = -1.;

  //res = r-Ax
  this->Ax(o_x,o_res);
  elliptic->scaledAddPfloatKernel(Nrows, one, o_r, mone, o_res);

  //smooth the fine problem x = x + S(r-Ax)
  this->smoother(o_res, o_res, xIsZero);
  elliptic->scaledAddPfloatKernel(Nrows, one, o_res, one, o_x);
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
void MGLevel::smoothChebyshev (occa::memory &o_r, occa::memory &o_x, bool xIsZero)
{
  if(ChebyshevIterations == 1){
    smoothChebyshevOneIteration(o_r,o_x,xIsZero);
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

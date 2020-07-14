/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#include "parAlmond.hpp"

namespace parAlmond {

void agmgLevel::smoothJacobi(dfloat *r, dfloat *x,
                             const bool x_is_zero) {

  // x = x + inv(D)*(b-A*x)
  if(x_is_zero){
    vectorDotStar(Nrows,1.0,A->diagInv,r,0.0,x);
    return;
  }

  dfloat *res = (dfloat *) scratch;

  A->SpMV(-1.0, x, 1.0, r, res);
  vectorDotStar(Nrows, 1.0, A->diagInv, res, 1.0, x);
}


void agmgLevel::smoothDampedJacobi(dfloat *r, dfloat *x,
                                   const bool x_is_zero) {

  // x = x + alpha*inv(D)*(b-A*x)
  if(x_is_zero){
    vectorDotStar(Nrows,lambda,A->diagInv,r,0.0,x);
    return;
  }

  dfloat *res = (dfloat *) scratch;

  A->SpMV(-1.0, x, 1.0, r, res);
  vectorDotStar(Nrows, lambda, A->diagInv, res, 1.0, x);
}

void agmgLevel::smoothChebyshev(dfloat *r, dfloat *x,
                                const bool x_is_zero) {

  const dfloat theta = 0.5*(lambda1+lambda0);
  const dfloat delta = 0.5*(lambda1-lambda0);
  const dfloat invTheta = 1.0/theta;
  const dfloat sigma = theta/delta;
  dfloat rho_n = 1./sigma;
  dfloat rho_np1;

  dfloat *res = ((dfloat*) scratch) + 0*Ncols;
  dfloat *Ad  = ((dfloat*) scratch) + 1*Ncols;
  dfloat *d   = ((dfloat*) scratch) + 2*Ncols;

  if(x_is_zero){ //skip the Ax if x is zero
    //res = D^{-1}r
    vectorDotStar(Nrows, 1.0, A->diagInv, r, 0.0, res);
    vectorSet(Nrows, 0.0, x);
    //d = invTheta*res
    vectorAdd(Nrows, invTheta, res, 0.0, d);
  } else {
    //res = D^{-1}(r-Ax)
    A->SpMV(-1.0, x, 1.0, r, res);
    vectorDotStar(Nrows, A->diagInv, res);

    //d = invTheta*res
    vectorAdd(Nrows, invTheta, res, 0.0, d);
  }

  for (int k=0;k<ChebyshevIterations;k++) {
    //x_k+1 = x_k + d_k
    vectorAdd(Nrows, 1.0, d, 1.0, x);

    //r_k+1 = r_k - D^{-1}Ad_k
    A->SpMV(1.0, d, 0.0, Ad);
    vectorDotStar(Nrows, -1.0, A->diagInv, Ad, 1.0, res);

    rho_np1 = 1.0/(2.*sigma-rho_n);

    //d_k+1 = rho_k+1*rho_k*d_k  + 2*rho_k+1*r_k+1/delta
    vectorAdd(Nrows, 2.0*rho_np1/delta, res, rho_np1*rho_n, d);
    rho_n = rho_np1;
  }
  //x_k+1 = x_k + d_k
  vectorAdd(Nrows, 1.0, d, 1.0, x);
}

void agmgLevel::smoothJacobi(occa::memory o_r, occa::memory o_x,
                             bool x_is_zero) {

  // occaTimerTic(parAlmond->device,"device smoothJacobi");
  if(x_is_zero){
    vectorDotStar(Nrows, 1.0, o_A->o_diagInv, o_r, 0.0, o_x);
    // occaTimerToc(parAlmond->device,"device smoothJacobi");
    return;
  }

  occa::memory o_res = o_scratch;

  // res = r-A*x
  o_A->SpMV(-1.0, o_x, 1.0, o_r, o_res);

  // x = x + alpha*inv(D)*res
  vectorDotStar(Nrows, 1.0, o_A->o_diagInv, o_res, 1.0, o_x);
  // occaTimerToc(parAlmond->device,"hyb smoothJacobi");
}

void agmgLevel::smoothDampedJacobi(occa::memory o_r, occa::memory o_x,
                                   bool x_is_zero){

  // occaTimerTic(parAlmond->device,"device smoothDampedJacobi");
  if(x_is_zero){
    vectorDotStar(Nrows, lambda, o_A->o_diagInv, o_r, 0.0, o_x);
    // occaTimerToc(parAlmond->device,"device smoothDampedJacobi");
    return;
  }

  occa::memory o_res = o_scratch;

  // res = r-A*x
  o_A->SpMV(-1.0, o_x, 1.0, o_r, o_res);

  // x = x + alpha*inv(D)*res
  vectorDotStar(Nrows, lambda, o_A->o_diagInv, o_res, 1.0, o_x);
  // occaTimerToc(parAlmond->device,"device smoothDampedJacobi");
}

void agmgLevel::smoothChebyshev(occa::memory o_r, occa::memory o_x,
                                bool x_is_zero) {

  const dfloat theta = 0.5*(lambda1+lambda0);
  const dfloat delta = 0.5*(lambda1-lambda0);
  const dfloat invTheta = 1.0/theta;
  const dfloat sigma = theta/delta;
  dfloat rho_n = 1./sigma;
  dfloat rho_np1;

  occa::memory o_res = o_scratch + 0*Ncols*sizeof(dfloat);
  occa::memory o_Ad  = o_scratch + 1*Ncols*sizeof(dfloat);
  occa::memory o_d   = o_scratch + 2*Ncols*sizeof(dfloat);

  // occaTimerTic(parAlmond->device,"device smoothChebyshev");

  if(x_is_zero){ //skip the Ax if x is zero
    //res = D^{-1}r
    vectorDotStar(Nrows, 1.0, o_A->o_diagInv, o_r, 0.0, o_res);
    vectorSet(Nrows, 0.0, o_x);
    //d = invTheta*res
    vectorAdd(Nrows, invTheta, o_res, 0.0, o_d);
  } else {
    //res = D^{-1}(r-Ax)
    o_A->SpMV(-1.0, o_x, 1.0, o_r, o_res);
    vectorDotStar(Nrows, o_A->o_diagInv, o_res);

    //d = invTheta*res
    vectorAdd(Nrows, invTheta, o_res, 0.0, o_d);
  }

  for (int k=0;k<ChebyshevIterations;k++) {
    //x_k+1 = x_k + d_k
    vectorAdd(Nrows, 1.0, o_d, 1.0, o_x);

    //r_k+1 = r_k - D^{-1}Ad_k
    o_A->SpMV(1.0, o_d, 0.0, o_Ad);
    vectorDotStar(Nrows, -1.0, o_A->o_diagInv, o_Ad, 1.0, o_res);

    rho_np1 = 1.0/(2.*sigma-rho_n);

    //d_k+1 = rho_k+1*rho_k*d_k  + 2*rho_k+1*r_k+1/delta
    vectorAdd(Nrows, 2.0*rho_np1/delta, o_res, rho_np1*rho_n, o_d);
    rho_n = rho_np1;
  }
  //x_k+1 = x_k + d_k
  vectorAdd(Nrows, 1.0, o_d, 1.0, o_x);

  // occaTimerToc(parAlmond->device,"device smoothChebyshev");
}

} //namespace parAlmond

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
#include "timer.hpp"

int ellipticSolve(elliptic_t* elliptic,
                  occa::memory &o_r, occa::memory &o_x)
{
  mesh_t* mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  int maxIter = 1000;
  options.getArgs("MAXIMUM ITERATIONS", maxIter);
  dfloat tol = 1e-6;
  options.getArgs("SOLVER TOLERANCE", tol);

  elliptic->resNormFactor = 1 / (elliptic->Nfields * mesh->volume);

  if(elliptic->var_coeff && options.compareArgs("PRECONDITIONER", "JACOBI"))
    ellipticUpdateJacobi(elliptic);

  // compute initial residual r = rhs - x0
  ellipticAx(elliptic, mesh->NglobalGatherElements, mesh->o_globalGatherElementList, o_x, elliptic->o_Ap, dfloatString);
  ellipticAx(elliptic, mesh->NlocalGatherElements, mesh->o_localGatherElementList, o_x, elliptic->o_Ap, dfloatString);
  ellipticScaledAdd(elliptic, -1.f, elliptic->o_Ap, 1.f, o_r);
  if(elliptic->allNeumann) ellipticZeroMean(elliptic, o_r);
  oogs::startFinish(o_r, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, ogsAdd, elliptic->oogs);
  if(elliptic->Nmasked) mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, o_r);

  if(options.compareArgs("RESIDUAL PROJECTION","TRUE")) {
    timer::tic("pre",1);
    elliptic->o_x0.copyFrom(o_x, elliptic->Nfields * elliptic->Ntotal * sizeof(dfloat));
    elliptic->residualProjection->pre(o_r);
    timer::toc("pre");
  }

  dlong Niter;
  if(!options.compareArgs("KRYLOV SOLVER", "NONBLOCKING")) {
    Niter = pcg (elliptic, o_r, o_x, tol, maxIter);
  }else{
    printf("NONBLOCKING Krylov solvers currently not supported!");
    ABORT(EXIT_FAILURE);
/*
    if(!options.compareArgs("KRYLOV SOLVER", "FLEXIBLE"))
      Niter = nbpcg (elliptic, o_r, o_x, tol, maxIter);
    else
      Niter = nbfpcg (elliptic, o_r, o_x, tol, maxIter);
 */
  }

  if(options.compareArgs("RESIDUAL PROJECTION","TRUE")) {
    ellipticScaledAdd(elliptic, -1.f, elliptic->o_x0, 1.f, o_x);
    timer::tic("post",1);
    elliptic->residualProjection->post(o_x);
    timer::toc("post");
    ellipticScaledAdd(elliptic, 1.f, elliptic->o_x0, 1.f, o_x);
  }

  if(elliptic->allNeumann)
    ellipticZeroMean(elliptic, o_x);

  return Niter;
}

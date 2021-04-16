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
#include "platform.hpp"
#include "timer.hpp"
#include "linAlg.hpp"

void ellipticSolve(elliptic_t* elliptic, occa::memory &o_r, occa::memory &o_x)
{
  mesh_t* mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  int maxIter = 999;
  options.getArgs("MAXIMUM ITERATIONS", maxIter);
  const int verbose = options.compareArgs("VERBOSE", "TRUE");
  elliptic->resNormFactor = 1 / (elliptic->Nfields * mesh->volume);

  if(verbose) {
    const dfloat rhsNorm = 
      platform->linAlg->weightedNorm2Many(
        mesh->Nlocal,
        elliptic->Nfields,
        elliptic->Ntotal,
        elliptic->o_invDegree,
        o_r,
        platform->comm.mpiComm
      )
      * sqrt(elliptic->resNormFactor); 
    if(platform->comm.mpiRank == 0) printf("RHS norm: %.15e\n", rhsNorm);
  }

  if(options.compareArgs("KRYLOV SOLVER", "PGMRES")){
    elliptic->o_rtmp.copyFrom(o_r, elliptic->Nfields * elliptic->Ntotal * sizeof(dfloat));
    if(elliptic->allNeumann) ellipticZeroMean(elliptic, elliptic->o_rtmp);
    oogs::startFinish(elliptic->o_rtmp, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, ogsAdd, elliptic->oogs);
    if(elliptic->Nmasked) mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, elliptic->o_rtmp);
  }

  if(elliptic->var_coeff && options.compareArgs("PRECONDITIONER", "JACOBI"))
    ellipticUpdateJacobi(elliptic);

  // compute initial residual r = rhs - Ax0
  ellipticAx(elliptic, mesh->NglobalGatherElements, mesh->o_globalGatherElementList, o_x, elliptic->o_Ap, dfloatString);
  ellipticAx(elliptic, mesh->NlocalGatherElements, mesh->o_localGatherElementList, o_x, elliptic->o_Ap, dfloatString);
  platform->linAlg->axpbyMany(
    mesh->Nlocal,
    elliptic->Nfields,
    elliptic->Ntotal,
    -1.0,
    elliptic->o_Ap,
    1.0,
    o_r
  );
  if(elliptic->allNeumann) ellipticZeroMean(elliptic, o_r);
  oogs::startFinish(o_r, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, ogsAdd, elliptic->oogs);
  if(elliptic->Nmasked) mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, o_r);

  if(options.compareArgs("RESIDUAL PROJECTION","TRUE")) {
    platform->timer.tic(elliptic->name + " proj pre",1);
    elliptic->o_x0.copyFrom(o_x, elliptic->Nfields * elliptic->Ntotal * sizeof(dfloat));
    elliptic->res00Norm = 
      platform->linAlg->weightedNorm2Many(
        mesh->Nlocal,
        elliptic->Nfields,
        elliptic->Ntotal,
        elliptic->o_invDegree,
        o_r,
        platform->comm.mpiComm
      )
      * sqrt(elliptic->resNormFactor); 
    if(std::isnan(elliptic->res00Norm)) {
      if(platform->comm.mpiRank == 0) printf("Unreasonable res00Norm!\n");
      ABORT(EXIT_FAILURE);
    }
    elliptic->residualProjection->pre(o_r);
    platform->timer.toc(elliptic->name + " proj pre");
  }

  elliptic->res0Norm = 
    platform->linAlg->weightedNorm2Many(
      mesh->Nlocal,
      elliptic->Nfields,
      elliptic->Ntotal /* offset */,
      elliptic->o_invDegree,
      o_r,
      platform->comm.mpiComm
    )
    * sqrt(elliptic->resNormFactor); 
  if(std::isnan(elliptic->res0Norm)) {
    if(platform->comm.mpiRank == 0) printf("Unreasonable res0Norm!\n");
    ABORT(EXIT_FAILURE);
  }

  dfloat tol = 1e-6;
  options.getArgs("SOLVER TOLERANCE", tol);
  if(options.compareArgs("LINEAR SOLVER STOPPING CRITERION", "RELATIVE")) 
    tol *= elliptic->res0Norm;

  if(!options.compareArgs("KRYLOV SOLVER", "NONBLOCKING")) {
    elliptic->resNorm = elliptic->res0Norm;
    if(options.compareArgs("KRYLOV SOLVER", "PCG"))
      elliptic->Niter = pcg (elliptic, o_r, o_x, tol, maxIter, elliptic->resNorm);
    else if(options.compareArgs("KRYLOV SOLVER", "PGMRES"))
      elliptic->Niter = pgmres (elliptic, o_r, o_x, tol, maxIter, elliptic->resNorm);
    else{
      if(platform->comm.mpiRank == 0) printf("Linear solver %s is not supported!\n", options.getArgs("KRYLOV SOLVER").c_str());
      ABORT(EXIT_FAILURE);
    }
  }else{
    if(platform->comm.mpiRank == 0) printf("NONBLOCKING Krylov solvers currently not supported!");
    ABORT(EXIT_FAILURE);
  }


  if(options.compareArgs("RESIDUAL PROJECTION","TRUE")) { 
    platform->linAlg->axpbyMany(
      mesh->Nlocal,
      elliptic->Nfields,
      elliptic->Ntotal,
      -1.0,
      elliptic->o_x0,
      1.0,
      o_x
    );
    platform->timer.tic(elliptic->name + " proj post",1);
    elliptic->residualProjection->post(o_x);
    platform->timer.toc(elliptic->name + " proj post");
    platform->linAlg->axpbyMany(
      mesh->Nlocal,
      elliptic->Nfields,
      elliptic->Ntotal,
      1.0,
      elliptic->o_x0,
      1.0,
      o_x
    );
  } else {
    elliptic->res00Norm = elliptic->res0Norm;
  }

  if(elliptic->allNeumann)
    ellipticZeroMean(elliptic, o_x);
}

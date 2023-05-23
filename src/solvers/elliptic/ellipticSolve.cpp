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
#include "ellipticPrecon.h"
#include "platform.hpp"
#include "linAlg.hpp"

void ellipticSolve(elliptic_t* elliptic, occa::memory &o_r, occa::memory &o_x)
{
  ellipticUpdateWorkspace(elliptic); // in case o_wrk has changed

  setupAide& options = elliptic->options;
  precon_t *precon = elliptic->precon;

  mesh_t* mesh = elliptic->mesh;

  std::string name = elliptic->name;
  if(name.find("scalar") != std::string::npos){
    name = "scalar";
  }

  int maxIter = 999;
  options.getArgs("MAXIMUM ITERATIONS", maxIter);
  const int verbose = platform->options.compareArgs("VERBOSE", "TRUE");
  elliptic->resNormFactor = 1 / mesh->volume;

  if(verbose) {
    const dfloat rhsNorm = 
      platform->linAlg->weightedNorm2Many(
        mesh->Nlocal,
        elliptic->Nfields,
        elliptic->fieldOffset,
        elliptic->o_invDegree,
        o_r,
        platform->comm.mpiComm
      )
      * sqrt(elliptic->resNormFactor); 
    if(platform->comm.mpiRank == 0) printf("%s RHS norm: %.15e\n", elliptic->name.c_str(), rhsNorm);
  }

  if(verbose) {
    const dfloat rhsNorm = 
      platform->linAlg->weightedNorm2Many(
        mesh->Nlocal,
        elliptic->Nfields,
        elliptic->fieldOffset,
        elliptic->o_invDegree,
        o_x,
        platform->comm.mpiComm
      )
      * sqrt(elliptic->resNormFactor); 
    if(platform->comm.mpiRank == 0) printf("%s x0 norm: %.15e\n", elliptic->name.c_str(), rhsNorm);
  }

  if(options.compareArgs("ELLIPTIC PRECO COEFF FIELD", "TRUE")) {
    if(options.compareArgs("PRECONDITIONER", "MULTIGRID")) {
      ellipticMultiGridUpdateLambda(elliptic);
    }

    if(options.compareArgs("PRECONDITIONER", "JACOBI") ||
       options.compareArgs("MULTIGRID SMOOTHER","DAMPEDJACOBI")) {
      ellipticUpdateJacobi(elliptic);
    }
  }

  // compute initial residual r = rhs - Ax0
  ellipticAx(elliptic, mesh->Nelements, mesh->o_elementList, o_x, elliptic->o_Ap, dfloatString);
  platform->linAlg->axpbyMany(
    mesh->Nlocal,
    elliptic->Nfields,
    elliptic->fieldOffset,
    -1.0,
    elliptic->o_Ap,
    1.0,
    o_r
  );
  if(elliptic->allNeumann) ellipticZeroMean(elliptic, o_r);
  ellipticApplyMask(elliptic, o_r, dfloatString);
  oogs::startFinish(o_r, elliptic->Nfields, elliptic->fieldOffset, ogsDfloat, ogsAdd, elliptic->oogs);

  elliptic->o_x0.copyFrom(o_x, elliptic->Nfields * elliptic->fieldOffset * sizeof(dfloat));
  platform->linAlg->fill(elliptic->fieldOffset * elliptic->Nfields, 0.0, o_x);
  if(options.compareArgs("INITIAL GUESS","PROJECTION") ||
     options.compareArgs("INITIAL GUESS","PROJECTION-ACONJ")) {
    
    platform->timer.tic(name + " proj pre",1);
    elliptic->res00Norm = 
      platform->linAlg->weightedNorm2Many(
        mesh->Nlocal,
        elliptic->Nfields,
        elliptic->fieldOffset,
        elliptic->o_invDegree,
        o_r,
        platform->comm.mpiComm
      )
      * sqrt(elliptic->resNormFactor); 

    nrsCheck(std::isnan(elliptic->res00Norm), MPI_COMM_SELF, EXIT_FAILURE,
             "%s unreasonable res00Norm!\n", name.c_str());

    elliptic->solutionProjection->pre(o_r);

    platform->timer.toc(name + " proj pre");
  }

  elliptic->res0Norm = 
    platform->linAlg->weightedNorm2Many(
      mesh->Nlocal,
      elliptic->Nfields,
      elliptic->fieldOffset /* offset */,
      elliptic->o_invDegree,
      o_r,
      platform->comm.mpiComm
    )
    * sqrt(elliptic->resNormFactor); 

  nrsCheck(std::isnan(elliptic->res0Norm), MPI_COMM_SELF, EXIT_FAILURE,
           "%s unreasonable res00Norm!\n", name.c_str());

  dfloat tol = 1e-6;
  options.getArgs("SOLVER TOLERANCE", tol);
  if(options.compareArgs("LINEAR SOLVER STOPPING CRITERION", "RELATIVE")) 
    tol *= elliptic->res0Norm;

  if(!options.compareArgs("SOLVER", "NONBLOCKING")) {
    elliptic->resNorm = elliptic->res0Norm;

    if(options.compareArgs("SOLVER", "PCG")) {
      elliptic->Niter = pcg (elliptic, o_r, o_x, tol, maxIter, elliptic->resNorm);
    } else if(options.compareArgs("SOLVER", "PGMRES")) {
      elliptic->Niter = pgmres (elliptic, o_r, o_x, tol, maxIter, elliptic->resNorm);
    } else{
      nrsAbort(platform->comm.mpiComm, EXIT_FAILURE,
               "Linear solver %s is not supported!\n", options.getArgs("SOLVER").c_str());
    }

    if(elliptic->Niter == maxIter && platform->comm.mpiRank == 0)
      printf("iteration limit of %s linear solver reached!\n", name.c_str());

  }else{
    nrsAbort(platform->comm.mpiComm, EXIT_FAILURE,
             "%s\n", "NONBLOCKING Krylov solvers currently not supported!");
  }

  if(options.compareArgs("INITIAL GUESS","PROJECTION") ||
     options.compareArgs("INITIAL GUESS","PROJECTION-ACONJ")) { 
    platform->timer.tic(name + " proj post",1);
    elliptic->solutionProjection->post(o_x);
    platform->timer.toc(name + " proj post");
  } else {
    elliptic->res00Norm = elliptic->res0Norm;
  }

  platform->linAlg->axpbyMany(
    mesh->Nlocal,
    elliptic->Nfields,
    elliptic->fieldOffset,
    1.0,
    elliptic->o_x0,
    1.0,
    o_x
  );

  if(elliptic->allNeumann)
    ellipticZeroMean(elliptic, o_x);
}

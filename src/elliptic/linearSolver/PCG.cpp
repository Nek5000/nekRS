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
#include "linAlg.hpp"

int pcg(elliptic_t* elliptic, occa::memory &o_r, occa::memory &o_x,
        const dfloat tol, const int MAXIT, dfloat &rdotr)
{
  
  mesh_t* mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  const int flexible = options.compareArgs("KRYLOV SOLVER", "FLEXIBLE");
  const int verbose = options.compareArgs("VERBOSE", "TRUE");
  const int fixedIterationCountFlag = options.compareArgs("FIXED ITERATION COUNT", "TRUE");

  dfloat rdotz1;
  dfloat alpha;

  /*aux variables */
  occa::memory &o_p  = elliptic->o_p;
  occa::memory &o_z  = elliptic->o_z;
  occa::memory &o_Ap = elliptic->o_Ap;
  occa::memory &o_weight = elliptic->o_invDegree;
  platform->linAlg->fill(elliptic->Nfields * elliptic->Ntotal, 0.0, o_p);

  if(platform->comm.mpiRank == 0 && verbose)
    printf("CG: initial res norm %.15e WE NEED TO GET TO %e \n", rdotr, tol);

  int iter; 
  for(iter = 1; iter <= MAXIT; ++iter) {
    ellipticPreconditioner(elliptic, o_r, o_z);

    const dfloat rdotz2 = rdotz1;
    rdotz1 = platform->linAlg->weightedInnerProdMany(
      mesh->Nlocal,
      elliptic->Nfields,
      elliptic->Ntotal,
      o_weight,
      o_r,
      o_z,
      platform->comm.mpiComm);

    //printf("norm rdotz1: %.15e\n", rdotz1);

    dfloat beta = 0;
    if(iter > 1) {
      beta = rdotz1/rdotz2;
      if(flexible) {
        const dfloat zdotAp = platform->linAlg->weightedInnerProdMany(
          mesh->Nlocal,
          elliptic->Nfields,
          elliptic->Ntotal,
          o_weight,
          o_z,
          o_Ap,
          platform->comm.mpiComm);
        beta = -alpha * zdotAp/rdotz2;
        //printf("norm zdotAp: %.15e\n", zdotAp);
      }
    }

    platform->linAlg->axpbyMany(
      mesh->Nlocal,
      elliptic->Nfields,
      elliptic->Ntotal,
      1.0,
      o_z,
      beta,
      o_p);

    ellipticOperator(elliptic, o_p, o_Ap, dfloatString);
    const dfloat pAp = platform->linAlg->weightedInnerProdMany(
      mesh->Nlocal,
      elliptic->Nfields,
      elliptic->Ntotal,
      o_weight,
      o_p,
      o_Ap,
      platform->comm.mpiComm);
    alpha = rdotz1 / pAp;

    //printf("norm pAp: %.15e\n", pAp);

    //  x <= x + alpha*p
    //  r <= r - alpha*A*p
    //  dot(r,r)
    rdotr = sqrt(ellipticUpdatePCG(elliptic, o_p, o_Ap, alpha, o_x, o_r) * elliptic->resNormFactor);
      
    if (verbose && (platform->comm.mpiRank == 0))
      printf("it %d r norm %.15e\n", iter, rdotr);

    if(rdotr <= tol && !fixedIterationCountFlag) break;
  }

  return iter;
}

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

//#define DEBUG

static dfloat update(elliptic_t* elliptic,
                     occa::memory &o_p, occa::memory &o_Ap, const dfloat alpha,
                     occa::memory &o_x, occa::memory &o_r)
{
  mesh_t* mesh = elliptic->mesh;

  const bool serial = platform->serial;

  // r <= r - alpha*A*p
  // dot(r,r)
  elliptic->updatePCGKernel(mesh->Nlocal,
                            elliptic->fieldOffset,
                            elliptic->o_invDegree,
                            o_Ap,
                            alpha,
                            o_r,
                            elliptic->o_tmpNormr);

  dfloat rdotr1 = 0;
#ifdef ELLIPTIC_ENABLE_TIMER
    //platform->timer.tic("dotp",1);
#endif
  if(serial) {
    rdotr1 = *((dfloat *) elliptic->o_tmpNormr.ptr());
  } else {
    const dlong Nblock = (mesh->Nlocal + BLOCKSIZE - 1) / BLOCKSIZE;
    elliptic->o_tmpNormr.copyTo(elliptic->tmpNormr, Nblock*sizeof(dfloat));
    for(int n = 0; n < Nblock; ++n)
      rdotr1 += elliptic->tmpNormr[n];
  }

  // x <= x + alpha*p
  platform->linAlg->axpbyMany(
    mesh->Nlocal,
    elliptic->Nfields,
    elliptic->fieldOffset,
    alpha,
    o_p,
    1.0,
    o_x);

  MPI_Allreduce(MPI_IN_PLACE, &rdotr1, 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);
#ifdef ELLIPTIC_ENABLE_TIMER
    //platform->timer.toc("dotp");
#endif

  platform->flopCounter->add(elliptic->name + " ellipticUpdatePC",
                             elliptic->Nfields * static_cast<double>(mesh->Nlocal) * 6 + mesh->Nlocal);

  return rdotr1;
}

int pcg(elliptic_t* elliptic, occa::memory &o_r, occa::memory &o_x,
        const dfloat tol, const int MAXIT, dfloat &rdotr)
{
  
  mesh_t* mesh = elliptic->mesh;
  setupAide& options = elliptic->options;

  const int flexible = options.compareArgs("SOLVER", "FLEXIBLE");
  const int verbose = platform->options.compareArgs("VERBOSE", "TRUE");
  const int fixedIteration = false;

  dfloat rdotz1;
  dfloat alpha;

  /*aux variables */
  occa::memory &o_p  = elliptic->o_p;
  occa::memory &o_z = (!options.compareArgs("PRECONDITIONER", "NONE")) ? elliptic->o_z : o_r;
  occa::memory &o_Ap = elliptic->o_Ap;
  occa::memory &o_weight = elliptic->o_invDegree;
  platform->linAlg->fill(elliptic->Nfields * elliptic->fieldOffset, 0.0, o_p);

  if(platform->comm.mpiRank == 0 && verbose) {
    if(flexible) 
      printf("PFCG ");	  
    else
      printf("PCG ");	  
    printf("%s: initial res norm %.15e WE NEED TO GET TO %e \n", elliptic->name.c_str(), rdotr, tol);
  }

  int iter = 0;
  do {
    iter++;
    const dfloat rdotz2 = rdotz1;
    if(!options.compareArgs("PRECONDITIONER", "NONE")) {
      ellipticPreconditioner(elliptic, o_r, o_z);

      rdotz1 = platform->linAlg->weightedInnerProdMany(
        mesh->Nlocal,
        elliptic->Nfields,
        elliptic->fieldOffset,
        o_weight,
        o_r,
        o_z,
        platform->comm.mpiComm);
    } else {
      rdotz1 = rdotr; 
    }

#ifdef DEBUG
    printf("norm rdotz1: %.15e\n", rdotz1);
#endif

    dfloat beta = 0;
    if(iter > 1) {
      beta = rdotz1/rdotz2;
      if(flexible) {
        const dfloat zdotAp = platform->linAlg->weightedInnerProdMany(
          mesh->Nlocal,
          elliptic->Nfields,
          elliptic->fieldOffset,
          o_weight,
          o_z,
          o_Ap,
          platform->comm.mpiComm);
        beta = -alpha * zdotAp/rdotz2;
#ifdef DEBUG
        printf("norm zdotAp: %.15e\n", zdotAp);
#endif
      }
    }

#ifdef DEBUG
        printf("beta: %.15e\n", beta);
#endif


    platform->linAlg->axpbyMany(
      mesh->Nlocal,
      elliptic->Nfields,
      elliptic->fieldOffset,
      1.0,
      o_z,
      beta,
      o_p);

    ellipticOperator(elliptic, o_p, o_Ap, dfloatString);
    const dfloat pAp = platform->linAlg->weightedInnerProdMany(
      mesh->Nlocal,
      elliptic->Nfields,
      elliptic->fieldOffset,
      o_weight,
      o_p,
      o_Ap,
      platform->comm.mpiComm);
    alpha = rdotz1 / (pAp + 1e-300);

#ifdef DEBUG
    printf("alpha: %.15e\n", alpha);
    printf("norm pAp: %.15e\n", pAp);
#endif

    //  x <= x + alpha*p
    //  r <= r - alpha*A*p
    //  dot(r,r)
    rdotr = sqrt(update(elliptic, o_p, o_Ap, alpha, o_x, o_r) * elliptic->resNormFactor);
#ifdef DEBUG
    printf("rdotr: %.15e\n", rdotr);
#endif
    if (platform->comm.mpiRank == 0)
      nrsCheck(std::isnan(rdotr), MPI_COMM_SELF, EXIT_FAILURE,
               "%s\n", "Detected invalid resiual norm while running linear solver!");

    if (verbose && (platform->comm.mpiRank == 0))
      printf("it %d r norm %.15e\n", iter, rdotr);
  }
  while (rdotr > tol && iter < MAXIT);

  return iter;
}

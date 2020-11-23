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

dfloat ellipticUpdatePCG(elliptic_t* elliptic,
                         occa::memory &o_p, occa::memory &o_Ap, const dfloat alpha,
                         occa::memory &o_x, occa::memory &o_r)
{
  setupAide &options = elliptic->options;

  int fixedIterationCountFlag = 0;
  int enableGatherScatters = 1;
  int enableReductions = 1;
  int flexible = options.compareArgs("KRYLOV SOLVER", "FLEXIBLE");
  int verbose = options.compareArgs("VERBOSE", "TRUE");
  int serial = options.compareArgs("THREAD MODEL", "SERIAL");
  int continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");
  int ipdg = options.compareArgs("DISCRETIZATION", "IPDG");

  mesh_t* mesh = elliptic->mesh;
  const dlong Nlocal = mesh->Np * mesh->Nelements;

  dfloat rdotr1 = 0;

  if(serial == 1 && continuous == 1) {
    elliptic->updatePCGKernel(Nlocal,
                              elliptic->Ntotal,
                              elliptic->o_invDegree,
                              o_p,
                              o_Ap,
                              alpha,
                              o_x,
                              o_r,
                              elliptic->o_tmpNormr);

#ifdef ELLIPTIC_ENABLE_TIMER
  timer::tic("dotp",1);
#endif
    elliptic->o_tmpNormr.copyTo(&rdotr1, sizeof(dfloat));
    dfloat globalrdotr1 = 0;
    if(enableReductions)
      MPI_Allreduce(&rdotr1, &globalrdotr1, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
    else
      globalrdotr1 = 1;
#ifdef ELLIPTIC_ENABLE_TIMER
  timer::toc("dotp");
#endif

    return globalrdotr1;
  }

  if(!continuous) {
    // x <= x + alpha*p
    ellipticScaledAdd(elliptic,  alpha, o_p,  1.f, o_x);

    // [
    // r <= r - alpha*A*p
    ellipticScaledAdd(elliptic, -alpha, o_Ap, 1.f, o_r);

    // dot(r,r)
    if(enableReductions)
      rdotr1 = ellipticWeightedNorm2(elliptic, elliptic->o_invDegree, o_r);
    else
      rdotr1 = 1;
  }else{
    // x <= x + alpha*p
    // r <= r - alpha*A*p
    // dot(r,r)
    if(elliptic->blockSolver)
      elliptic->updatePCGKernel(Nlocal,
                                elliptic->Ntotal,
                                elliptic->NblocksUpdatePCG,
                                elliptic->o_invDegree,
                                o_p,
                                o_Ap,
                                alpha,
                                o_x,
                                o_r,
                                elliptic->o_tmpNormr);

    else
      elliptic->updatePCGKernel(Nlocal,
                                elliptic->NblocksUpdatePCG,
                                elliptic->o_invDegree,
                                o_p,
                                o_Ap,
                                alpha,
                                o_x,
                                o_r,
                                elliptic->o_tmpNormr);

    elliptic->o_tmpNormr.copyTo(elliptic->tmpNormr);

    rdotr1 = 0;
    for(int n = 0; n < elliptic->NblocksUpdatePCG; ++n)
      rdotr1 += elliptic->tmpNormr[n];

    dfloat globalrdotr1 = 0;
    MPI_Allreduce(&rdotr1, &globalrdotr1, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

    rdotr1 = globalrdotr1;
  }

  return rdotr1;
}

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

dfloat ellipticUpdatePCG(elliptic_t* elliptic,
                         occa::memory &o_p, occa::memory &o_Ap, const dfloat alpha,
                         occa::memory &o_x, occa::memory &o_r)
{
  mesh_t* mesh = elliptic->mesh;

  const bool serial = platform->serial;

  // r <= r - alpha*A*p
  // dot(r,r)
  elliptic->updatePCGKernel(mesh->Nlocal,
                            elliptic->Ntotal,
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
    elliptic->Ntotal,
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

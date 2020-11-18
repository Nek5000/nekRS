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

typedef union intorfloat {
  int ier;
  float w;
} ierw_t;

#if 0

dfloat ellipticWeightedNorm2(elliptic_t* elliptic, occa::memory &o_w, occa::memory &o_a)
{
  mesh_t* mesh = elliptic->mesh;
  dfloat* tmp = elliptic->tmp;
  dlong Nblock = elliptic->Nblock;
  dlong Nblock2 = elliptic->Nblock2;
  dlong Ntotal = mesh->Nelements * mesh->Np;

  if(elliptic->options.compareArgs("THREAD MODEL", "SERIAL")) {
    const dfloat* __restrict__ cpu_w = (dfloat*)__builtin_assume_aligned(o_w.ptr(),
                                                                         USE_OCCA_MEM_BYTE_ALIGN);
    const dfloat* __restrict__ cpu_a = (dfloat*)__builtin_assume_aligned(o_a.ptr(),
                                                                         USE_OCCA_MEM_BYTE_ALIGN);

    // w'*(a.a)
    dfloat wa2 = 0;

    const hlong M = mesh->Nelements * mesh->Np;

    for(hlong i = 0; i < M; ++i) {
      const dfloat ai = cpu_a[i];
      wa2 += ai * ai * cpu_w[i];
    }

    dfloat globalwa2 = 0;
    MPI_Allreduce(&wa2, &globalwa2, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

    return globalwa2;
  }

  occa::memory &o_tmp = elliptic->o_tmp;
  occa::memory &o_tmp2 = elliptic->o_tmp2;

  if(elliptic->options.compareArgs("DISCRETIZATION","CONTINUOUS"))
    elliptic->weightedNorm2Kernel(Ntotal, o_w, o_a, o_tmp);
  else
    elliptic->innerProductKernel(Ntotal, o_a, o_a, o_tmp);

  /* add a second sweep if Nblock>Ncutoff */
  dlong Ncutoff = 100;
  dlong Nfinal;

  if(Nblock > Ncutoff) {
    mesh->sumKernel(Nblock, o_tmp, o_tmp2);

    o_tmp2.copyTo(tmp);

    Nfinal = Nblock2;
  }else {
    o_tmp.copyTo(tmp);

    Nfinal = Nblock;
  }

  dfloat wab = 0;
  for(dlong n = 0; n < Nfinal; ++n)
    wab += tmp[n];

  dfloat globalwab = 0;
  MPI_Allreduce(&wab, &globalwab, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

  return globalwab;
}

#endif

dfloat ellipticInnerProduct(elliptic_t* elliptic, occa::memory &o_a, occa::memory &o_b)
{
  mesh_t* mesh = elliptic->mesh;
  dfloat* tmp = elliptic->tmp;
  dlong Nblock = elliptic->Nblock;
  dlong Nlocal = mesh->Nelements * mesh->Np;

  occa::memory &o_tmp = elliptic->o_tmp;

#ifdef ELLIPTIC_ENABLE_TIMER
  timer::tic("dotp",1);
#endif

  if(elliptic->blockSolver)
    elliptic->innerProductKernel(Nlocal, elliptic->Ntotal, o_a, o_b, o_tmp);
  else
    elliptic->innerProductKernel(Nlocal, o_a, o_b, o_tmp);

  o_tmp.copyTo(tmp);

  dfloat ab = 0;
  for(dlong n = 0; n < Nblock; ++n)
    ab += tmp[n];

  dfloat globalab = 0;
  MPI_Allreduce(&ab, &globalab, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

#ifdef ELLIPTIC_ENABLE_TIMER
  timer::toc("dotp");
#endif

  return globalab;
}

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

dfloat ellipticWeightedInnerProduct(elliptic_t* elliptic,
                                    occa::memory &o_w,
                                    occa::memory &o_a,
                                    occa::memory &o_b)
{
#ifdef ELLIPTIC_ENABLE_TIMER
  timer::tic("dotp",1);
#endif
  setupAide &options = elliptic->options;

  const int continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");
  const int serial = options.compareArgs("THREAD MODEL", "SERIAL");
  int enableReductions = 1;
  options.getArgs("DEBUG ENABLE REDUCTIONS", enableReductions);

  mesh_t* mesh = elliptic->mesh;
  dfloat* tmp = elliptic->tmp;
  dlong Nblock = elliptic->Nblock;
  dlong Nblock2 = elliptic->Nblock2;

  occa::memory &o_tmp = elliptic->o_tmp;
  occa::memory &o_tmp2 = elliptic->o_tmp2;

  const dlong Nlocal = mesh->Np * mesh->Nelements;

  if(continuous == 1) {
    if(elliptic->blockSolver)
      elliptic->weightedInnerProduct2Kernel(Nlocal, elliptic->Ntotal, o_w, o_a, o_b, o_tmp);
    else
      elliptic->weightedInnerProduct2Kernel(Nlocal, o_w, o_a, o_b, o_tmp);
  }else {
    elliptic->innerProductKernel(Nlocal, o_a, o_b, o_tmp);
  }

  if(serial == 1 && continuous == 1) {
    dfloat wab;
    o_tmp.copyTo(&wab, sizeof(dfloat));
    dfloat globalwab = 0;
    MPI_Allreduce(&wab, &globalwab, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

    return globalwab;
  }

  /* add a second sweep if Nblock>Ncutoff */
  dlong Ncutoff = 4000;
  dlong Nfinal;
  if(Nblock >= Ncutoff) {
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

#ifdef ELLIPTIC_ENABLE_TIMER
  timer::toc("dotp");
#endif

  return globalwab;
}

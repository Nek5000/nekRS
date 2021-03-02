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
#include "linAlg.hpp"

dfloat ellipticWeightedNorm2(elliptic_t* elliptic, occa::memory &o_w, occa::memory &o_a)
{
  mesh_t* mesh = elliptic->mesh;

#ifdef ELLIPTIC_ENABLE_TIMER
  platform->timer.tic("dotp",1);
#endif

  const dfloat result = platform->linAlg->weightedNorm2Many(mesh->Nlocal, elliptic->Nfields, elliptic->Ntotal, o_w, o_a, platform->comm.mpiComm);
#ifdef ELLIPTIC_ENABLE_TIMER
  platform->timer.toc("dotp");
#endif

  return result;

}

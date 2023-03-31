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

void ellipticZeroMean(elliptic_t* elliptic, occa::memory &o_q)
{
  mesh_t* mesh = elliptic->mesh;
  const hlong Nglobal = elliptic->NelementsGlobal * mesh->Np; 

  if(elliptic->blockSolver) {
    nrsAbort(platform->comm.mpiComm, EXIT_FAILURE,
             "NULL space handling for Block solver current not supported!\n", "");
  } else {
    dfloat qmeanGlobal = platform->linAlg->sum(mesh->Nlocal, o_q, platform->comm.mpiComm);
    qmeanGlobal /= (dfloat) Nglobal;
    platform->linAlg->add(mesh->Nlocal, -qmeanGlobal, o_q);
  }
}

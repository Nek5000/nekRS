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
#include "ogsInterface.h"

#define USE_WEIGHTED 1

void ellipticZeroMean(elliptic_t* elliptic, occa::memory &o_q)
{
  dfloat qmeanLocal;
  dfloat qmeanGlobal;

  dlong Nblock = elliptic->Nblock;
  dfloat* tmp = elliptic->tmp;
  mesh_t* mesh = elliptic->mesh;

  occa::memory &o_tmp = elliptic->o_tmp;
  const dlong Nlocal =  mesh->Np * mesh->Nelements;

  if(elliptic->blockSolver) {
    // check field by field
    for(int fld = 0; fld < elliptic->Nfields; fld++)
      // check this field for all Neumann
      if(elliptic->allBlockNeumann[fld]) {
        elliptic->innerProductFieldKernel(Nlocal,
                                          fld,
                                          elliptic->Ntotal,
                                          elliptic->o_invDegree,
                                          o_q,
                                          o_tmp);
#ifdef ELLIPTIC_ENABLE_TIMER
        timer::tic("dotp",1);
#endif
        // finish reduction
        o_tmp.copyTo(tmp);
        qmeanLocal = 0;
        for(dlong n = 0; n < Nblock; ++n)
          qmeanLocal += tmp[n];

        // globalize reduction
        MPI_Allreduce(&qmeanLocal, &qmeanGlobal, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
#ifdef ELLIPTIC_ENABLE_TIMER
        timer::toc("dotp");
#endif

        qmeanGlobal *= elliptic->nullProjectBlockWeightGlobal[fld];

        // q[n] = q[n] - qmeanGlobal for field id :fld
        elliptic->addScalarBlockFieldKernel(Nlocal, fld, elliptic->Ntotal,  -qmeanGlobal, o_q);
      }
  }else{
#if USE_WEIGHTED == 1
    elliptic->innerProductKernel(mesh->Nelements * mesh->Np, elliptic->o_invDegree, o_q, o_tmp);
#else
    mesh->sumKernel(mesh->Nelements * mesh->Np, o_q, o_tmp);
#endif

#ifdef ELLIPTIC_ENABLE_TIMER
  timer::tic("dotp",1);
#endif
    o_tmp.copyTo(tmp);

    // finish reduction
    qmeanLocal = 0;
    for(dlong n = 0; n < Nblock; ++n)
      qmeanLocal += tmp[n];

    // globalize reduction
    MPI_Allreduce(&qmeanLocal, &qmeanGlobal, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
#ifdef ELLIPTIC_ENABLE_TIMER
  timer::toc("dotp");
#endif

    // normalize
#if USE_WEIGHTED == 1
    qmeanGlobal *= elliptic->nullProjectWeightGlobal;
#else
    qmeanGlobal /= ((dfloat) elliptic->NelementsGlobal * (dfloat)mesh->Np);
#endif
    // q[n] = q[n] - qmeanGlobal
    mesh->addScalarKernel(mesh->Nelements * mesh->Np, -qmeanGlobal, o_q);
  }
}

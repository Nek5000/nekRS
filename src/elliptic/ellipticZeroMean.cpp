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
  platform_t* platform = platform_t::getInstance();
  dfloat qmeanLocal;
  dfloat qmeanGlobal;

  dlong Nblock = elliptic->Nblock;
  dfloat* tmp = elliptic->tmp;
  mesh_t* mesh = elliptic->mesh;
  linAlg_t* linAlg = linAlg_t::getInstance();

  occa::memory &o_tmp = elliptic->o_tmp;
  const dlong Nlocal =  mesh->Np * mesh->Nelements;

  if(elliptic->blockSolver) {
    // check field by field
    for(int fld = 0; fld < elliptic->Nfields; fld++)
      // check this field for all Neumann
      if(elliptic->allBlockNeumann[fld]) {
#ifdef ELLIPTIC_ENABLE_TIMER
        timer::tic("dotp",1);
#endif
        dfloat qmeanGlobal =
          linAlg->innerProd(Nlocal,
            elliptic->o_invDegree,
            o_q,
            platform->comm,
            fld * elliptic->Ntotal);
#ifdef ELLIPTIC_ENABLE_TIMER
        timer::toc("dotp");
#endif

        qmeanGlobal *= elliptic->nullProjectBlockWeightGlobal[fld];

        // q[n] = q[n] - qmeanGlobal for field id :fld
        //elliptic->addScalarBlockFieldKernel(Nlocal, fld, elliptic->Ntotal,  -qmeanGlobal, o_q);
        linAlg->add(Nlocal, -qmeanGlobal, o_q, fld * elliptic->Ntotal);
      }
  }else{

#ifdef ELLIPTIC_ENABLE_TIMER
    timer::tic("dotp",1);
#endif
  dfloat qmeanGlobal = 0.0;
#if USE_WEIGHTED == 1
    qmeanGlobal = linAlg->innerProd(mesh->Nlocal, elliptic->o_invDegree, o_q, platform->comm);
#else
    qmeanGlobal = linAlg->sum(mesh->Nlocal, o_q, platform->comm);
#endif
#ifdef ELLIPTIC_ENABLE_TIMER
    timer::toc("dotp");
#endif

    // normalize
    qmeanGlobal /= ((dfloat) elliptic->NelementsGlobal * (dfloat)mesh->Np);
    // q[n] = q[n] - qmeanGlobal
    linAlg->add(Nlocal, -qmeanGlobal, o_q);
  }
}

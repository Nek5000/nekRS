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

#include <vector>

GmresData::GmresData(elliptic_t* elliptic)
: restart(
    [&](){
      int _restart;
      elliptic->options.getArgs("PGMRES RESTART", _restart);
      return _restart;
    }()
  ),
  o_V(elliptic->Nfields * elliptic->Ntotal, restart, sizeof(dfloat)),
  o_Z(elliptic->Nfields * elliptic->Ntotal, restart, sizeof(dfloat)),
  o_y(platform->device.malloc(restart, sizeof(dfloat))),
  H((dfloat *) calloc((restart+1)*(restart+1), sizeof(dfloat))),
  sn((dfloat *) calloc(restart, sizeof(dfloat))),
  cs((dfloat *) calloc(restart, sizeof(dfloat))),
  s((dfloat *) calloc(restart+1, sizeof(dfloat))),
  y((dfloat *) calloc(restart, sizeof(dfloat)))
{
}

void initializeGmresData(elliptic_t* elliptic)
{
  GmresData* gmresData = new GmresData(elliptic);
  elliptic->gmresData = gmresData;
}

namespace{
void gmresUpdate(elliptic_t* elliptic,
                 occa::memory o_x,
                 int I){
  const int restart = elliptic->gmresData->restart;
  mesh_t* mesh = elliptic->mesh;
  dfloat* y = elliptic->gmresData->y;
  dfloat* H = elliptic->gmresData->H;
  dfloat* s = elliptic->gmresData->s;
  deviceVector_t& o_V = elliptic->gmresData->o_V;
  deviceVector_t& o_Z = elliptic->gmresData->o_Z;
  occa::memory& o_y = elliptic->gmresData->o_y;
  occa::memory& o_z = elliptic->o_z;

  for(int k=I-1; k>=0; --k){
    y[k] = s[k];

    for(int m=k+1; m<I; ++m)
      y[k] -= H[k + m*(restart+1)]*y[m];

    y[k] /= H[k + k*(restart+1)];
  }

  o_y.copyFrom(y, I * sizeof(dfloat));

  elliptic->updatePGMRESSolutionKernel(
    mesh->Nlocal,
    elliptic->Ntotal,
    I,
    o_y,
    o_Z,
    o_x
  );
}
}

// Ax=r
int pgmres(elliptic_t* elliptic, occa::memory &o_r, occa::memory &o_x,
        const dfloat tol, const int MAXIT, dfloat &rdotr)
{

  mesh_t* mesh = elliptic->mesh;
  linAlg_t& linAlg = *(platform->linAlg);

  occa::memory& o_w = elliptic->o_p;
  occa::memory& o_b = elliptic->o_rtmp;

  occa::memory& o_z = elliptic->o_z;
  occa::memory& o_Ax = elliptic->o_Ap;
  deviceVector_t& o_V = elliptic->gmresData->o_V;
  deviceVector_t& o_Z = elliptic->gmresData->o_Z;

  dfloat* y = elliptic->gmresData->y;
  dfloat* H = elliptic->gmresData->H;
  dfloat* sn = elliptic->gmresData->sn;
  dfloat* cs = elliptic->gmresData->cs;
  dfloat* s = elliptic->gmresData->s;

  const int restart = elliptic->gmresData->restart;

  const bool verbose = platform->options.compareArgs("VERBOSE", "TRUE");

  dfloat nr = rdotr / sqrt(elliptic->resNormFactor);
  dfloat error = rdotr;
  const dfloat TOL = tol;

  if (verbose&&(platform->comm.mpiRank==0))
    printf("PGMRES: initial res norm %12.12f \n", rdotr);

  int iter=0;

  for(iter=0;iter<MAXIT;){

    s[0] = nr;

    // V(:,0) = r/nr
    linAlg.axpbyMany(
      mesh->Nlocal,
      elliptic->Nfields,
      elliptic->Ntotal,
      1.0 / nr,
      o_r,
      0.0,
      o_V);

    //Construct orthonormal basis via Gram-Schmidt
    for(int i=0;i<restart;++i){
      // z := M^{-1} V(:,i)
      ellipticPreconditioner(elliptic, o_V.at(i), o_Z.at(i));

      // w := A z
      ellipticOperator(elliptic, o_Z.at(i), o_w, dfloatString);

      for(int k=0; k<=i; ++k){
        dfloat hki = linAlg.weightedInnerProdMany(
          mesh->Nlocal,
          elliptic->Nfields,
          elliptic->Ntotal,
          elliptic->o_invDegree,
          o_w, o_V.at(k), platform->comm.mpiComm);

        // w = w - hki*V[k]
        linAlg.axpbyMany(
          mesh->Nlocal,
          elliptic->Nfields,
          elliptic->Ntotal,
          -hki, o_V.at(k), 1.0, o_w);

        // H(k,i) = hki
        H[k + i*(restart+1)] = hki;
      }

      dfloat nw = linAlg.weightedNorm2Many(
        mesh->Nlocal,
        elliptic->Nfields,
        elliptic->Ntotal,
        elliptic->o_invDegree, o_w, platform->comm.mpiComm);

      // H(i+1,i) = ||w||_2
      H[i+1 + i*(restart+1)] = nw;

      // V(:,i+1) = w/nw
      if (i<restart-1)
        linAlg.axpbyMany(
          mesh->Nlocal,
          elliptic->Nfields,
          elliptic->Ntotal,
          (1./nw), o_w, 0., o_V.at(i+1));

      //apply Givens rotation
      for(int k=0; k<i; ++k){
        const dfloat h1 = H[k +     i*(restart+1)];
        const dfloat h2 = H[k + 1 + i*(restart+1)];

        H[k +     i*(restart+1)] =  cs[k]*h1 + sn[k]*h2;
        H[k + 1 + i*(restart+1)] = -sn[k]*h1 + cs[k]*h2;
      }

      // form i-th rotation matrix
      const dfloat h1 = H[i+    i*(restart+1)];
      const dfloat h2 = H[i+1 + i*(restart+1)];
      const dfloat hr = sqrt(h1*h1 + h2*h2);
      cs[i] = h1/hr;
      sn[i] = h2/hr;

      H[i   + i*(restart+1)] = cs[i]*h1 + sn[i]*h2;
      H[i+1 + i*(restart+1)] = 0;

      //approximate residual norm
      s[i+1] = -sn[i]*s[i];
      s[i]   =  cs[i]*s[i];

      iter++;
      error = fabs(s[i+1]) * sqrt(elliptic->resNormFactor);
      rdotr = error;

      if (verbose&&(platform->comm.mpiRank==0)) {
        printf("GMRES: it %d, approx residual norm %12.12le \n", iter, error);
      }

      if(error < TOL || iter==MAXIT) {
        //update approximation
        gmresUpdate(elliptic, o_x, i+1);
        break;
      }
    }

    //exit if tolerance is reached
    if(error < TOL || iter==MAXIT) break;

    //update approximation
    gmresUpdate(elliptic, o_x, restart);

    // restart GMRES
    // compute A*x
    ellipticOperator(elliptic, o_x, o_Ax, dfloatString);

    // subtract r = b - A*x
    linAlg.axpbyzMany(
      mesh->Nlocal,
      elliptic->Nfields,
      elliptic->Ntotal,
      -1.0, o_Ax, 1.0, o_b, o_r);

    nr = linAlg.weightedNorm2Many(
      mesh->Nlocal,
      elliptic->Nfields,
      elliptic->Ntotal,
      elliptic->o_invDegree, o_r, platform->comm.mpiComm);

    error = nr * sqrt(elliptic->resNormFactor);
    rdotr = nr * sqrt(elliptic->resNormFactor);
    //exit if tolerance is reached
    if(error<=TOL) return iter;
  }

  return iter;
}
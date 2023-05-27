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

#include <stdio.h>
#include <stdlib.h>
#include "mesh3D.h"
#include "platform.hpp"
#include "linAlg.hpp"

void mesh_t::geometricFactors()
{
  const auto tmpSize = Nlocal * sizeof(dfloat);
  occa::memory o_tmp = (platform->o_mempool.o_ptr.size() >= tmpSize) ?
                       platform->o_mempool.slice0 : 
                       platform->device.malloc(tmpSize); 

  geometricFactorsKernel(Nelements,
                         o_D,
                         o_gllw,
                         o_x,
                         o_y,
                         o_z,
                         o_LMM,
                         o_vgeo,
                         o_ggeo,
                         o_tmp);


  const dfloat minJ =
      platform->linAlg->min(Nlocal, o_tmp, platform->comm.mpiComm);
  const dfloat maxJ =
      platform->linAlg->max(Nlocal, o_tmp, platform->comm.mpiComm);

  nrsCheck(minJ < 0 || maxJ < 0, platform->comm.mpiComm, EXIT_FAILURE,
           "%s\n", "Invalid element Jacobian < 0 found!");

  double flopsCubatureGeometricFactors = 0.0;
  if (cubNq > 1) {
    cubatureGeometricFactorsKernel(Nelements, o_cubD, o_x, o_y, o_z, o_cubInterpT, o_cubw, o_cubvgeo);

    flopsCubatureGeometricFactors += 18 * Np * Nq;                                             // deriv
    flopsCubatureGeometricFactors += 18 * (cubNq * Np + cubNq * cubNq * Nq * Nq + cubNp * Nq); // c->f interp
    flopsCubatureGeometricFactors += 55 * cubNp; // geometric factor computation
    flopsCubatureGeometricFactors *= static_cast<double>(Nelements);
  }

  double flopsGeometricFactors = 18 * Np * Nq + 91 * Np;
  flopsGeometricFactors *= static_cast<double>(Nelements);
  platform->flopCounter->add("mesh_t::update", flopsGeometricFactors + flopsCubatureGeometricFactors);
}

void interpolateHex3D(dfloat* I, dfloat* x, int N, dfloat* Ix, int M)
{
  dfloat* Ix1 = (dfloat*) calloc(N * N * M, sizeof(dfloat));
  dfloat* Ix2 = (dfloat*) calloc(N * M * M, sizeof(dfloat));

  for(int k = 0; k < N; ++k)
    for(int j = 0; j < N; ++j)
      for(int i = 0; i < M; ++i) {
        dfloat tmp = 0;
        for(int n = 0; n < N; ++n)
          tmp += I[i * N + n] * x[k * N * N + j * N + n];
        Ix1[k * N * M + j * M + i] = tmp;
      }

  for(int k = 0; k < N; ++k)
    for(int j = 0; j < M; ++j)
      for(int i = 0; i < M; ++i) {
        dfloat tmp = 0;
        for(int n = 0; n < N; ++n)
          tmp += I[j * N + n] * Ix1[k * N * M + n * M + i];
        Ix2[k * M * M + j * M + i] = tmp;
      }

  for(int k = 0; k < M; ++k)
    for(int j = 0; j < M; ++j)
      for(int i = 0; i < M; ++i) {
        dfloat tmp = 0;
        for(int n = 0; n < N; ++n)
          tmp += I[k * N + n] * Ix2[n * M * M + j * M + i];
        Ix[k * M * M + j * M + i] = tmp;
      }

  free(Ix1);
  free(Ix2);
}

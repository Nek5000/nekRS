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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "platform.hpp"
#include "mesh3D.h"

void mesh_t::surfaceGeometricFactors()
{
  surfaceGeometricFactorsKernel(Nelements, o_gllw, o_faceNodes, o_vgeo, o_sgeo);

  double flopsSurfaceGeometricFactors = 32 * Nq * Nq;
  flopsSurfaceGeometricFactors *= static_cast<double>(Nelements);

  platform->flopCounter->add("mesh_t::update", flopsSurfaceGeometricFactors);
}


void interpolateFaceHex3D(int* faceNodes, dfloat* I, dfloat* x, int N, dfloat* Ix, int M)
{
  dfloat* Ix0 = (dfloat*) calloc(N * N, sizeof(dfloat));
  dfloat* Ix1 = (dfloat*) calloc(N * M, sizeof(dfloat));

  for(int j = 0; j < N; ++j)
    for(int i = 0; i < N; ++i)
      Ix0[j * N + i] = x[faceNodes[j * N + i]];

  for(int j = 0; j < N; ++j)
    for(int i = 0; i < M; ++i) {
      dfloat tmp = 0;
      for(int n = 0; n < N; ++n)
        tmp += I[i * N + n] * Ix0[j * N + n];
      Ix1[j * M + i] = tmp;
    }

  for(int j = 0; j < M; ++j)
    for(int i = 0; i < M; ++i) {
      dfloat tmp = 0;
      for(int n = 0; n < N; ++n)
        tmp += I[j * N + n] * Ix1[n * M + i];
      Ix[j * M + i] = tmp;
    }

  free(Ix0);
  free(Ix1);
}

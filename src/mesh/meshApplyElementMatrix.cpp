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

#include "mesh.h"

void meshApplyElementMatrix(mesh_t* mesh, dfloat* A, dfloat* q, dfloat* Aq)
{
  dfloat* Aqn = (dfloat*) calloc(mesh->Np,sizeof(dfloat));
  for (dlong e = 0; e < mesh->Nelements; e++) {
    for (int n = 0; n < mesh->Np; n++) {
      Aqn[n] = 0;
      for (int k = 0; k < mesh->Np; k++)
        Aqn[n] += A[k + n * mesh->Np] * q[k + e * mesh->Np];
    }
    for (int n = 0; n < mesh->Np; n++) Aq[n + e * mesh->Np] = Aqn[n];
  }
  free(Aqn);
}

void meshApplyVectorElementMatrix(mesh_t* mesh,
                                  int Nfield,
                                  const dlong offset,
                                  dfloat* A,
                                  dfloat* q,
                                  dfloat* Aq)
{
  dfloat* Aqn = (dfloat*) calloc(mesh->Np,sizeof(dfloat));
  for(int fld = 0; fld < Nfield; fld++)
    for (dlong e = 0; e < mesh->Nelements; e++) {
      for (int n = 0; n < mesh->Np; n++) {
        Aqn[n] = 0;
        for (int k = 0; k < mesh->Np; k++)
          Aqn[n] += A[k + n * mesh->Np] * q[k + e * mesh->Np + fld * offset];
      }
      for (int n = 0; n < mesh->Np; n++) Aq[n + e * mesh->Np + fld * offset] = Aqn[n];
    }
  free(Aqn);
}
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
#define NODE_GEN

void meshLoadReferenceNodesHex3D(mesh3D* mesh, int N, int cubN)
{
  mesh->N = N;
  mesh->Nq = N + 1;
  mesh->cubNq = cubN + 1;
  mesh->Nfp = (N + 1) * (N + 1);
  mesh->Np = (N + 1) * (N + 1) * (N + 1);
  mesh->Nverts = 8;

  int Nrows, Ncols;

  mesh->r = (dfloat*) malloc(mesh->Np * sizeof(dfloat));
  mesh->s = (dfloat*) malloc(mesh->Np * sizeof(dfloat));
  mesh->t = (dfloat*) malloc(mesh->Np * sizeof(dfloat));
  NodesHex3D(mesh->N, mesh->r, mesh->s, mesh->t);

  mesh->faceNodes = (int*) malloc(mesh->Nfaces * mesh->Nfp * sizeof(int));
  FaceNodesHex3D(mesh->N, mesh->r, mesh->s, mesh->t, mesh->faceNodes);

  //GLL quadrature
  mesh->gllz = (dfloat*) malloc((mesh->N + 1) * sizeof(dfloat));
  mesh->gllw = (dfloat*) malloc((mesh->N + 1) * sizeof(dfloat));
  JacobiGLL(mesh->N, mesh->gllz, mesh->gllw);

  mesh->D = (dfloat*) malloc(mesh->Nq * mesh->Nq * sizeof(dfloat));
  Dmatrix1D(mesh->N, mesh->Nq, mesh->gllz, mesh->Nq, mesh->gllz, mesh->D);

  mesh->DW = (dfloat*) malloc(mesh->Nq * mesh->Nq * sizeof(dfloat));
  DWmatrix1D(mesh->N, mesh->D, mesh->DW);

  mesh->interpRaise = (dfloat* ) calloc(mesh->Nq * (mesh->Nq + 1),sizeof(dfloat));
  mesh->interpLower = (dfloat* ) calloc((mesh->Nq - 1) * (mesh->Nq),sizeof(dfloat));
  DegreeRaiseMatrix1D(mesh->N, mesh->N + 1, mesh->interpRaise);
  DegreeRaiseMatrix1D(mesh->N - 1, mesh->N, mesh->interpLower);

  mesh->cubNfp = mesh->cubNq * mesh->cubNq;
  mesh->cubNp = mesh->cubNq * mesh->cubNq * mesh->cubNq;
  // cubN+1 point Gauss-Legendre quadrature
  mesh->cubr = (dfloat*) malloc(mesh->cubNq * sizeof(dfloat));
  mesh->cubw = (dfloat*) malloc(mesh->cubNq * sizeof(dfloat));
  JacobiGLL(mesh->cubNq - 1, mesh->cubr, mesh->cubw);

  mesh->cubInterp = (dfloat*) calloc(mesh->Nq * mesh->cubNq, sizeof(dfloat));
  InterpolationMatrix1D(mesh->N, mesh->Nq, mesh->r, mesh->cubNq, mesh->cubr, mesh->cubInterp); //uses the fact that r = gllz for 1:Nq

  //cubature project cubProject = cubInterp^T
  mesh->cubProject = (dfloat*) calloc(mesh->cubNq * mesh->Nq, sizeof(dfloat));
  matrixTranspose(mesh->cubNq, mesh->Nq, mesh->cubInterp, mesh->Nq, mesh->cubProject, mesh->cubNq);

  //cubature derivates matrix, cubD: differentiate on cubature nodes
  mesh->cubD = (dfloat*) malloc(mesh->cubNq * mesh->cubNq * sizeof(dfloat));
  Dmatrix1D(mesh->cubNq - 1, mesh->cubNq, mesh->cubr, mesh->cubNq, mesh->cubr, mesh->cubD);
  // weak cubature derivative = cubD^T
  mesh->cubDW  = (dfloat*) calloc(mesh->cubNq * mesh->cubNq, sizeof(dfloat));
  for(int i = 0; i < mesh->cubNq; ++i)
    for(int j = 0; j < mesh->cubNq; ++j)
      mesh->cubDW[j + i * mesh->cubNq] = mesh->cubD[i + j * mesh->cubNq];

  mesh->intNfp = 0;
  mesh->intLIFT = NULL;
  mesh->max_EL_nnz = 0;
  mesh->intNfp = 0;

  // find node indices of vertex nodes
  dfloat NODETOL = 1e-6;
  mesh->vertexNodes = (int*) calloc(mesh->Nverts, sizeof(int));
  for(int n = 0; n < mesh->Np; ++n) {
    if( (mesh->r[n] + 1) * (mesh->r[n] + 1) + (mesh->s[n] + 1) * (mesh->s[n] + 1) +
        (mesh->t[n] + 1) * (mesh->t[n] + 1) < NODETOL)
      mesh->vertexNodes[0] = n;
    if( (mesh->r[n] - 1) * (mesh->r[n] - 1) + (mesh->s[n] + 1) * (mesh->s[n] + 1) +
        (mesh->t[n] + 1) * (mesh->t[n] + 1) < NODETOL)
      mesh->vertexNodes[1] = n;
    if( (mesh->r[n] - 1) * (mesh->r[n] - 1) + (mesh->s[n] - 1) * (mesh->s[n] - 1) +
        (mesh->t[n] + 1) * (mesh->t[n] + 1) < NODETOL)
      mesh->vertexNodes[2] = n;
    if( (mesh->r[n] + 1) * (mesh->r[n] + 1) + (mesh->s[n] - 1) * (mesh->s[n] - 1) +
        (mesh->t[n] + 1) * (mesh->t[n] + 1) < NODETOL)
      mesh->vertexNodes[3] = n;
    if( (mesh->r[n] + 1) * (mesh->r[n] + 1) + (mesh->s[n] + 1) * (mesh->s[n] + 1) +
        (mesh->t[n] - 1) * (mesh->t[n] - 1) < NODETOL)
      mesh->vertexNodes[4] = n;
    if( (mesh->r[n] - 1) * (mesh->r[n] - 1) + (mesh->s[n] + 1) * (mesh->s[n] + 1) +
        (mesh->t[n] - 1) * (mesh->t[n] - 1) < NODETOL)
      mesh->vertexNodes[5] = n;
    if( (mesh->r[n] - 1) * (mesh->r[n] - 1) + (mesh->s[n] - 1) * (mesh->s[n] - 1) +
        (mesh->t[n] - 1) * (mesh->t[n] - 1) < NODETOL)
      mesh->vertexNodes[6] = n;
    if( (mesh->r[n] + 1) * (mesh->r[n] + 1) + (mesh->s[n] - 1) * (mesh->s[n] - 1) +
        (mesh->t[n] - 1) * (mesh->t[n] - 1) < NODETOL)
      mesh->vertexNodes[7] = n;
  }

  mesh->max_EL_nnz = 0;
}

/*

   The MIT License (MIT)

   Copyright (c) 2020 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mesh.h"

extern "C" {
void dgesv_ ( int* N, int* NRHS, double* A,
              int* LDA,
              int* IPIV,
              double* B,
              int* LDB,
              int* INFO );
}

// C = A/B  = trans(trans(B)\trans(A))
// assume row major
void matrixRightSolve(int NrowsA, int NcolsA, dfloat* A, int NrowsB, int NcolsB, dfloat* B, dfloat* C)
{
  int info;

  int NrowsX = NcolsB;
  int NcolsX = NrowsB;

  int NrowsY = NcolsA;
  int NcolsY = NrowsA;

  int lwork = NrowsX * NcolsX;

  // compute inverse mass matrix
  double* tmpX = (double*) calloc(NrowsX * NcolsX, sizeof(double));
  double* tmpY = (double*) calloc(NrowsY * NcolsY, sizeof(double));

  int* ipiv = (int*) calloc(NrowsX, sizeof(int));
  double* work = (double*) calloc(lwork, sizeof(double));

  for(int n = 0; n < NrowsX * NcolsX; ++n)
    tmpX[n] = B[n];

  for(int n = 0; n < NrowsY * NcolsY; ++n)
    tmpY[n] = A[n];

  dgesv_(&NrowsX, &NcolsY, tmpX, &NrowsX, ipiv, tmpY, &NrowsY, &info); // ?

  for(int n = 0; n < NrowsY * NcolsY; ++n)
    C[n] = tmpY[n];

  if(info)
    printf("matrixRightSolve: dgesv reports info = %d when inverting matrix\n", info);

  free(work);
  free(ipiv);
  free(tmpX);
  free(tmpY);
}

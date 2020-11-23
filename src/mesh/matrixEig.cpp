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
void dgeev_(char* JOBVL, char* JOBVR, int* N, double* A, int* LDA, double* WR, double* WI,
            double* VL, int* LDVL, double* VR, int* LDVR, double* WORK, int* LWORK, int* INFO );
}

// compute right eigenvectors
void matrixEig(int N, dfloat* A, dfloat* VR, dfloat* WR, dfloat* WI)
{
  char JOBVL = 'N';
  char JOBVR = 'V';
  int LDA = N;
  int LDVL = N;
  int LDVR = N;
  int LWORK = 8 * N;

  double* tmpA  = (double*) calloc(N * N,sizeof(double));
  double* tmpWR = (double*) calloc(N,sizeof(double));
  double* tmpWI = (double*) calloc(N,sizeof(double));
  double* tmpVR = (double*) calloc(N * N,sizeof(double));
  double* tmpVL = NULL;
  double* WORK  = (double*) calloc(LWORK,sizeof(double));

  int info;

  for(int n = 0; n < N; ++n)
    for(int m = 0; m < N; ++m)
      tmpA[n + m * N] = A[n * N + m];

  dgeev_ (&JOBVL, &JOBVR, &N, tmpA, &LDA, tmpWR, tmpWI, tmpVL, &LDVL, tmpVR, &LDVR, WORK, &LWORK, &info);

  for(int n = 0; n < N; ++n) {
    WR[n] = tmpWR[n];
    WI[n] = tmpWI[n];
    for(int m = 0; m < N; ++m)
      VR[n + m * N] = tmpVR[n * N + m];
  }
}

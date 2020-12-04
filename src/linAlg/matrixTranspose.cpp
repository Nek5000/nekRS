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

#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mesh.h>

void matrixTranspose(const int M, const int N,
                     const dfloat* A, const int LDA,
                     dfloat* AT, const int LDAT)
{
  //A & A^T - Row major ordering
  //M = number of rows of A, columns of A^T
  //N = number of columns of A, rows of A^T
  //LDA  - leading dimension of A (>=M)
  //LDAT - leading dimension of A^T (>=N)

  //quick return
  if (N < 1 || M < 1) return;

  //check for weird input
  if (LDA < N || LDAT < M) {
    printf("Bad input to matrixTranspose\n");
    return;
  }

  for (int n = 0; n < N; n++) //for all cols of A^T
    for (int m = 0; m < M; m++) //for all rows of A^T
      AT[n * LDAT + m] = A[m * LDA + n];
}
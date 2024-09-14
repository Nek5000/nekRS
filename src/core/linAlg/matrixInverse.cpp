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

#include "nekrsSys.hpp"
#include "linAlg.hpp"


extern "C" {
  void dgetrf_(int *M, int *N, double *A, int *lda, int *IPIV, int *INFO);
  void dgetri_(int *N, double *A, int *lda, int *IPIV, double *WORK, int *lwork, int *INFO);
}

std::vector<dfloat> linAlg_t::matrixInverse(const int N_, const std::vector<dfloat>& A)
{
  int N = N_;
  int lwork = N * N;
  int info;

  std::vector<double> invA(N * N, 0.0);

  int* ipiv = (int*) calloc(N, sizeof(int));
  double* work = (double*) calloc(lwork, sizeof(double));

  for(int n = 0; n < N * N; ++n)
    invA[n] = A[n];

  dgetrf_ (&N, &N, invA.data(), &N, ipiv, &info);
  dgetri_ (&N, invA.data(), &N, ipiv, work, &lwork, &info);

  if(info)
    printf("inv: dgetrf/dgetri reports info = %d when inverting matrix\n", info);

  free(work);
  free(ipiv);
 
  std::vector<dfloat> out(invA.size());
  for (int i = 0; i < invA.size(); ++i) out[i] = invA[i];
 
  return out;
}

std::vector<dfloat> linAlg_t::matrixPseudoInverse(const int N, const std::vector<dfloat>& A)
{
  const int M = A.size() / N; // rows of A

  // inv(trans(A)*A)
  std::vector<dfloat> V(M * N);
  for(int n=0; n<N; n++){
    for(int m=0; m<N; m++){
      double tmp = 0; 
      for(int i = 0; i<M; i++){
       tmp += A[i*N +n]*A[i*N +m];
      }
      V[n*N + m] = tmp; 
    }
  }

  auto invV = matrixInverse(N, V); 

  // // (A^T*A)^-1 * A^T
  std::vector<dfloat> invA(N * M);
  for(int n=0; n<N; n++){
    for(int m=0; m<M; m++){
      double tmp = 0; 
      for(int i=0; i<N; i++){
        tmp += invV[n*N + i]*A[m*N + i];
      }
      invA[n*M + m] = tmp; 
    }
  }
 
  return invA;
}

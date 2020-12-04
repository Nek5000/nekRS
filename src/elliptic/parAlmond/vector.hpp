/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#ifndef PARALMOND_VECTOR_HPP
#define PARALMOND_VECTOR_HPP

namespace parAlmond {

//------------------------------------------------------------------------
//
//  Host vector operations
//
//------------------------------------------------------------------------

void vectorSet(const dlong m, const dfloat alpha, dfloat *a);

void vectorRandomize(const dlong m, dfloat *a);

void vectorScale(const dlong m, const dfloat alpha, dfloat *a);

void vectorAddScalar(const dlong m, const dfloat alpha, dfloat *a);

// y = beta*y + alpha*x
void vectorAdd(const dlong n, const dfloat alpha, const dfloat *x,
               const dfloat beta, dfloat *y);

// z = beta*y + alpha*x
void vectorAdd(const dlong n, const dfloat alpha, const dfloat *x,
               const dfloat beta, const dfloat *y, dfloat *z);

// b = a*b
void vectorDotStar(const dlong m, const dfloat *a, dfloat *b);

// c = alpha*a*b + beta*c
void vectorDotStar(const dlong m, const dfloat alpha, const dfloat *a,
                   const dfloat *b, const dfloat beta,  dfloat *c);

dfloat vectorNorm(const dlong n, const dfloat *a, MPI_Comm comm);

dfloat vectorInnerProd(const dlong n, const dfloat *a, const dfloat *b,
                       MPI_Comm comm);

dfloat vectorMaxAbs(const dlong n, const dfloat *a, MPI_Comm comm);

// returns aDotbc[0] = a\dot b, aDotbc[1] = a\dot c, aDotbc[2] = b\dot b,
void kcycleCombinedOp1(const dlong n, dfloat *aDotbc, const dfloat *a,
                      const dfloat *b, const dfloat *c, const dfloat* w,
                      const bool weighted, MPI_Comm comm);

// returns aDotbcd[0] = a\dot b, aDotbcd[1] = a\dot c, aDotbcd[2] = a\dot d,
void kcycleCombinedOp2(const dlong n, dfloat *aDotbcd, const dfloat *a,
                       const dfloat *b, const dfloat *c, const dfloat* d,
                       const dfloat *w, const bool weighted, MPI_Comm comm);

// y = beta*y + alpha*x, and return y\dot y
dfloat vectorAddInnerProd(const dlong n, const dfloat alpha, const dfloat *x,
                          const dfloat beta, dfloat *y,
                          const dfloat *w, const bool weighted, MPI_Comm comm);

//------------------------------------------------------------------------
//
//  Device vector operations
//
//------------------------------------------------------------------------

void vectorSet(const dlong N, const dfloat alpha, occa::memory o_a);

void vectorRandomize(const dlong m, occa::memory o_a);

void vectorScale(const dlong N, const dfloat alpha, occa::memory o_a);

void vectorAddScalar(const dlong N, const dfloat alpha, occa::memory o_a);

void vectorAdd(const dlong N, const dfloat alpha, occa::memory o_x,
               const dfloat beta, occa::memory o_y);

void vectorAdd(const dlong N, const dfloat alpha, occa::memory o_x,
               const dfloat beta, occa::memory o_y, occa::memory o_z);

void vectorDotStar(const dlong N, occa::memory o_a, occa::memory o_b);

void vectorDotStar(const dlong N, const dfloat alpha, occa::memory o_a,
                   occa::memory o_b, const dfloat beta, occa::memory o_c);

dfloat vectorNorm(const dlong n, occa::memory o_a, MPI_Comm comm);

dfloat vectorInnerProd(const dlong N, occa::memory o_x, occa::memory o_y,
                       MPI_Comm comm);

dfloat vectorMaxAbs(const dlong n, occa::memory o_a, MPI_Comm comm);

// returns aDotbc[0] = a\dot b, aDotbc[1] = a\dot c, aDotbc[2] = b\dot b,
void kcycleCombinedOp1(const dlong N, dfloat *aDotbc, occa::memory o_a,
                       occa::memory o_b, occa::memory o_c, occa::memory o_w,
                       const bool weighted, MPI_Comm comm);

// returns aDotbcd[0] = a\dot b, aDotbcd[1] = a\dot c, aDotbcd[2] = a\dot d,
void kcycleCombinedOp2(const dlong N, dfloat *aDotbcd,
                        occa::memory o_a, occa::memory o_b,
                        occa::memory o_c, occa::memory o_d,
                        occa::memory o_w, const bool weighted, MPI_Comm comm);

// y = beta*y + alpha*x, and return y\dot y
dfloat vectorAddInnerProd(const dlong N, const dfloat alpha, occa::memory o_x,
                          const dfloat beta, occa::memory o_y,
                          occa::memory o_w, const bool weighted, MPI_Comm comm);

} //namespace parAlmond

#endif
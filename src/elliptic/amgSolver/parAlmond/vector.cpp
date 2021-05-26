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

#include "parAlmond.hpp"

namespace parAlmond {

//------------------------------------------------------------------------
//
//  Host vector operations
//
//------------------------------------------------------------------------

void vectorSet(const dlong m, const dfloat alpha, dfloat *a){
  // #pragma omp parallel for
  for(dlong i=0; i<m; i++)
    a[i] = alpha;
}

void vectorRandomize(const dlong m, dfloat *a){
  // #pragma omp parallel for
  for(dlong i=0; i<m; i++)
    a[i] = (dfloat) drand48();
}

void vectorScale(const dlong m, const dfloat alpha, dfloat *a){
  // #pragma omp parallel for
  for(dlong i=0; i<m; i++)
    a[i] *= alpha;
}

void vectorAddScalar(const dlong m, const dfloat alpha, dfloat *a){
  // #pragma omp parallel for
  for(dlong i=0; i<m; i++)
    a[i] += alpha;
}

// y = beta*y + alpha*x
void vectorAdd(const dlong n, const dfloat alpha, const dfloat *x,
               const dfloat beta, dfloat *y){
  if (beta) {
    // #pragma omp parallel for
    for(dlong i=0; i<n; i++)
      y[i] = beta*y[i] + alpha*x[i];
  } else {
    // #pragma omp parallel for
    for(dlong i=0; i<n; i++)
      y[i] = alpha*x[i];
  }
}

// z = beta*y + alpha*x
void vectorAdd(const dlong n, const dfloat alpha, const dfloat *x,
               const dfloat beta, const dfloat *y, dfloat *z){
  // #pragma omp parallel for
  for(dlong i=0; i<n; i++)
    z[i] = beta*y[i] + alpha*x[i];
}

// b = a*b
void vectorDotStar(const dlong m, const dfloat *a, dfloat *b){
  // #pragma omp parallel for
  for(dlong i=0; i<m; i++)
    b[i] *= a[i];
}

// c = alpha*a*b + beta*c
void vectorDotStar(const dlong m, const dfloat alpha, const dfloat *a,
                   const dfloat *b, const dfloat beta,  dfloat *c){
  if (beta) {
    // #pragma omp parallel for
    for(dlong i=0; i<m; i++)
      c[i] = beta*c[i]+ alpha*a[i]*b[i];
  } else {
    // #pragma omp parallel for
    for(dlong i=0; i<m; i++)
      c[i] = alpha*a[i]*b[i];
  }
}

dfloat vectorNorm(const dlong n, const dfloat *a, MPI_Comm comm){
  dfloat result = 0., gresult = 0.;
  // #pragma omp parallel for reduction(+:result)
  for(dlong i=0; i<n; i++)
    result += a[i]*a[i];

  MPI_Allreduce(&result, &gresult, 1, MPI_DFLOAT, MPI_SUM, comm);
  return sqrt(gresult);
}

dfloat vectorInnerProd(const dlong n, const dfloat *a, const dfloat *b,
                       MPI_Comm comm){
  dfloat result = 0., gresult = 0.;
  // #pragma omp parallel for reduction(+:result)
  for(dlong i=0; i<n; i++)
    result += a[i]*b[i];

  MPI_Allreduce(&result, &gresult, 1, MPI_DFLOAT, MPI_SUM, comm);
  return gresult;
}

dfloat vectorMaxAbs(const dlong n, const dfloat *a, MPI_Comm comm){
  dfloat maxVal=0.0;
  dfloat gmaxVal=0.0;

  //  #pragma omp parallel for reduction(max:maxVal)
  for(dlong i=0; i<n; i++){
    dfloat a2 = (a[i] < 0) ? -a[i] : a[i];
    if(maxVal < a2){
      maxVal = a2;
    }
  }

  MPI_Allreduce(&maxVal, &gmaxVal, 1, MPI_DFLOAT, MPI_MAX, comm);
  return gmaxVal;
}

// returns aDotbc[0] = a\dot b, aDotbc[1] = a\dot c, aDotbc[2] = b\dot b,
void kcycleCombinedOp1(const dlong n, dfloat *aDotbc, const dfloat *a,
                      const dfloat *b, const dfloat *c, const dfloat* w,
                      const bool weighted, MPI_Comm comm) {
  dfloat result[3] = {0.,0.,0.};
  if (weighted) {
    // #pragma omp parallel for reduction(+:aDotb) reduction(+:aDotc) reduction(+:bDotb)
    for(dlong i=0; i<n; i++) {
      result[0] += w[i]*a[i]*b[i];
      result[1] += w[i]*a[i]*c[i];
      result[2] += w[i]*b[i]*b[i];
    }
  } else {
    // #pragma omp parallel for reduction(+:aDotb) reduction(+:aDotc) reduction(+:bDotb)
    for(dlong i=0; i<n; i++) {
      result[0] += a[i]*b[i];
      result[1] += a[i]*c[i];
      result[2] += b[i]*b[i];
    }
  }
  MPI_Allreduce(result,aDotbc,3,MPI_DFLOAT,MPI_SUM,comm);
}

// returns aDotbcd[0] = a\dot b, aDotbcd[1] = a\dot c, aDotbcd[2] = a\dot d,
void kcycleCombinedOp2(const dlong n, dfloat *aDotbcd, const dfloat *a,
                       const dfloat *b, const dfloat *c, const dfloat* d,
                       const dfloat *w, const bool weighted, MPI_Comm comm) {
  dfloat result[3] = {0.,0.,0.};
  if (weighted) {
    // #pragma omp parallel for reduction(+:aDotb) reduction(+:aDotc) reduction(+:aDotd)
    for(dlong i=0; i<n; i++) {
      result[0] += w[i]*a[i]*b[i];
      result[1] += w[i]*a[i]*c[i];
      result[2] += w[i]*a[i]*d[i];
    }
  } else {
    // #pragma omp parallel for reduction(+:aDotb) reduction(+:aDotc) reduction(+:aDotd)
    for(dlong i=0; i<n; i++) {
      result[0] += a[i]*b[i];
      result[1] += a[i]*c[i];
      result[2] += a[i]*d[i];
    }
  }
  MPI_Allreduce(result,aDotbcd,3,MPI_DFLOAT,MPI_SUM,comm);
}



// y = beta*y + alpha*x, and return y\dot y
dfloat vectorAddInnerProd(const dlong n, const dfloat alpha, const dfloat *x,
                          const dfloat beta, dfloat *y,
                          const dfloat *w, const bool weighted, MPI_Comm comm){
  dfloat result = 0.;
  dfloat gresult = 0.;
  if (weighted) {
    if (beta) {
      // #pragma omp parallel for reduction(+:result)
      for(dlong i=0; i<n; i++) {
        y[i] = beta*y[i] + alpha*x[i];
        result += w[i]*y[i]*y[i];
      }
    } else {
      // #pragma omp parallel for reduction(+:result)
      for(dlong i=0; i<n; i++) {
        y[i] = alpha*x[i];
        result += w[i]*y[i]*y[i];
      }
    }
  } else {
    if (beta) {
      // #pragma omp parallel for reduction(+:result)
      for(dlong i=0; i<n; i++) {
        y[i] = beta*y[i] + alpha*x[i];
        result += y[i]*y[i];
      }
    } else {
      // #pragma omp parallel for reduction(+:result)
      for(dlong i=0; i<n; i++) {
        y[i] = alpha*x[i];
        result += y[i]*y[i];
      }
    }
  }
  MPI_Allreduce(&result,&gresult,1,MPI_DFLOAT,MPI_SUM,comm);
  return gresult;
}


//------------------------------------------------------------------------
//
//  Device vector operations
//
//------------------------------------------------------------------------

void vectorSet(const dlong N, const dfloat alpha, occa::memory o_a){
  if (N) vectorSetKernel(N, alpha, o_a);
}

//void vectorRandomize(const dlong m, occa::memory o_a)

void vectorScale(const dlong N, const dfloat alpha, occa::memory o_a){
  if (N) vectorScaleKernel(N, alpha, o_a);
}

void vectorAddScalar(const dlong N, const dfloat alpha, occa::memory o_a){
  if (N) vectorAddScalarKernel(N, alpha, o_a);
}

void vectorAdd(const dlong N, const dfloat alpha, occa::memory o_x,
               const dfloat beta, occa::memory o_y){
  if (N) vectorAddKernel1(N, alpha, beta, o_x, o_y);
}

void vectorAdd(const dlong N, const dfloat alpha, occa::memory o_x,
               const dfloat beta, occa::memory o_y, occa::memory o_z){
  if (N) vectorAddKernel2(N, alpha, beta, o_x, o_y, o_z);
}

void vectorDotStar(const dlong N, occa::memory o_a, occa::memory o_b){
  if (N) vectorDotStarKernel1(N, o_a, o_b);
}

void vectorDotStar(const dlong N, const dfloat alpha, occa::memory o_a,
                   occa::memory o_b, const dfloat beta, occa::memory o_c){
  if (N) vectorDotStarKernel2(N, alpha, beta, o_a, o_b, o_c);
}

//dfloat vectorNorm(const dlong n, occa::memory o_a, MPI_Comm comm)

dfloat vectorInnerProd(const dlong N, occa::memory o_x, occa::memory o_y,
                       MPI_Comm comm){

  dlong numBlocks = (N < NBLOCKS) ? N : NBLOCKS;

  vectorInnerProdKernel(numBlocks,N,o_x,o_y,o_reductionScratch);
  o_reductionScratch.copyTo(reductionScratch,numBlocks*sizeof(dfloat),0);

  dfloat result =0., gresult = 0.;
  //#pragma omp parallel for reduction(+:result)
  for (dlong i=0; i<numBlocks; i++) {
    result += ((dfloat*)reductionScratch)[i];
  }
  MPI_Allreduce(&result, &gresult, 1, MPI_DFLOAT, MPI_SUM, comm);
  return gresult;
}

//dfloat vectorMaxAbs(const dlong n, occa::memory o_a)

// returns aDotbc[0] = a\dot b, aDotbc[1] = a\dot c, aDotbc[2] = b\dot b,
void kcycleCombinedOp1(const dlong N, dfloat *aDotbc, occa::memory o_a,
                       occa::memory o_b, occa::memory o_c, occa::memory o_w,
                       const bool weighted, MPI_Comm comm) {

  dfloat result[3] = {0.,0.,0.};
  dlong numBlocks = (N < NBLOCKS) ? N : NBLOCKS;

  if (weighted) {
    kcycleWeightedCombinedOp1Kernel(numBlocks,N,o_a,o_b,o_c,o_w,o_reductionScratch);
  } else {
    kcycleCombinedOp1Kernel(numBlocks,N,o_a,o_b,o_c,o_reductionScratch);
  }
  o_reductionScratch.copyTo(reductionScratch,3*numBlocks*sizeof(dfloat),0);

  // #pragma omp parallel for reduction(+:aDotb) reduction(+:aDotc) reduction(+:bDotb)
  for(dlong i=0; i<numBlocks; i++) {
    result[0] += ((dfloat*)reductionScratch)[3*i+0];
    result[1] += ((dfloat*)reductionScratch)[3*i+1];
    result[2] += ((dfloat*)reductionScratch)[3*i+2];
  }
  MPI_Allreduce(result,aDotbc,3,MPI_DFLOAT,MPI_SUM,comm);
}

// returns aDotbcd[0] = a\dot b, aDotbcd[1] = a\dot c, aDotbcd[2] = a\dot d,
void kcycleCombinedOp2(const dlong N, dfloat *aDotbcd,
                        occa::memory o_a, occa::memory o_b,
                        occa::memory o_c, occa::memory o_d,
                        occa::memory o_w, const bool weighted, MPI_Comm comm) {

  dfloat result[3] = {0.,0.,0.};
  dlong numBlocks = (N < NBLOCKS) ? N : NBLOCKS;

  if (weighted) {
    kcycleWeightedCombinedOp2Kernel(numBlocks,N,o_a,o_b,o_c,o_d,o_w,o_reductionScratch);
  } else {
    kcycleCombinedOp2Kernel(numBlocks,N,o_a,o_b,o_c,o_d,o_reductionScratch);
  }
  o_reductionScratch.copyTo(reductionScratch,3*numBlocks*sizeof(dfloat),0);

  dfloat aDotb = 0., aDotc = 0., aDotd = 0.;
  // #pragma omp parallel for reduction(+:aDotb) reduction(+:aDotc) reduction(+:aDotd)
  for(dlong i=0; i<numBlocks; i++) {
    result[0] += ((dfloat*)reductionScratch)[3*i+0];
    result[1] += ((dfloat*)reductionScratch)[3*i+1];
    result[2] += ((dfloat*)reductionScratch)[3*i+2];
  }
  MPI_Allreduce(result,aDotbcd,3,MPI_DFLOAT,MPI_SUM,comm);
}

// y = beta*y + alpha*x, and return y\dot y
dfloat vectorAddInnerProd(const dlong N, const dfloat alpha, occa::memory o_x,
                          const dfloat beta, occa::memory o_y,
                          occa::memory o_w, const bool weighted, MPI_Comm comm){

  dfloat result = 0.;
  dfloat gresult = 0.;
  dlong numBlocks = (N < NBLOCKS) ? N : NBLOCKS;

  if (weighted) {
    vectorAddWeightedInnerProdKernel(numBlocks,N,alpha,beta,o_x,o_y,o_w,o_reductionScratch);
  } else {
    vectorAddInnerProdKernel(numBlocks,N,alpha,beta,o_x,o_y,o_reductionScratch);
  }
  o_reductionScratch.copyTo(reductionScratch,numBlocks*sizeof(dfloat),0);

  // #pragma omp parallel for reduction(+:result)
  for (dlong i=0; i<numBlocks; i++) {
    result += ((dfloat*)reductionScratch)[i];
  }
  MPI_Allreduce(&result,&gresult,1,MPI_DFLOAT,MPI_SUM,comm);
  return gresult;
}



} //namespace parAlmond
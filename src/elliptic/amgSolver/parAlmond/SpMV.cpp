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
//  CSR matrix
//
//------------------------------------------------------------------------
void CSR::SpMV(const dfloat alpha, dfloat *x,
               const dfloat beta, dfloat *y) {
  // y[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  if (beta) {
    // #pragma omp parallel for
    for(dlong i=0; i<Nrows; i++){ //local
      dfloat result = 0.0;
      for(dlong jj=rowStarts[i]; jj<rowStarts[i+1]; jj++)
        result += vals[jj]*x[cols[jj]];

      y[i] = alpha*result + beta*y[i];
    }
  } else {
    // #pragma omp parallel for
    for(dlong i=0; i<Nrows; i++){ //local
      dfloat result = 0.0;
      for(dlong jj=rowStarts[i]; jj<rowStarts[i+1]; jj++)
        result += vals[jj]*x[cols[jj]];

      y[i] = alpha*result;
    }
  }
}

void CSR::SpMV(const dfloat alpha, dfloat *x,
               const dfloat beta, const dfloat *y, dfloat *z) {
  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  // #pragma omp parallel for
  for(dlong i=0; i<Nrows; i++){ //local
    dfloat result = 0.0;
    for(dlong jj=rowStarts[i]; jj<rowStarts[i+1]; jj++)
      result += vals[jj]*x[cols[jj]];

    z[i] = alpha*result + beta*y[i];
  }
}

void CSR::SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta,
              occa::memory o_y){
  // y[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  // occaTimerTic(device,"SpMV CSR");
  if (Nrows)
    SpMVcsrKernel1(Nrows, alpha, beta, o_rowStarts, o_cols, o_vals,
                          o_x, o_y);
  // occaTimerToc(device,"SpMV CSR");
}

void CSR::SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta,
              occa::memory o_y, occa::memory o_z){
  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  // occaTimerTic(device,"SpMV CSR");
  if (Nrows)
    SpMVcsrKernel2(Nrows, alpha, beta, o_rowStarts, o_cols, o_vals,
                          o_x, o_y, o_z);
  // occaTimerToc(device,"SpMV CSR");
}


//------------------------------------------------------------------------
//
//  ELL matrix
//
//------------------------------------------------------------------------
void ELL::SpMV(const dfloat alpha, dfloat *x,
               const dfloat beta, dfloat *y) {
  // y[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  if (beta) {
    // #pragma omp parallel for
    for(dlong i=0; i<Nrows; i++){ //local
      dfloat result = 0.0;
      for(dlong c=0; c<nnzPerRow; c++) {
        dlong col = cols[c+nnzPerRow*i];
        if (col>-1) {
          result += vals[c+nnzPerRow*i]*x[col];
        }
      }
      y[i] = alpha*result + beta*y[i];
    }
  } else {
    // #pragma omp parallel for
    for(dlong i=0; i<Nrows; i++){ //local
      dfloat result = 0.0;
      for(dlong c=0; c<nnzPerRow; c++) {
        dlong col = cols[c+nnzPerRow*i];
        if (col>-1) {
          result += vals[c+nnzPerRow*i]*x[col];
        }
      }
      y[i] = alpha*result;
    }
  }
}

void ELL::SpMV(const dfloat alpha, dfloat *x,
               const dfloat beta, const dfloat *y, dfloat *z) {
  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  // #pragma omp parallel for
  for(dlong i=0; i<Nrows; i++){ //local
    dfloat result = 0.0;
    for(dlong c=0; c<nnzPerRow; c++) {
      dlong col = cols[c+nnzPerRow*i];
      if (col>-1) {
        result += vals[c+nnzPerRow*i]*x[col];
      }
    }
    z[i] = alpha*result + beta*y[i];
  }
}

void ELL::SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta,
                occa::memory o_y) {
  // y[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  if(nnzPerRow){
    // occaTimerTic(device,"SpMV ELL");
    SpMVellKernel1(Nrows, nnzPerRow,
                             alpha, beta, o_cols, o_vals, o_x, o_y);
    // occaTimerToc(device,"SpMV ELL");
  }
}

void ELL::SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta,
                occa::memory o_y, occa::memory o_z) {
  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  if(nnzPerRow){
    // occaTimerTic(device,"SpMV ELL");
    SpMVellKernel2(Nrows, nnzPerRow,
                             alpha, beta, o_cols, o_vals, o_x, o_y, o_z);
    // occaTimerToc(device,"SpMV ELL");
  }
}

//------------------------------------------------------------------------
//
//  MCSR matrix
//
//------------------------------------------------------------------------
void MCSR::SpMV(const dfloat alpha, dfloat *x,
                const dfloat beta, dfloat *y){
  // y[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  if (beta) {
    // #pragma omp parallel for
    for(dlong i=0; i<actualRows; i++){ //local
      dlong row = rows[i];
      dfloat result = 0.0;
      for(dlong jj=rowStarts[i]; jj<rowStarts[i+1]; jj++)
        result += vals[jj]*x[cols[jj]];

      y[row] = alpha*result + beta*y[row];
    }
  } else {
    // #pragma omp parallel for
    for(dlong i=0; i<actualRows; i++){ //local
      dlong row = rows[i];
      dfloat result = 0.0;
      for(dlong jj=rowStarts[i]; jj<rowStarts[i+1]; jj++)
        result += vals[jj]*x[cols[jj]];

      y[row] = alpha*result;
    }
  }
}

void MCSR::SpMV(const dfloat alpha, dfloat *x,
                const dfloat beta, const dfloat *y, dfloat *z){
  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  // #pragma omp parallel for
  for(dlong i=0; i<actualRows; i++){ //local
    dlong row = rows[i];
    dfloat result = 0.0;
    for(dlong jj=rowStarts[i]; jj<rowStarts[i+1]; jj++)
      result += vals[jj]*x[cols[jj]];

    z[row] = alpha*result + beta*y[row];
  }
}

void MCSR::SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta,
                occa::memory o_y) {
  // y[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  // occaTimerTic(device,"SpMV MCSR");
  if (actualRows)
    SpMVmcsrKernel1(actualRows, alpha, beta,
                    o_rowStarts, o_rows, o_cols, o_vals,
                    o_x, o_y);
  // occaTimerToc(device,"SpMV MCSR");
}

void MCSR::SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta,
                occa::memory o_y, occa::memory o_z) {
  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  // occaTimerTic(device,"SpMV MCSR");
  if (actualRows)
    SpMVmcsrKernel2(actualRows, alpha, beta,
                    o_rowStarts, o_rows, o_cols, o_vals,
                    o_x, o_y, o_z);
  // occaTimerToc(device,"SpMV MCSR");
}


//------------------------------------------------------------------------
//
//  parCSR matrix
//
//------------------------------------------------------------------------
void parCSR::SpMV(const dfloat alpha, dfloat *x,
                  const dfloat beta, dfloat *y) {

  this->haloExchangeStart(x);

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  diag->SpMV(alpha, x, beta, y);

  this->haloExchangeFinish(x);

  offd->SpMV(alpha, x, 1.0, y);

  //rank 1 correction if there is a nullspace
  if (nullSpace) {
    dfloat gamma = vectorInnerProd(Nrows, null, x, comm)*nullSpacePenalty;
    vectorAdd(Nrows, alpha*gamma, null, 1.0, y);
  }
}

void parCSR::SpMV(const dfloat alpha, dfloat *x,
                  const dfloat beta, const dfloat *y, dfloat *z) {

  this->haloExchangeStart(x);

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  diag->SpMV(alpha, x, beta, y, z);

  this->haloExchangeFinish(x);

  offd->SpMV(alpha, x, 1.0, z);

  //rank 1 correction if there is a nullspace
  if (nullSpace) {
    dfloat gamma = vectorInnerProd(Nrows, null, x, comm)*nullSpacePenalty;
    vectorAdd(Nrows, alpha*gamma, null, 1.0, z);
  }
}

void parCSR::SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta,
                  occa::memory o_y) {

  this->haloExchangeStart(o_x);

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  diag->SpMV(alpha, o_x, beta, o_y);

  this->haloExchangeFinish(o_x);

  offd->SpMV(alpha, o_x, 1.0, o_y);

  //rank 1 correction if there is a nullspace
  if (nullSpace) {
    dfloat gamma = vectorInnerProd(Nrows, o_null, o_x, comm)*nullSpacePenalty;
    vectorAdd(Nrows, alpha*gamma, o_null, 1.0, o_y);
  }
}

void parCSR::SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta,
                  occa::memory o_y, occa::memory o_z) {

  this->haloExchangeStart(o_x);

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  diag->SpMV(alpha, o_x, beta, o_y, o_z);

  this->haloExchangeFinish(o_x);

  offd->SpMV(alpha, o_x, 1.0, o_z);

  //rank 1 correction if there is a nullspace
  if (nullSpace) {
    dfloat gamma = vectorInnerProd(Nrows, o_null, o_x, comm)*nullSpacePenalty;
    vectorAdd(Nrows, alpha*gamma, o_null, 1.0, o_z);
  }
}

//------------------------------------------------------------------------
//
//  parHYB matrix
//
//------------------------------------------------------------------------
void parHYB::SpMV(const dfloat alpha, dfloat *x,
                  const dfloat beta, dfloat *y) {

  this->haloExchangeStart(x);

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  E->SpMV(alpha, x, beta, y);

  this->haloExchangeFinish(x);

  C->SpMV(alpha, x, 1.0, y);

  //rank 1 correction if there is a nullspace
  if (nullSpace) {
    dfloat gamma = vectorInnerProd(Nrows, null, x, comm)*nullSpacePenalty;
    vectorAdd(Nrows, alpha*gamma, null, 1.0, y);
  }
}

void parHYB::SpMV(const dfloat alpha, dfloat *x,
                  const dfloat beta, const dfloat *y, dfloat *z) {

  this->haloExchangeStart(x);

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  E->SpMV(alpha, x, beta, y, z);

  this->haloExchangeFinish(x);

  C->SpMV(alpha, x, 1.0, z);

  //rank 1 correction if there is a nullspace
  if (nullSpace) {
    dfloat gamma = vectorInnerProd(Nrows, null, x, comm)*nullSpacePenalty;
    vectorAdd(Nrows, alpha*gamma, null, 1.0, z);
  }
}

void parHYB::SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta,
                  occa::memory o_y) {

  this->haloExchangeStart(o_x);

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  E->SpMV(alpha, o_x, beta, o_y);

  this->haloExchangeFinish(o_x);

  C->SpMV(alpha, o_x, 1.0, o_y);

  //rank 1 correction if there is a nullspace
  if (nullSpace) {
    dfloat gamma = vectorInnerProd(Nrows, o_null, o_x, comm)*nullSpacePenalty;
    vectorAdd(Nrows, alpha*gamma, o_null, 1.0, o_y);
  }
}

void parHYB::SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta,
                  occa::memory o_y, occa::memory o_z) {

  this->haloExchangeStart(o_x);

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  E->SpMV(alpha, o_x, beta, o_y, o_z);

  this->haloExchangeFinish(o_x);

  C->SpMV(alpha, o_x, 1.0, o_z);

  //rank 1 correction if there is a nullspace
  if (nullSpace) {
    dfloat gamma = vectorInnerProd(Nrows, o_null, o_x, comm)*nullSpacePenalty;
    vectorAdd(Nrows, alpha*gamma, o_null, 1.0, o_z);
  }
}


} //namespace parAlmond
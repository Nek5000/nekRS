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

parCSR* strongGraph(parCSR *A){

  const dlong N = A->Nrows;
  const dlong M = A->Ncols;

  parCSR *C = new parCSR(N, M);

  C->diag->rowStarts = (dlong *) calloc(N+1,sizeof(dlong));
  C->offd->rowStarts = (dlong *) calloc(N+1,sizeof(dlong));

  dfloat *maxOD;
  if (N) maxOD = (dfloat *) calloc(N,sizeof(dfloat));

  dfloat *diagA = A->diagA;

  // #pragma omp parallel for
  for(dlong i=0; i<N; i++){
    const int sign = (diagA[i] >= 0) ? 1:-1;
    //    const dfloat Aii = fabs(diagA[i]);

    //find maxOD
    //local entries
    dlong Jstart = A->diag->rowStarts[i];
    dlong Jend   = A->diag->rowStarts[i+1];
    for(dlong jj= Jstart; jj<Jend; jj++){
      dlong col = A->diag->cols[jj];
      if (col==i) continue;
      //      dfloat Ajj = fabs(diagA[col]);
      dfloat OD = -sign*A->diag->vals[jj]; // /(sqrt(Aii)*sqrt(Ajj));
      if(OD > maxOD[i]) maxOD[i] = OD;
    }
    //non-local entries
    Jstart = A->offd->rowStarts[i],
    Jend   = A->offd->rowStarts[i+1];
    for(dlong jj= Jstart; jj<Jend; jj++){
      //      dlong col = A->offd->cols[jj];
      //      dfloat Ajj = fabs(diagA[col]);
      dfloat OD = -sign*A->offd->vals[jj]; // /(sqrt(Aii)*sqrt(Ajj));
      if(OD > maxOD[i]) maxOD[i] = OD;
    }

    int diag_strong_per_row = 1; // diagonal entry
    //local entries
    Jstart = A->diag->rowStarts[i],
    Jend   = A->diag->rowStarts[i+1];
    for(dlong jj = Jstart; jj<Jend; jj++){
      dlong col = A->diag->cols[jj];
      if (col==i) continue;
      //      dfloat Ajj = fabs(diagA[col]);
      dfloat OD = -sign*A->diag->vals[jj]; // /(sqrt(Aii)*sqrt(Ajj));
      if(OD > COARSENTHREASHOLD*maxOD[i]) diag_strong_per_row++;
    }
    int offd_strong_per_row = 0;
    //non-local entries
    Jstart = A->offd->rowStarts[i], Jend = A->offd->rowStarts[i+1];
    for(dlong jj= Jstart; jj<Jend; jj++){
      dlong col = A->offd->cols[jj];
      //      dfloat Ajj = fabs(diagA[col]);
      dfloat OD = -sign*A->offd->vals[jj]; // /(sqrt(Aii)*sqrt(Ajj));
      if(OD > COARSENTHREASHOLD*maxOD[i]) offd_strong_per_row++;
    }

    C->diag->rowStarts[i+1] = diag_strong_per_row;
    C->offd->rowStarts[i+1] = offd_strong_per_row;
  }

  // cumulative sum
  for(dlong i=1; i<N+1 ; i++) {
    C->diag->rowStarts[i] += C->diag->rowStarts[i-1];
    C->offd->rowStarts[i] += C->offd->rowStarts[i-1];
  }
  C->diag->nnz = C->diag->rowStarts[N];
  C->offd->nnz = C->offd->rowStarts[N];

  C->diag->cols = (dlong *) malloc(C->diag->nnz*sizeof(dlong));
  C->offd->cols = (dlong *) malloc(C->offd->nnz*sizeof(dlong));
  // C->diag->vals = (dfloat *) malloc(0);
  // C->offd->vals = (dfloat *) malloc(0);

  // fill in the columns for strong connections
  // #pragma omp parallel for
  for(dlong i=0; i<N; i++){
    const int sign = (diagA[i] >= 0) ? 1:-1;
    //    const dfloat Aii = fabs(diagA[i]);

    dlong diagCounter = C->diag->rowStarts[i];
    dlong offdCounter = C->offd->rowStarts[i];

    //local entries
    dlong Jstart = A->diag->rowStarts[i];
    dlong Jend   = A->diag->rowStarts[i+1];
    for(dlong jj = Jstart; jj<Jend; jj++){
      dlong col = A->diag->cols[jj];
      if (col==i) {
        C->diag->cols[diagCounter++] = col;// diag entry
        continue;
      }
      //      dfloat Ajj = fabs(diagA[col]);
      dfloat OD = -sign*A->diag->vals[jj]; // /(sqrt(Aii)*sqrt(Ajj));
      if(OD > COARSENTHREASHOLD*maxOD[i])
        C->diag->cols[diagCounter++] = col;
    }
    Jstart = A->offd->rowStarts[i], Jend = A->offd->rowStarts[i+1];
    for(dlong jj = Jstart; jj<Jend; jj++){
      dlong col = A->offd->cols[jj];
      //      dfloat Ajj = fabs(diagA[col]);
      dfloat OD = -sign*A->offd->vals[jj]; // /(sqrt(Aii)*sqrt(Ajj));
      if(OD > COARSENTHREASHOLD*maxOD[i])
        C->offd->cols[offdCounter++] = col;
    }
  }
  if(N) free(maxOD);

  return C;
}

} //namespace parAlmond

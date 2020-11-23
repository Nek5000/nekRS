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

parCSR *constructProlongation(parCSR *A, hlong *FineToCoarse,
                            hlong *globalAggStarts, dfloat **nullCoarseA){
  // MPI info
  int rank, size;
  MPI_Comm_rank(A->comm, &rank);
  MPI_Comm_size(A->comm, &size);

  const dlong N = A->Nrows;

  const hlong globalAggOffset = globalAggStarts[rank];
  const dlong NCoarse = (dlong) (globalAggStarts[rank+1]-globalAggStarts[rank]); //local num agg

  parCSR* P = new parCSR(N, NCoarse, A->comm, A->device);

  P->globalRowStarts = A->globalRowStarts;
  P->globalColStarts = globalAggStarts;

  P->diag->rowStarts = (dlong *) calloc(N+1, sizeof(dlong));
  P->offd->rowStarts = (dlong *) calloc(N+1, sizeof(dlong));

  // each row has exactly one nonzero
  for(dlong i=0; i<N; i++) {
    hlong col = FineToCoarse[i];
    if ((col>globalAggOffset-1)&&(col<globalAggOffset+NCoarse)) {
      P->diag->rowStarts[i+1]++;
    } else {
      P->offd->rowStarts[i+1]++;
    }
  }
  for(dlong i=0; i<N; i++) {
    P->diag->rowStarts[i+1] += P->diag->rowStarts[i];
    P->offd->rowStarts[i+1] += P->offd->rowStarts[i];
  }
  P->diag->nnz = P->diag->rowStarts[N];
  P->offd->nnz = P->offd->rowStarts[N];

  // Halo setup
  hlong *colIds = (hlong *) malloc(P->offd->nnz*sizeof(hlong));
  dlong cnt=0;
  for (dlong i=0;i<N;i++) {
    hlong col = FineToCoarse[i];
    if ((col<globalAggOffset)||(col>globalAggOffset+NCoarse-1))
      colIds[cnt++] = col;
  }
  P->haloSetup(colIds);

  P->diag->cols = (dlong *)  calloc(P->diag->nnz, sizeof(dlong));
  P->diag->vals = (dfloat *) calloc(P->diag->nnz, sizeof(dfloat));
  P->offd->cols = (dlong *)  calloc(P->offd->nnz, sizeof(dlong));
  P->offd->vals = (dfloat *) calloc(P->offd->nnz, sizeof(dfloat));

  dlong diagCnt = 0;
  dlong offdCnt = 0;
  for(dlong i=0; i<N; i++) {
    hlong col = FineToCoarse[i];
    if ((col>globalAggStarts[rank]-1)&&(col<globalAggStarts[rank+1])) {
      P->diag->cols[diagCnt  ] = (dlong) (col - globalAggOffset); //local index
      P->diag->vals[diagCnt++] = A->null[i];
    } else {
      P->offd->cols[offdCnt  ] = colIds[offdCnt];
      P->offd->vals[offdCnt++] = A->null[i];
    }
  }

  // normalize the columns of P
  *nullCoarseA = (dfloat *) calloc(P->Ncols,sizeof(dfloat));

  //add local nonzeros
  for(dlong i=0; i<P->diag->nnz; i++)
    (*nullCoarseA)[P->diag->cols[i]] += P->diag->vals[i] * P->diag->vals[i];

  //add nonlocal nonzeros
  for(dlong i=0; i<P->offd->nnz; i++)
    (*nullCoarseA)[P->offd->cols[i]] += P->offd->vals[i] * P->offd->vals[i];

  ogsGatherScatter((*nullCoarseA), ogsDfloat,  ogsAdd, P->ogs);

  for(dlong i=0; i<NCoarse; i++)
    (*nullCoarseA)[i] = sqrt((*nullCoarseA)[i]);

  for(dlong i=NCoarse; i<P->Ncols; i++)
    (*nullCoarseA)[i] = 0.;

  ogsGatherScatter((*nullCoarseA), ogsDfloat,  ogsAdd, P->ogs);

  for(dlong i=0; i<P->diag->nnz; i++)
    P->diag->vals[i] /= (*nullCoarseA)[P->diag->cols[i]];
  for(dlong i=0; i<P->offd->nnz; i++)
    P->offd->vals[i] /= (*nullCoarseA)[P->offd->cols[i]];

  return P;
}

} //namespace parAlmond
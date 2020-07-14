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

parCSR *galerkinProd(parCSR *A, parCSR *P){

  // MPI info
  int rank, size;
  MPI_Comm_rank(A->comm, &rank);
  MPI_Comm_size(A->comm, &size);

  hlong *globalAggStarts = P->globalColStarts;
  hlong globalAggOffset = globalAggStarts[rank];

  //The galerkin product can be computed as
  // (P^T A P)_IJ = sum_{i in Agg_I} sum_{j in Agg_J} P_iI A_ij P_jJ
  // Since each row of P has only one entry, we can share the necessary
  // P entries, form the products, and send them to their destination rank

  const dlong N = A->Nrows;
  const dlong M = A->Ncols;

  //printf("Level has %d rows, and is making %d aggregates\n", N, globalAggStarts[rank+1]-globalAggStarts[rank]);

  hlong  *Pcols = (hlong  *) calloc(M,sizeof(hlong));
  dfloat *Pvals = (dfloat *) calloc(M,sizeof(dfloat));

  //record the entries of P that this rank has
  dlong cnt =0;
  for (dlong i=0;i<N;i++) {
    for (dlong j=P->diag->rowStarts[i];j<P->diag->rowStarts[i+1];j++) {
      Pcols[cnt] = P->diag->cols[j] + globalAggOffset; //global ID
      Pvals[cnt] = P->diag->vals[j];
      cnt++;
    }
    for (dlong j=P->offd->rowStarts[i];j<P->offd->rowStarts[i+1];j++) {
      Pcols[cnt] = P->colMap[P->offd->cols[j]]; //global ID
      Pvals[cnt] = P->offd->vals[j];
      cnt++;
    }
  }

  //fill the halo region
  ogsGatherScatter(Pcols, ogsHlong,  ogsAdd, A->ogs);
  ogsGatherScatter(Pvals, ogsDfloat, ogsAdd, A->ogs);



  dlong sendNtotal = A->diag->nnz+A->offd->nnz;
  nonzero_t *sendPTAP = (nonzero_t *) calloc(sendNtotal,sizeof(nonzero_t));

  // Make the MPI_NONZERO_T data type
  nonzero_t NZ;
  MPI_Datatype MPI_NONZERO_T;
  MPI_Datatype dtype[3] = {MPI_HLONG, MPI_HLONG, MPI_DFLOAT};
  int blength[3] = {1, 1, 1};
  MPI_Aint addr[3], displ[3];
  MPI_Get_address ( &(NZ.row), addr+0);
  MPI_Get_address ( &(NZ.col), addr+1);
  MPI_Get_address ( &(NZ.val), addr+2);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  MPI_Type_create_struct (3, blength, displ, dtype, &MPI_NONZERO_T);
  MPI_Type_commit (&MPI_NONZERO_T);

  //form the fine PTAP products
  cnt =0;
  for (dlong i=0;i<N;i++) {
    dlong start = A->diag->rowStarts[i];
    dlong end   = A->diag->rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong  col = A->diag->cols[j];
      const dfloat val = A->diag->vals[j];

      sendPTAP[cnt].row = Pcols[i];
      sendPTAP[cnt].col = Pcols[col];
      sendPTAP[cnt].val = val*Pvals[i]*Pvals[col];
      cnt++;
    }
    start = A->offd->rowStarts[i];
    end   = A->offd->rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong  col = A->offd->cols[j];
      const dfloat val = A->offd->vals[j];

      sendPTAP[cnt].row = Pcols[i];
      sendPTAP[cnt].col = Pcols[col];
      sendPTAP[cnt].val = val*Pvals[i]*Pvals[col];
      cnt++;
    }
  }

  free(Pcols);
  free(Pvals);

  //sort entries by the coarse row and col
  qsort(sendPTAP, sendNtotal, sizeof(nonzero_t), compareNonZeroByRow);

  //count number of non-zeros we're sending
  int *sendCounts = (int *) calloc(size,sizeof(int));
  int *recvCounts = (int *) calloc(size,sizeof(int));
  int *sendOffsets = (int *) calloc(size+1,sizeof(int));
  int *recvOffsets = (int *) calloc(size+1,sizeof(int));

  int r=0;
  for(dlong i=0;i<sendNtotal;++i) {
    hlong id = sendPTAP[i].row;
    while(id>=globalAggStarts[r+1]) r++;
    sendCounts[r]++;
  }

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(sendCounts, 1, MPI_INT,
               recvCounts, 1, MPI_INT, A->comm);

  // find send and recv offsets for gather
  for(int r=0;r<size;++r){
    sendOffsets[r+1] = sendOffsets[r] + sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r] + recvCounts[r];
  }
  dlong recvNtotal = recvOffsets[size];

  nonzero_t *recvPTAP = (nonzero_t *) calloc(recvNtotal,sizeof(nonzero_t));

  MPI_Alltoallv(sendPTAP, sendCounts, sendOffsets, MPI_NONZERO_T,
                recvPTAP, recvCounts, recvOffsets, MPI_NONZERO_T,
                A->comm);

  //clean up
  MPI_Barrier(A->comm);
  free(sendPTAP);
  free(sendCounts); free(recvCounts);
  free(sendOffsets); free(recvOffsets);

  //sort entries by the coarse row and col
  qsort(recvPTAP, recvNtotal, sizeof(nonzero_t), compareNonZeroByRow);

  //count total number of nonzeros;
  dlong nnz =0;
  if (recvNtotal) nnz++;
  for (dlong i=1;i<recvNtotal;i++)
    if ((recvPTAP[i].row!=recvPTAP[i-1].row)||
        (recvPTAP[i].col!=recvPTAP[i-1].col)) nnz++;

  nonzero_t *PTAP = (nonzero_t *) calloc(nnz,sizeof(nonzero_t));

  //compress nonzeros
  nnz = 0;
  if (recvNtotal) PTAP[nnz++] = recvPTAP[0];
  for (dlong i=1;i<recvNtotal;i++) {
    if ((recvPTAP[i].row!=recvPTAP[i-1].row)||
        (recvPTAP[i].col!=recvPTAP[i-1].col)) {
      PTAP[nnz++] = recvPTAP[i];
    } else {
      PTAP[nnz-1].val += recvPTAP[i].val;
    }
  }

  //clean up
  MPI_Barrier(A->comm);
  free(recvPTAP);

  dlong numAggs = (dlong) (globalAggStarts[rank+1]-globalAggStarts[rank]); //local number of aggregates

  parCSR *Ac = new parCSR(numAggs, numAggs, A->comm, A->device);

  Ac->globalRowStarts = globalAggStarts;
  Ac->globalColStarts = globalAggStarts;

  Ac->diag->rowStarts = (dlong *) calloc(numAggs+1, sizeof(dlong));
  Ac->offd->rowStarts = (dlong *) calloc(numAggs+1, sizeof(dlong));

  for (dlong n=0;n<nnz;n++) {
    dlong row = (dlong) (PTAP[n].row - globalAggOffset);
    if ((PTAP[n].col > globalAggStarts[rank]-1)&&
        (PTAP[n].col < globalAggStarts[rank+1])) {
      Ac->diag->rowStarts[row+1]++;
    } else {
      Ac->offd->rowStarts[row+1]++;
    }
  }

  // cumulative sum
  for(dlong i=0; i<numAggs; i++) {
    Ac->diag->rowStarts[i+1] += Ac->diag->rowStarts[i];
    Ac->offd->rowStarts[i+1] += Ac->offd->rowStarts[i];
  }
  Ac->diag->nnz = Ac->diag->rowStarts[numAggs];
  Ac->offd->nnz = Ac->offd->rowStarts[numAggs];

  // Halo setup
  hlong *colIds = (hlong *) malloc(Ac->offd->nnz*sizeof(hlong));
  cnt=0;
  for (dlong n=0;n<nnz;n++) {
    if ((PTAP[n].col <= (globalAggStarts[rank]-1))||
        (PTAP[n].col >= globalAggStarts[rank+1])) {
      colIds[cnt++] = PTAP[n].col;
    }
  }
  Ac->haloSetup(colIds);

  //fill the CSR matrices
  Ac->diagA   = (dfloat *) calloc(Ac->Ncols, sizeof(dfloat));
  Ac->diagInv = (dfloat *) calloc(Ac->Ncols, sizeof(dfloat));
  Ac->diag->cols = (dlong *)  calloc(Ac->diag->nnz, sizeof(dlong));
  Ac->offd->cols = (dlong *)  calloc(Ac->offd->nnz, sizeof(dlong));
  Ac->diag->vals = (dfloat *) calloc(Ac->diag->nnz, sizeof(dfloat));
  Ac->offd->vals = (dfloat *) calloc(Ac->offd->nnz, sizeof(dfloat));
  dlong diagCnt = 0;
  dlong offdCnt = 0;
  for (dlong n=0;n<nnz;n++) {
    if ((PTAP[n].col > globalAggStarts[rank]-1)&&
        (PTAP[n].col < globalAggStarts[rank+1])) {
      Ac->diag->cols[diagCnt] = (dlong) (PTAP[n].col - globalAggOffset);
      Ac->diag->vals[diagCnt] = PTAP[n].val;

      //record the diagonal
      dlong row = (dlong) (PTAP[n].row - globalAggOffset);
      if (row==Ac->diag->cols[diagCnt])
        Ac->diagA[row] = Ac->diag->vals[diagCnt];

      diagCnt++;
    } else {
      Ac->offd->cols[offdCnt] = colIds[offdCnt];
      Ac->offd->vals[offdCnt] = PTAP[n].val;
      offdCnt++;
    }
  }

  //compute the inverse diagonal
  for (dlong n=0;n<Ac->Nrows;n++) Ac->diagInv[n] = 1.0/Ac->diagA[n];

  //propagate nullspace flag
  Ac->nullSpace = A->nullSpace;
  Ac->nullSpacePenalty = A->nullSpacePenalty;

  //clean up
  MPI_Barrier(A->comm);
  MPI_Type_free(&MPI_NONZERO_T);
  free(colIds);
  free(PTAP);

  return Ac;
}

} //namespace parAlmond
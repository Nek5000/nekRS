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

static void formAggregatesDefault(parCSR *A, parCSR *C,
                                  hlong* FineToCoarse,
                                  hlong* globalAggStarts);

static void formAggregatesLPSCN(parCSR *A, parCSR *C,
                                hlong* FineToCoarse,
                                hlong* globalAggStarts,
                                setupAide options);

void formAggregates(parCSR *A, parCSR *C,
                    hlong* FineToCoarse,
                    hlong* globalAggStarts,
                    setupAide options) {

  if (options.compareArgs("PARALMOND AGGREGATION STRATEGY", "DEFAULT")) {
    formAggregatesDefault(A, C, FineToCoarse, globalAggStarts);
  } else if (options.compareArgs("PARALMOND AGGREGATION STRATEGY", "LPSCN")) {
    formAggregatesLPSCN(A, C, FineToCoarse, globalAggStarts, options);
  } else {
    printf("WARNING:  Missing or bad value for option PARALMOND AGGREGATION STRATEGY.  Using default.\n");
    formAggregatesDefault(A, C, FineToCoarse, globalAggStarts);
  }
}

/*****************************************************************************/
// Default aggregation algorithm
//
// Parallel Distance 2 Maximal Independant Set (MIS-2) graph partitioning
//
/*****************************************************************************/

static void formAggregatesDefault(parCSR *A, parCSR *C,
                                  hlong* FineToCoarse,
                                  hlong* globalAggStarts){

  int rank, size;
  MPI_Comm_rank(A->comm, &rank);
  MPI_Comm_size(A->comm, &size);

  const dlong N   = C->Nrows;
  const dlong M   = C->Ncols;
  const dlong diagNNZ = C->diag->nnz;
  const dlong offdNNZ = C->offd->nnz;

  dfloat *rands = (dfloat *) calloc(M, sizeof(dfloat));
  int   *states = (int *)    calloc(M, sizeof(int));

  dfloat *Tr = (dfloat *) calloc(M, sizeof(dfloat));
  int    *Ts = (int *)    calloc(M, sizeof(int));
  hlong  *Ti = (hlong *)  calloc(M, sizeof(hlong));
  hlong  *Tc = (hlong *)  calloc(M, sizeof(hlong));

  hlong *globalRowStarts = A->globalRowStarts;

  for(dlong i=0; i<N; i++)
    rands[i] = (dfloat) drand48();

  // add the number of non-zeros in each column
  int *colCnt = (int *) calloc(M,sizeof(int));
  for(dlong i=0; i<diagNNZ; i++)
    colCnt[C->diag->cols[i]]++;

  for(dlong i=0; i<offdNNZ; i++)
    colCnt[C->offd->cols[i]]++;

  //gs for total column counts
  ogsGatherScatter(colCnt, ogsInt, ogsAdd, A->ogs);

  //add random pertubation
  for(int i=0;i<N;++i)
    rands[i] += colCnt[i];

  //gs to fill halo region
  ogsGatherScatter(rands, ogsDfloat, ogsAdd, A->ogs);

  hlong done = 0;
  while(!done){
    // first neighbours
    // #pragma omp parallel for
    for(dlong i=0; i<N; i++){

      int smax = states[i];
      dfloat rmax = rands[i];
      hlong imax = i + globalRowStarts[rank];

      if(smax != 1){
        //local entries
        for(dlong jj=C->diag->rowStarts[i];jj<C->diag->rowStarts[i+1];jj++){
          const dlong col = C->diag->cols[jj];
          if (col==i) continue;
          if(customLess(smax, rmax, imax, states[col], rands[col], col + globalRowStarts[rank])){
            smax = states[col];
            rmax = rands[col];
            imax = col + globalRowStarts[rank];
          }
        }
        //nonlocal entries
        for(dlong jj=C->offd->rowStarts[i];jj<C->offd->rowStarts[i+1];jj++){
          const dlong col = C->offd->cols[jj];
          if(customLess(smax, rmax, imax, states[col], rands[col], A->colMap[col])) {
            smax = states[col];
            rmax = rands[col];
            imax = A->colMap[col];
          }
        }
      }
      Ts[i] = smax;
      Tr[i] = rmax;
      Ti[i] = imax;
    }

    //share results
    for (dlong n=N;n<M;n++) {
      Tr[n] = 0.;
      Ts[n] = 0;
      Ti[n] = 0;
    }
    ogsGatherScatter(Tr, ogsDfloat, ogsAdd, A->ogs);
    ogsGatherScatter(Ts, ogsInt,    ogsAdd, A->ogs);
    ogsGatherScatter(Ti, ogsHlong,  ogsAdd, A->ogs);

    // second neighbours
    // #pragma omp parallel for
    for(dlong i=0; i<N; i++){
      int    smax = Ts[i];
      dfloat rmax = Tr[i];
      hlong  imax = Ti[i];

      //local entries
      for(dlong jj=C->diag->rowStarts[i];jj<C->diag->rowStarts[i+1];jj++){
        const dlong col = C->diag->cols[jj];
        if (col==i) continue;
        if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
          smax = Ts[col];
          rmax = Tr[col];
          imax = Ti[col];
        }
      }
      //nonlocal entries
      for(dlong jj=C->offd->rowStarts[i];jj<C->offd->rowStarts[i+1];jj++){
        const dlong col = C->offd->cols[jj];
        if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
          smax = Ts[col];
          rmax = Tr[col];
          imax = Ti[col];
        }
      }

      // if I am the strongest among all the 1 and 2 ring neighbours
      // I am an MIS node
      if((states[i] == 0) && (imax == (i + globalRowStarts[rank])))
        states[i] = 1;

      // if there is an MIS node within distance 2, I am removed
      if((states[i] == 0) && (smax == 1))
        states[i] = -1;
    }

    //share results
    for (dlong n=N;n<M;n++) states[n] = 0;
    ogsGatherScatter(states, ogsInt, ogsAdd, A->ogs);

    // if number of undecided nodes = 0, algorithm terminates
    hlong cnt = std::count(states, states+N, 0);
    MPI_Allreduce(&cnt,&done,1,MPI_HLONG, MPI_SUM,A->comm);
    done = (done == 0) ? 1 : 0;
  }

  dlong numAggs = 0;
  dlong *gNumAggs = (dlong *) calloc(size,sizeof(dlong));

  // count the coarse nodes/aggregates
  for(dlong i=0; i<N; i++)
    if(states[i] == 1) numAggs++;

  MPI_Allgather(&numAggs,1,MPI_DLONG,gNumAggs,1,MPI_DLONG,A->comm);

  globalAggStarts[0] = 0;
  for (int r=0;r<size;r++)
    globalAggStarts[r+1] = globalAggStarts[r] + gNumAggs[r];

  numAggs = 0;
  // enumerate the coarse nodes/aggregates
  for(dlong i=0; i<N; i++) {
    if(states[i] == 1) {
      FineToCoarse[i] = globalAggStarts[rank] + numAggs++;
    } else {
      FineToCoarse[i] = -1;
    }
  }
  for(dlong i=N; i<M; i++) FineToCoarse[i] = 0;

  //share the initial aggregate flags
  ogsGatherScatter(FineToCoarse, ogsHlong, ogsAdd, A->ogs);

  // form the aggregates
  // #pragma omp parallel for
  for(dlong i=0; i<N; i++){
    int   smax = states[i];
    dfloat rmax = rands[i];
    hlong  imax = i + globalRowStarts[rank];
    hlong  cmax = FineToCoarse[i];

    if(smax != 1){
      //local entries
      for(dlong jj=C->diag->rowStarts[i];jj<C->diag->rowStarts[i+1];jj++){
        const dlong col = C->diag->cols[jj];
        if (col==i) continue;
        if(customLess(smax, rmax, imax, states[col], rands[col], col + globalRowStarts[rank])){
          smax = states[col];
          rmax = rands[col];
          imax = col + globalRowStarts[rank];
          cmax = FineToCoarse[col];
        }
      }
      //nonlocal entries
      for(dlong jj=C->offd->rowStarts[i];jj<C->offd->rowStarts[i+1];jj++){
        const dlong col = C->offd->cols[jj];
        if(customLess(smax, rmax, imax, states[col], rands[col], A->colMap[col])){
          smax = states[col];
          rmax = rands[col];
          imax = A->colMap[col];
          cmax = FineToCoarse[col];
        }
      }
    }
    Ts[i] = smax;
    Tr[i] = rmax;
    Ti[i] = imax;
    Tc[i] = cmax;

    if((states[i] == -1) && (smax == 1) && (cmax > -1))
      FineToCoarse[i] = cmax;
  }

  //share results
  for (dlong n=N;n<M;n++) {
    FineToCoarse[n] = 0;
    Tr[n] = 0.;
    Ts[n] = 0;
    Ti[n] = 0;
    Tc[n] = 0;
  }
  ogsGatherScatter(FineToCoarse, ogsHlong,  ogsAdd, A->ogs);
  ogsGatherScatter(Tr,     ogsDfloat, ogsAdd, A->ogs);
  ogsGatherScatter(Ts,     ogsInt,    ogsAdd, A->ogs);
  ogsGatherScatter(Ti,     ogsHlong,  ogsAdd, A->ogs);
  ogsGatherScatter(Tc,     ogsHlong,  ogsAdd, A->ogs);

  // second neighbours
  // #pragma omp parallel for
  for(dlong i=0; i<N; i++){
    int    smax = Ts[i];
    dfloat rmax = Tr[i];
    hlong  imax = Ti[i];
    hlong  cmax = Tc[i];

    //local entries
    for(dlong jj=C->diag->rowStarts[i];jj<C->diag->rowStarts[i+1];jj++){
      const dlong col = C->diag->cols[jj];
      if (col==i) continue;
      if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
        smax = Ts[col];
        rmax = Tr[col];
        imax = Ti[col];
        cmax = Tc[col];
      }
    }
    //nonlocal entries
    for(dlong jj=C->offd->rowStarts[i];jj<C->offd->rowStarts[i+1];jj++){
      const dlong col = C->offd->cols[jj];
      if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
        smax = Ts[col];
        rmax = Tr[col];
        imax = Ti[col];
        cmax = Tc[col];
      }
    }

    if((states[i] == -1) && (smax == 1) && (cmax > -1))
      FineToCoarse[i] = cmax;
  }

  //share results
  for (dlong n=N;n<M;n++) FineToCoarse[n] = 0;
  ogsGatherScatter(FineToCoarse, ogsHlong,  ogsAdd, A->ogs);

  free(rands);
  free(states);
  free(Tr);
  free(Ts);
  free(Ti);
  free(Tc);

  delete C;
}

/*****************************************************************************/
// Alan's "locally partial strong connected nodes" (LPSCN) algorithm.
/*****************************************************************************/

typedef struct{
  dlong  index;
  dlong  Nnbs;
  dlong LNN;
} nbs_t;


int compareNBSmaxLPSCN(const void *a, const void *b)
{
  nbs_t *pa = (nbs_t *)a;
  nbs_t *pb = (nbs_t *)b;

  if (pa->Nnbs + pa->LNN < pb->Nnbs + pb->LNN)  return +1;
  if (pa->Nnbs + pa->LNN > pb->Nnbs + pb->LNN)  return -1;

  if (pa->index < pa->index ) return +1;
  if (pa->index > pa->index ) return -1;

  return 0;
}


int compareNBSminLPSCN(const void *a, const void *b)
{
  nbs_t *pa = (nbs_t *)a;
  nbs_t *pb = (nbs_t *)b;

  if (pa->Nnbs + pa->LNN < pb->Nnbs + pb->LNN)  return -1;
  if (pa->Nnbs + pa->LNN > pb->Nnbs + pb->LNN)  return +1;

  if (pa->index < pa->index ) return +1;
  if (pa->index > pa->index ) return -1;

  return 0;
}

static void formAggregatesLPSCN(parCSR *A, parCSR *C,
                                hlong* FineToCoarse,
                                hlong* globalAggStarts,
                                setupAide options) {
  int rank, size;
  MPI_Comm_rank(A->comm, &rank);
  MPI_Comm_size(A->comm, &size);

  const dlong N   = C->Nrows;
  const dlong M   = C->Ncols;
  const dlong diagNNZ = C->diag->nnz;

  dlong   *states = (dlong *)    calloc(M, sizeof(dlong));    //  M > N
  for(dlong i=0; i<N; i++)   // initialize states to -1
    states[i] = -1;

  hlong *globalRowStarts = A->globalRowStarts;

  // construct the local neigbors
  nbs_t *V = (nbs_t *) calloc(N,sizeof(nbs_t));

  for(dlong i=0; i<N; i++){
    V[i].index = i;
    V[i].Nnbs  = C->diag->rowStarts[i+1] - C->diag->rowStarts[i];
    dlong dist = C->diag->cols[C->diag->rowStarts[i+1]-1] - C->diag->cols[C->diag->rowStarts[i]];

    if (dist  > 0 )
      V[i].LNN   =  V[i].Nnbs/dist;
    else
      V[i].LNN = 0;
  }

  // sort the fake nbs
  if (options.compareArgs("PARALMOND LPSCN ORDERING", "MAX")) {
    qsort(V,N,sizeof(nbs_t),compareNBSmaxLPSCN);
  } else if (options.compareArgs("PARALMOND LPSCN ORDERING", "MIN")) {
    qsort(V,N,sizeof(nbs_t),compareNBSminLPSCN);
  } else { //default NONE
   //do nothing, use the native ordering
  }

  // First aggregates
  dlong Agg_num = 0;
  //#pragma omp parallel for
  for(dlong i=0; i<N; i++){

    if(V[i].Nnbs<=1) continue; //skip

    dlong id = V[i].index;

    if (states[id] == -1){
      dlong start = C->diag->rowStarts[id];
      dlong end   = C->diag->rowStarts[id+1];

      // verify that all NBS are free
      int flag = 0;
      for(dlong j=start; j<end;j++) {
        if (states[C->diag->cols[j]]>-1) {
          flag=1;
          break;
        }
      }
      if (flag) continue; //skip

      // construct the aggregate
      states[id] = Agg_num;
      for(dlong j=start; j<end;j++)
        states[C->diag->cols[j]] = Agg_num;

      Agg_num++;
    }
  }

  dlong *aggCnts = (dlong *) calloc(N,sizeof(dlong));

  // count the number of nodes at each aggregate
  int maxNeighbors = 0;
  for (dlong i=0;i<N;i++) {
    if (states[i]>-1)
      aggCnts[states[i]]++;

    int Nneighbors = C->diag->rowStarts[i+1] - C->diag->rowStarts[i];
    maxNeighbors = (maxNeighbors > Nneighbors) ? maxNeighbors : Nneighbors;
  }

  dlong *neighborAgg = (dlong *) calloc(maxNeighbors, sizeof(dlong));

  // #pragma omp parallel for
  for(dlong i=0; i<N; i++){

    dlong id = V[i].index;

    //visit unaggregated nodes
    if (states[id] == -1){

      dlong start = C->diag->rowStarts[id];
      dlong end   = C->diag->rowStarts[id+1];
      int Nneighbors = end - start;

      if (Nneighbors > 1){    // at least one neighbor
        dlong Agg;
        int maxConnections=0;

        int Nactive = Nneighbors;

        //there must be at least one neighbor that has been aggregated, else this node would have been aggregated already
        for (int j=start;j<end;j++) {
          if (states[C->diag->cols[j]]==-1) {
            neighborAgg[j-start] = -1; //ignore unaggregated neighbors
            Nactive--;
          } else {
            // flag the aggregated neighbors, we'll visit them one at a time.
            neighborAgg[j-start] = states[C->diag->cols[j]];
          }
        }

        int pos = 0;
        //find the first unflagged aggregated neighbor
        while (Nactive>0) {
          //find the next aggregate id to check
          while (neighborAgg[pos]==-1) pos++;

          dlong curAgg = neighborAgg[pos]; //current aggregate under consideration

          //count the number of neighbours aggregating to this aggregate (unflag them as we go)
          int cnt = 0;
          for (int j=pos;j<Nneighbors;j++) {
            if (neighborAgg[j]==curAgg) {
              cnt++;
              Nactive--;
              neighborAgg[j]=-1;
            }
          }

          if (cnt > maxConnections) {
            //more connections to this aggregate, aggregate here
            maxConnections = cnt;
            Agg = curAgg;
          } else if (cnt == maxConnections) {
            //tied for connections, pick to the smaller aggregate
            if (aggCnts[curAgg]<aggCnts[Agg])
              Agg = curAgg;
          }
        }
        states[id] = Agg;
        aggCnts[Agg]++;

      } else { //isolated node

        states[id] = Agg_num;  // becomes a new aggregate
        aggCnts[Agg_num]++;
        Agg_num++;

      }
    }
  }

  dlong *gNumAggs = (dlong *) calloc(size,sizeof(dlong));

  // count the coarse nodes/aggregates in each rannk
  MPI_Allgather(&Agg_num,1,MPI_DLONG,gNumAggs,1,MPI_DLONG,A->comm);

  globalAggStarts[0] = 0;
  for (int r=0;r<size;r++)
    globalAggStarts[r+1] = globalAggStarts[r] + gNumAggs[r];

  // enumerate the coarse nodes/aggregates
  for(dlong i=0; i<N; i++)
    FineToCoarse[i] = globalAggStarts[rank] + states[i];

  // share results
  for (dlong n=N;n<M;n++) FineToCoarse[n] = 0;
  ogsGatherScatter(FineToCoarse, ogsHlong,  ogsAdd, A->ogs);

  free(states);
  free(V);
  free(aggCnts);

  delete C;
}

} //namespace parAlmond

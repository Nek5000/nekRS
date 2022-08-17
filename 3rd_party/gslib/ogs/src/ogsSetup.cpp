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

#include "ogstypes.h"
#include "ogs.hpp"
#include "ogsKernels.hpp"
#include "ogsInterface.h"

typedef struct{

  dlong localId;    // local node id
  hlong baseId;     // original global index

  dlong newId;         // new global id
  int owned;

}parallelNode_t;

// compare on baseId then by localId
int compareBaseId(const void *a, const void *b){

  parallelNode_t *fa = (parallelNode_t*) a;
  parallelNode_t *fb = (parallelNode_t*) b;

  if(abs(fa->baseId) < abs(fb->baseId)) return -1; //group by abs(baseId)
  if(abs(fa->baseId) > abs(fb->baseId)) return +1;

  if(fa->localId < fb->localId) return -1; //sort by local id
  if(fa->localId > fb->localId) return +1;

  return 0;
}

// compare on haloOwned then localId
int compareLocalId(const void *a, const void *b){

  parallelNode_t *fa = (parallelNode_t*) a;
  parallelNode_t *fb = (parallelNode_t*) b;

  if(fa->localId < fb->localId) return -1;
  if(fa->localId > fb->localId) return +1;

  return 0;
}

void setupRowBlocks(ogs_t *ogs, occa::device &device) 
{
  dlong blockSum=0;
  ogs->NrowBlocks=0;
  if (ogs->NlocalGather) ogs->NrowBlocks++;
  for (dlong i=0;i<ogs->NlocalGather;i++) {
    dlong rowSize = ogs->localGatherOffsets[i+1]-ogs->localGatherOffsets[i];

    if (rowSize > ogs::gatherNodesPerBlock) {
      //this row is pathalogically big. We can't currently run this
      std::cout << "Multiplicity of global node id: " << i << "in ogsSetup is too large.";
      exit(1);
    }

    if (blockSum+rowSize > ogs::gatherNodesPerBlock) { //adding this row will exceed the nnz per block
      ogs->NrowBlocks++; //count the previous block
      blockSum=rowSize; //start a new row block
    } else {
      blockSum+=rowSize; //add this row to the block
    }
  }

  dlong* blockRowStarts  = (dlong*) calloc(ogs->NrowBlocks+1,sizeof(dlong));

  blockSum=0;
  ogs->NrowBlocks=0;
  if (ogs->NlocalGather) ogs->NrowBlocks++;
  for (dlong i=0;i<ogs->NlocalGather;i++) {
    dlong rowSize = ogs->localGatherOffsets[i+1]-ogs->localGatherOffsets[i];

    if (blockSum+rowSize > ogs::gatherNodesPerBlock) { //adding this row will exceed the nnz per block
      blockRowStarts[ogs->NrowBlocks++] = i; //mark the previous block
      blockSum=rowSize; //start a new row block
    } else {
      blockSum+=rowSize; //add this row to the block
    }
  }
  blockRowStarts[ogs->NrowBlocks] = ogs->NlocalGather;
  ogs->o_blockRowStarts = device.malloc((ogs->NrowBlocks+1)*sizeof(dlong), blockRowStarts);
  free(blockRowStarts);
}

ogs_t *ogsSetup(dlong N, hlong *ids, MPI_Comm &comm,
                int verbose, occa::device device){

  //  ogs_t *ogs = (ogs_t*) calloc(1, sizeof(ogs_t));
  ogs_t *ogs = new ogs_t[1];

  ogs->NhaloGather = 0;
  ogs->Ngather = 0;
  ogs->Nlocal = 0;
  ogs->NlocalGather = 0;
  ogs->Nhalo = 0;
  ogs->NhaloGather = 0;
  ogs->NownedHalo = 0;

  ogs->localGatherOffsets = NULL;
  ogs->localGatherIds = NULL;

  ogs->haloGatherOffsets = NULL;
  ogs->haloGatherIds = NULL;

  ogs->hostGsh = NULL;
  ogs->haloGshSym = NULL;
  ogs->haloGshNonSym = NULL;

  ogs->invDegree = NULL;
  ogs->gatherInvDegree = NULL;
  
  ogs::Nrefs++;

  ogs->N = N;
  ogs->comm = comm;
  
  int rank, size;
  MPI_Comm_rank(ogs->comm, &rank);
  MPI_Comm_size(ogs->comm, &size);

  //make a host gs handle (calls gslib)
  ogs->hostGsh = ogsHostSetup(comm, N, ids, 0, 1);

  //use the host gs to find what nodes are local to this rank
  int *minRank = (int *) calloc(N,sizeof(int));
  int *maxRank = (int *) calloc(N,sizeof(int));
  hlong *flagIds   = (hlong *) calloc(N,sizeof(hlong));
  for (dlong i=0;i<N;i++) {
    minRank[i] = rank;
    maxRank[i] = rank;
    flagIds[i] = ids[i];
  }

  ogsHostGatherScatter(minRank, ogsInt, ogsMin, ogs->hostGsh); //minRank[n] contains the smallest rank taking part in the gather of node n
  ogsHostGatherScatter(maxRank, ogsInt, ogsMax, ogs->hostGsh); //maxRank[n] contains the largest rank taking part in the gather of node n
  ogsGsUnique(flagIds, N, comm); //one unique node in each group is 'flagged' kept positive while others are turned negative.

  //count local and halo nodes
  ogs->Nlocal=0; ogs->Nhalo=0; ogs->NownedHalo=0;
  for (dlong i=0;i<N;i++) {
    if (ids[i]==0) continue;

    if ((minRank[i]!=rank)||(maxRank[i]!=rank)) { // if I am not involved
      ogs->Nhalo++;
      if (flagIds[i]>0) ogs->NownedHalo++;
    } else {
      ogs->Nlocal++;
    }
  }

  //set up the local gatherScatter
  parallelNode_t *localNodes;
  
  if (ogs->Nlocal) {
    localNodes = (parallelNode_t*) calloc(ogs->Nlocal,sizeof(parallelNode_t));

    dlong cnt=0;
    for (dlong i=0;i<N;i++) {
      if (ids[i]==0) continue;

      if ((minRank[i]==rank)&&(maxRank[i]==rank)) {
        localNodes[cnt].localId = i;
        localNodes[cnt].baseId  = ids[i];
        localNodes[cnt].owned   = 0;
        cnt++;
      }
    }

    // sort based on base ids then local id
    qsort(localNodes, ogs->Nlocal, sizeof(parallelNode_t), compareBaseId);

    ogs->NlocalGather = 0;
    localNodes[0].newId = 0;
    localNodes[0].owned = 1;
    for (dlong i=1;i<ogs->Nlocal;i++) {
      int s = 0;
      if (localNodes[i].baseId!=localNodes[i-1].baseId) {
        ogs->NlocalGather++;
        s = 1;
      }
      localNodes[i].newId = ogs->NlocalGather;
      localNodes[i].owned = s;
    }
    ogs->NlocalGather++;

    // sort based on local ids
    qsort(localNodes, ogs->Nlocal, sizeof(parallelNode_t), compareLocalId);

    //tally up how many nodes are being gathered to each gatherNode and
    //  map to a local ordering
    dlong *localGatherCounts = (dlong*) calloc(ogs->NlocalGather,sizeof(dlong));
    dlong *localGatherMap    = (dlong*) calloc(ogs->NlocalGather,sizeof(dlong));
    cnt = 0;
    for (dlong i=0;i<ogs->Nlocal;i++) {
      dlong newId = localNodes[i].newId; //get the ordered id

      if (localNodes[i].owned)
        localGatherMap[newId] = cnt++; //record a new index if this is a new gatherNode

      localNodes[i].newId = localGatherMap[newId]; //reorder
      localGatherCounts[localGatherMap[newId]]++;  //tally
    }
    free(localGatherMap);

    ogs->localGatherOffsets = (dlong*) calloc(ogs->NlocalGather+1,sizeof(dlong));
    for (dlong i=0;i<ogs->NlocalGather;i++) {
      ogs->localGatherOffsets[i+1] = ogs->localGatherOffsets[i] + localGatherCounts[i];
      localGatherCounts[i] = 0;
    }

    ogs->localGatherIds = (dlong*) calloc(ogs->Nlocal,sizeof(dlong));
    for (dlong i=0;i<ogs->Nlocal;i++) {
      dlong gatherId = localNodes[i].newId;
      dlong offset = ogs->localGatherOffsets[gatherId];
      int index  = localGatherCounts[gatherId];

      ogs->localGatherIds[offset+index] = localNodes[i].localId;
      localGatherCounts[gatherId]++;
    }
    free(localGatherCounts);

    ogs->o_localGatherOffsets = device.malloc((ogs->NlocalGather+1)*sizeof(dlong), ogs->localGatherOffsets);
    ogs->o_localGatherIds     = device.malloc((ogs->Nlocal)*sizeof(dlong), ogs->localGatherIds);

    free(localNodes);
  }

  //set up the halo gatherScatter
  parallelNode_t *haloNodes;

    haloNodes = (parallelNode_t*) calloc(ogs->Nhalo+1,sizeof(parallelNode_t));

    dlong cnt=0;
    for (dlong i=0;i<N;i++) {
      if (ids[i]==0) continue;

      if ((minRank[i]!=rank)||(maxRank[i]!=rank)) {
        haloNodes[cnt].localId = i;
        haloNodes[cnt].baseId  = flagIds[i];
        haloNodes[cnt].owned   = 0;
        cnt++;
      }
    }

    // sort based on base ids then local id
    if(ogs->Nhalo){
      qsort(haloNodes, ogs->Nhalo, sizeof(parallelNode_t), compareBaseId);
      
      //move the flagged node to the lowest local index if present
      cnt = 0;
      ogs->NhaloGather=0;
      haloNodes[0].newId = 0;
      haloNodes[0].owned = 1;
      
      for (dlong i=1;i<ogs->Nhalo;i++) {
	int s = 0;
	if (abs(haloNodes[i].baseId)!=abs(haloNodes[i-1].baseId)) { //new gather node
	  s = 1;
	  cnt = i;
	  ogs->NhaloGather++;
	}
	
	haloNodes[i].owned = s;
	haloNodes[i].newId = ogs->NhaloGather;
	if (haloNodes[i].baseId>0) {
	  haloNodes[i].baseId   = -abs(haloNodes[i].baseId);
	  haloNodes[cnt].baseId =  abs(haloNodes[cnt].baseId);
	}
    }
      ogs->NhaloGather++;
      
      // sort based on local ids
      qsort(haloNodes, ogs->Nhalo, sizeof(parallelNode_t), compareLocalId);
    }
    
    //tally up how many nodes are being gathered to each gatherNode and
    //  map to a local ordering
    dlong *haloGatherCounts = (dlong*) calloc(ogs->NhaloGather+1,sizeof(dlong));
    dlong *haloGatherMap    = (dlong*) calloc(ogs->NhaloGather+1,sizeof(dlong));
    hlong *symIds    = (hlong *) calloc(ogs->NhaloGather+1,sizeof(hlong));
    hlong *nonSymIds = (hlong *) calloc(ogs->NhaloGather+1,sizeof(hlong));

    cnt = 0;
    dlong cnt2 = ogs->NownedHalo;
    for (dlong i=0;i<ogs->Nhalo;i++) {
      dlong newId = haloNodes[i].newId; //get the ordered id

      if (haloNodes[i].owned) {
        dlong c;
        if (haloNodes[i].baseId>0)
          c = cnt++;
        else
          c = cnt2++;

        symIds[c]    = abs(haloNodes[i].baseId); //record the base id
        nonSymIds[c] = haloNodes[i].baseId;      //record the base id
        haloGatherMap[newId] = c; //record a new index if this is a new gatherNode
      }

      haloNodes[i].newId = haloGatherMap[newId];  //reorder
      haloGatherCounts[haloGatherMap[newId]]++;  //tally
    }
    free(haloGatherMap);

    ogs->haloGatherOffsets = (dlong*) calloc(ogs->NhaloGather+1,sizeof(dlong));
    for (dlong i=0;i<ogs->NhaloGather;i++) {
      ogs->haloGatherOffsets[i+1] = ogs->haloGatherOffsets[i] + haloGatherCounts[i];
      haloGatherCounts[i] = 0;
    }

    ogs->haloGatherIds = (dlong*) calloc(ogs->Nhalo+1,sizeof(dlong));
    for (dlong i=0;i<ogs->Nhalo;i++) {
      dlong gatherId = haloNodes[i].newId;
      dlong offset = ogs->haloGatherOffsets[gatherId];
      int index  = haloGatherCounts[gatherId];

      ogs->haloGatherIds[offset+index] = haloNodes[i].localId;
      haloGatherCounts[gatherId]++;
    }
    free(haloGatherCounts);

    ogs->o_haloGatherOffsets = device.malloc((ogs->NhaloGather+1)*sizeof(dlong), ogs->haloGatherOffsets);
    ogs->o_haloGatherIds     = device.malloc((ogs->Nhalo+1)*sizeof(dlong), ogs->haloGatherIds);

    //make a host gs handle
    ogs->haloGshSym    = ogsHostSetup(comm, ogs->NhaloGather, symIds,    0,0);
    ogs->haloGshNonSym = ogsHostSetup(comm, ogs->NhaloGather, nonSymIds, 0,0);

    free(symIds); free(nonSymIds);
    free(haloNodes);

  free(minRank); free(maxRank); free(flagIds);

  //total number of owned gathered nodes
  ogs->Ngather = ogs->NlocalGather+ogs->NownedHalo;

  ogs->device = device;

  // build degree vectors
  ogs->invDegree = (dfloat*) calloc(N, sizeof(dfloat));
  ogs->gatherInvDegree = (dfloat*) calloc(ogs->Ngather, sizeof(dfloat));
  for(dlong n=0;n<N;++n) ogs->invDegree[n] = 1;

  ogs->o_invDegree = device.malloc(N*sizeof(dfloat), ogs->invDegree);
  ogs->o_gatherInvDegree = device.malloc(ogs->Ngather*sizeof(dfloat), ogs->gatherInvDegree);

  ogsGather(ogs->o_gatherInvDegree, ogs->o_invDegree, ogsDfloat, ogsAdd, ogs);

  if(ogs->Ngather)
    ogs->o_gatherInvDegree.copyTo(ogs->gatherInvDegree);

  ogsScatter(ogs->o_invDegree, ogs->o_gatherInvDegree, ogsDfloat, ogsAdd, ogs);

  if (N) ogs->o_invDegree.copyTo(ogs->invDegree);

  for(dlong n=0;n<ogs->N;++n)
    ogs->invDegree[n] = 1./ogs->invDegree[n];

  for(dlong n=0;n<ogs->Ngather;++n)
    ogs->gatherInvDegree[n] = 1./ogs->gatherInvDegree[n];

  if(ogs->Ngather)
    ogs->o_gatherInvDegree.copyFrom(ogs->gatherInvDegree);

  if(ogs->N)
    ogs->o_invDegree.copyFrom(ogs->invDegree);

  setupRowBlocks(ogs, device);

  return ogs;
}


void ogsFree(ogs_t *ogs) {

  if (ogs->Nlocal) {
    free(ogs->localGatherOffsets);
    free(ogs->localGatherIds);
    ogs->o_localGatherOffsets.free();
    ogs->o_localGatherIds.free();
  }

  if (ogs->Nhalo) {
    free(ogs->haloGatherOffsets);
    free(ogs->haloGatherIds);
    ogs->o_haloGatherOffsets.free();
    ogs->o_haloGatherIds.free();
    ogsHostFree(ogs->haloGshSym);
    ogsHostFree(ogs->haloGshNonSym);
  }

  if (ogs->N) {
    free(ogs->invDegree);
    ogs->o_invDegree.free();
    ogsHostFree(ogs->hostGsh);
  }

  if (ogs->Ngather) {
    free(ogs->gatherInvDegree);
    ogs->o_gatherInvDegree.free();
  }

  delete[] ogs;

  ogs::Nrefs--;
  if (!ogs::Nrefs) ogs::freeKernels();

}

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

//scratch space
size_t scratchSpaceBytes=0;
void *scratch=NULL;
occa::memory o_scratch;

size_t pinnedScratchSpaceBytes=0;
void *pinnedScratch=NULL;

occa::memory h_pinnedScratch;
occa::memory o_pinnedScratch;

size_t reductionScratchBytes=0;
void *reductionScratch=NULL;

occa::memory h_reductionScratch;
occa::memory o_reductionScratch;

 void *parAlmondHostMallocPinned(occa::device &device, size_t size, void *source, occa::memory &mem, occa::memory &h_mem){
   
#if 0
    
  mem = device.malloc(size, source);
  
  h_mem = device.mappedAlloc(size, source);
  
  void *ptr = h_mem.getMappedPointer();

#endif
  
  
  occa::properties props;
  props["mapped"] = true;
  
  if(source!=NULL)
    mem =  device.malloc(size, source);
  else
    mem =  device.malloc(size);

  h_mem =  device.malloc(size, props);
  
  void *ptr = h_mem.ptr(props);

  return ptr;

}
  
void allocateScratchSpace(size_t requiredBytes, occa::device device) {

  if (scratchSpaceBytes<requiredBytes) {
    if (scratchSpaceBytes!=0) {
      free(scratch);
      o_scratch.free();
    }
    scratch   = malloc(requiredBytes);
    memset(scratch, 0, requiredBytes);
    o_scratch = device.malloc(requiredBytes, scratch);
    scratchSpaceBytes = requiredBytes;
  }
  if (reductionScratchBytes==0) {
    reductionScratchBytes = 3*NBLOCKS*sizeof(dfloat);
    //    o_reductionScratch = device.mappedAlloc(reductionScratchBytes);
    //    reductionScratch = o_reductionScratch.getMappedPointer();
    reductionScratch = parAlmondHostMallocPinned(device, reductionScratchBytes, NULL, o_reductionScratch, h_reductionScratch);
  }
}

void freeScratchSpace() {
  if (scratchSpaceBytes!=0) {
    free(scratch);
    o_scratch.free();
  }
  scratchSpaceBytes=0;

  if (reductionScratchBytes!=0) {
    reductionScratchBytes = 0;
    o_reductionScratch.free();
  }
}

void allocatePinnedScratchSpace(size_t requiredBytes, occa::device device) {

  if (pinnedScratchSpaceBytes<requiredBytes) {
    if (pinnedScratchSpaceBytes!=0) {
      o_pinnedScratch.free();
    }
    //    o_pinnedScratch = device.mappedAlloc(requiredBytes);
    //    pinnedScratch = o_pinnedScratch.getMappedPointer();
    pinnedScratch = parAlmondHostMallocPinned(device, requiredBytes, NULL, o_pinnedScratch, h_pinnedScratch);
    pinnedScratchSpaceBytes = requiredBytes;
  }
}

void freePinnedScratchSpace() {
  if (pinnedScratchSpaceBytes!=0) {
    free(pinnedScratch);
    o_pinnedScratch.free();
  }
  pinnedScratchSpaceBytes=0;
}

// compare on global indices
int CompareGlobalId(const void *a, const void *b){

  parallelId_t *fa = (parallelId_t*) a;
  parallelId_t *fb = (parallelId_t*) b;

  if(fa->globalId < fb->globalId) return -1;
  if(fa->globalId > fb->globalId) return +1;

  if(fa->localId < fb->localId) return -1;
  if(fa->localId > fb->localId) return +1;

  return 0;
}

// compare on local indices
int CompareLocalId(const void *a, const void *b){

  parallelId_t *fa = (parallelId_t*) a;
  parallelId_t *fb = (parallelId_t*) b;

  if(fa->localId < fb->localId) return -1;
  if(fa->localId > fb->localId) return +1;

  if(fa->globalId < fb->globalId) return -1;
  if(fa->globalId > fb->globalId) return +1;

  return 0;
}

bool customLess(int smax, dfloat rmax, hlong imax, int s, dfloat r, hlong i){

  if(s > smax) return true;
  if(smax > s) return false;

  if(r > rmax) return true;
  if(rmax > r) return false;

  if(i > imax) return true;
  if(i < imax) return false;

  return false;
}

int compareOwner(const void *a, const void *b){
  parallelAggregate_t *pa = (parallelAggregate_t *) a;
  parallelAggregate_t *pb = (parallelAggregate_t *) b;

  if (pa->ownerRank < pb->ownerRank) return -1;
  if (pa->ownerRank > pb->ownerRank) return +1;

  return 0;
}

int compareAgg(const void *a, const void *b){
  parallelAggregate_t *pa = (parallelAggregate_t *) a;
  parallelAggregate_t *pb = (parallelAggregate_t *) b;

  if (pa->coarseId < pb->coarseId) return -1;
  if (pa->coarseId > pb->coarseId) return +1;

  if (pa->originRank < pb->originRank) return -1;
  if (pa->originRank > pb->originRank) return +1;

  return 0;
}

int compareOrigin(const void *a, const void *b){
  parallelAggregate_t *pa = (parallelAggregate_t *) a;
  parallelAggregate_t *pb = (parallelAggregate_t *) b;

  if (pa->originRank < pb->originRank) return -1;
  if (pa->originRank > pb->originRank) return +1;

  return 0;
}

int compareNonZeroByRow(const void *a, const void *b){
  nonzero_t *pa = (nonzero_t *) a;
  nonzero_t *pb = (nonzero_t *) b;

  if (pa->row < pb->row) return -1;
  if (pa->row > pb->row) return +1;

  if (pa->col < pb->col) return -1;
  if (pa->col > pb->col) return +1;

  return 0;
};


void matrixInverse(int N, dfloat *A){
  int lwork = N*N;
  int info;

  // compute inverse mass matrix
  double *tmpInvA = (double*) calloc(N*N, sizeof(double));

  int *ipiv = (int*) calloc(N, sizeof(int));
  double *work = (double*) calloc(lwork, sizeof(double));

  for(int n=0;n<N*N;++n){
    tmpInvA[n] = A[n];
  }

  dgetrf_ (&N, &N, tmpInvA, &N, ipiv, &info);
  dgetri_ (&N, tmpInvA, &N, ipiv, work, &lwork, &info);

  if(info)
    printf("inv: dgetrf/dgetri reports info = %d when inverting matrix\n", info);

  for(int n=0;n<N*N;++n)
    A[n] = tmpInvA[n];

  free(work);
  free(ipiv);
  free(tmpInvA);
}

void eig(const int Nrows, double *A, double *WR, double *WI){

  if(Nrows){
    int NB  = 256;
    char JOBVL  = 'V';
    char JOBVR  = 'V';
    int     N = Nrows;
    int   LDA = Nrows;
    int  LWORK  = (NB+2)*N;

    double *WORK  = new double[LWORK];
    double *VL  = new double[Nrows*Nrows];
    double *VR  = new double[Nrows*Nrows];

    int INFO = -999;

    dgeev_ (&JOBVL, &JOBVR, &N, A, &LDA, WR, WI,
      VL, &LDA, VR, &LDA, WORK, &LWORK, &INFO);


    // assert(INFO == 0);

    delete [] VL;
    delete [] VR;
    delete [] WORK;
  }
}


} //namespace parAlmond

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

#ifndef PARALMOND_UTILS_HPP
#define PARALMOND_UTILS_HPP

namespace parAlmond {

//scratch space
extern size_t scratchSpaceBytes;
extern void *scratch;
extern occa::memory o_scratch;

extern size_t pinnedScratchSpaceBytes;
extern void *pinnedScratch;
extern occa::memory o_pinnedScratch;

extern size_t reductionScratchBytes;
extern void *reductionScratch;
extern occa::memory o_reductionScratch;

void allocateScratchSpace(size_t requiredBytes, occa::device device);
void allocatePinnedScratchSpace(size_t requiredBytes, occa::device device);
void freeScratchSpace();
void freePinnedScratchSpace();

typedef struct {

  dlong localId;
  hlong globalId;

  dlong newId;

} parallelId_t;


typedef struct {

  dlong fineId;
  hlong coarseId;
  hlong newCoarseId;

  int originRank;
  int ownerRank;

} parallelAggregate_t;


typedef struct {

  hlong row;
  hlong col;
  dfloat val;

} nonzero_t;


int CompareGlobalId(const void *a, const void *b);
int CompareLocalId(const void *a, const void *b);

int compareOwner(const void *a, const void *b);
int compareAgg(const void *a, const void *b);
int compareOrigin(const void *a, const void *b);

int compareNonZeroByRow(const void *a, const void *b);

bool customLess(int smax, dfloat rmax, hlong imax, int s, dfloat r, hlong i);

extern "C"{
  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
  void dgeev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, double *WR, double *WI,
  double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO );
}

void eig(const int Nrows, double *A, double *WR, double *WI);

void matrixInverse(int N, dfloat *A);

} //namespace parAlmond

#endif
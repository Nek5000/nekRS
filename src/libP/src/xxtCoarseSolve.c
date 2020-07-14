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


#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "gslib.h"

typedef struct {
  struct crs_data *A;

  uint numLocalRows;
  uint nnz;

  ulong *rowIds;
  uint   *Ai; // row coordinates of non-zeros (local indexes)
  uint   *Aj; // column coordinates  of non-zeros (local indexes)
  double *Avals; // values of non-zeros

  double *x;
  double *rhs;

  char* intType;
  char* dfloatType;

} crs_t;

void * xxtSetup(uint  numLocalRows, 
                void* rowIds,
                uint  nnz, 
                void* Ai,
                void* Aj,
                void* Avals,
                int   nullSpace,
                const char* intType, 
                const char* dfloatType) {

  int n;
  struct comm com;
  crs_t *crsA = (crs_t*) calloc(1, sizeof(crs_t));

  comm_init(&com,(comm_ext) MPI_COMM_WORLD);

  crsA->numLocalRows = numLocalRows;
  crsA->nnz = nnz;

  if (!strcmp(dfloatType,"float")) { //float
    crsA->Avals = (double *) malloc(nnz*sizeof(double));
    crsA->x     = (double *) malloc(numLocalRows*sizeof(double));
    crsA->rhs   = (double *) malloc(numLocalRows*sizeof(double));
    float *AvalsFloat = (float *) Avals;
    for (n=0;n<nnz;n++) crsA->Avals[n] = (double) AvalsFloat[n];
  } else { //double
    crsA->Avals = (double*) Avals;
  }

  if (!strcmp(intType,"int")) { //int
    crsA->Ai = (uint*) calloc(nnz, sizeof(uint));
    crsA->Aj = (uint*) calloc(nnz, sizeof(uint));      
    
    crsA->rowIds = (ulong*) malloc(numLocalRows*sizeof(ulong));
    int *rowIdsInt = (int*) rowIds;
    for (n=0;n<numLocalRows;n++) crsA->rowIds[n] = (ulong) rowIdsInt[n];
    
    for(n=0;n<nnz;++n){
      crsA->Ai[n] = ((int*)Ai)[n];
      crsA->Aj[n] = ((int*)Aj)[n];
    }
    
  } else { //long
    printf("Exiting due to use of ulong in intType %s\n", intType);
    exit(-1);
#if 0
    crsA->Ai = (ulong*) calloc(nnz, sizeof(ulong));
    crsA->Aj = (ulong*) calloc(nnz, sizeof(ulong));  

    crsA->rowIds = (ulong *) rowIds;
    
    for(n=0;n<nnz;++n){
      crsA->Ai[n] = ((ulong*)Ai)[n];
      crsA->Aj[n] = ((ulong*)Aj)[n];
    }
#endif
  }

  crsA->intType = strdup(intType);
  crsA->dfloatType = strdup(dfloatType);

  crsA->A = crs_setup(crsA->numLocalRows,
		      crsA->rowIds,
		      crsA->nnz,
		      crsA->Ai,
		      crsA->Aj,
		      crsA->Avals,
		      nullSpace,
		      &com);

  crs_stats(crsA->A);

  return (void *) crsA;
}

void xxtSolve(void* x,
             void* A,
             void* rhs) {

  int n;
  float *xFloat;
  float *rhsFloat;
  
  crs_t *crsA = (crs_t *) A;

  if (!strcmp(crsA->dfloatType,"float")) {  
    xFloat   = (float *) x;
    rhsFloat = (float *) rhs;
    for (n=0;n<crsA->numLocalRows;n++) {
      crsA->x[n]   = (double) xFloat[n];
      crsA->rhs[n] = (double) rhsFloat[n];
    }
  } else {
    crsA->x   = (double*) x;
    crsA->rhs = (double*) rhs;
  }

  crs_solve(crsA->x,crsA->A,crsA->rhs);
  
  if (!strcmp(crsA->dfloatType,"float")) {
    
    float *xFloat   = (float *) x;
    float *rhsFloat = (float *) rhs;
    for (n=0;n<crsA->numLocalRows;n++) {
      xFloat[n] = (float) crsA->x[n];
    }
  }
}

void xxtFree(void* A) {
  crs_t *crsA = (crs_t *) A;

  crs_free(crsA->A);

  if (!strcmp(crsA->dfloatType,"float")) { 
    free(crsA->Avals);
    free(crsA->x);  
    free(crsA->rhs);
  }

  if (!strcmp(crsA->intType,"int")) { 
    free(crsA->rowIds);
  }
}

#include "genmap-impl.h"
#include "genmap-io.h"

#include <stdlib.h>
#include <stdio.h>
//
// GenmapInit
//
int GenmapInit(GenmapHandle *h, GenmapCommExternal ce) {
  GenmapMalloc(1, h);
  GenmapHandle h_ = *h;

  GenmapCreateComm(&h_->global, ce);
  GenmapCreateComm(&h_->local, ce);

  h_->elementArray.ptr = NULL;
  h_->elementArray.n = (*h)->elementArray.max = 0;

  h_->dbgLevel = 0;
  h_->printStat = 0;

  GenmapMalloc(1,&h_->histogram);
  return 0;
}
//
// GenmapFinalize
//
int GenmapFinalize(GenmapHandle h) {
  if(GenmapGetGlobalComm(h))
    GenmapDestroyComm(GenmapGetGlobalComm(h));
  if(GenmapGetLocalComm(h))
    GenmapDestroyComm(h->local);

  array_free(&(h->elementArray));

  GenmapFree(h->histogram);

  GenmapFree(h);

  return 0;
}
//
// GenmapMalloc, Realloc, Calloc and Free
//
int GenmapMallocArray(size_t n, size_t unit, void *p) {
  int ierr = posix_memalign((void **)p, GENMAP_ALIGN, n * unit);
  if(ierr)
    printf("GenmapMallocArray Failed: %s:%d\n", __FILE__, __LINE__);
  return ierr;
}

int GenmapCallocArray(size_t n, size_t unit, void *p) {
  int ierr = 0;
  *(void **)p = calloc(n, unit);
  if(n && unit && !*(void **)p) {
    ierr = 1;
    printf("GenmapCallocArray Failed: %s:%d\n", __FILE__, __LINE__);
  }
  return ierr;
}

int GenmapReallocArray(size_t n, size_t unit, void *p) {
  int ierr = 0;
  *(void **)p = realloc(*(void **)p, n * unit);
  if(n && unit && !*(void **)p) {
    ierr = 1;
    printf("GenmapReallocArray Failed: %s:%d\n", __FILE__, __LINE__);
  }
  return ierr;
}

int GenmapFree(void *p) {
  free(p);
  p = NULL;
  return 0;
}

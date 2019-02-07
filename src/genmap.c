#define _POSIX_C_SOURCE 200112
#include <genmap-impl.h>
#include <genmap-io.h>

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <time.h>
//
// Genmap Readers (FEM meshes, .map files, etc.)
//
static struct {
  char name[GENMAP_READER_LEN];
  int (*Create)(GenmapHandle h);
} GenmapReaders[GENMAP_MAX_READERS];
static size_t GenmapNumReaders = 0;
static int GenmapReadersRegistered = 0;
//
// Register readers -- each reader should call this
//
int GenmapRegisterReader(char *name, int (*Create)(GenmapHandle h)) {
  if(GenmapNumReaders >= sizeof(GenmapReaders) / sizeof(GenmapReaders[0])) {
    //TODO: GenmapError
    printf("Error: Too many readers.\n");
  }

  strncpy(GenmapReaders[GenmapNumReaders].name, name, GENMAP_READER_LEN);
  GenmapReaders[GenmapNumReaders].Create = Create;
  GenmapNumReaders++;

  return 0;
}
//
// GenmapRegister
//
int GenmapRegister() {
  int ierr;
  ierr = GenmapRegisterReader("interface", GenmapCreateHandle_interface);
  return ierr;
}
//
// GenmapInit
//
int GenmapInit(GenmapHandle *h, GenmapCommExternal ce, char *reader) {
  if(GenmapReadersRegistered == 0) {
    GenmapRegister();
    GenmapReadersRegistered = 1;
  }

  char *registeredReader;
  size_t matchLen = 0, matchIdx = 0;

  size_t i, j;
  for(i = 0; i < GenmapNumReaders; i++) {
    registeredReader = GenmapReaders[i].name;
    for(j = 0; reader[j]
        && (reader[j] == registeredReader[j]); j++) {
      if(j > matchLen) {
        matchLen = j;
        matchIdx = i;
      }
    }
  }

  if(!matchLen) {
    //TODO: GenmapError
    printf("Error: Reader not found.\n");
  }

  GenmapMalloc(1, h);
  (*h)->Create = GenmapReaders[matchIdx].Create;
  GenmapCreateHandle(*h);

  GenmapCreateComm(&(*h)->global, ce);
  GenmapCreateComm(&(*h)->local, ce);
  GenmapCreateHeader(&(*h)->header);

  return 0;
}
//
// GenmapFinalize
//
int GenmapFinalize(GenmapHandle h) {
  if(h->global)
    GenmapDestroyComm(h->global);
  if(h->local)
    GenmapDestroyComm(h->local);

  GenmapDestroyHeader(h->header);

  array_free(&(h->elementArray));

  GenmapDestroyHandle(h);

  return 0;
}
//
// GenmapHandle
//
int GenmapCreateHandle(GenmapHandle h) {
  // Datastructures
  h->global = NULL;
  h->local = NULL;
  h->header = NULL;
  h->elementArray.ptr = NULL;
  h->elementArray.n = h->elementArray.max = 0;

  // Function pointers
  h->Id = GenmapId;
  h->Np = GenmapNp;
  h->Gop = GenmapGop;
  h->Ax = GenmapAx;
  h->AxInit = GenmapAxInit;

  h->Create(h);

  return 0;
}

int GenmapDestroyHandle(GenmapHandle h) {
  GenmapFree(h);
  return 0;
}
//
// GenmapHeader: Create, Destroy
//
int GenmapCreateHeader(GenmapHeader *h) {
  GenmapMalloc(1, h);

  return 0;
}

int GenmapDestroyHeader(GenmapHeader h) {
  GenmapFree(h);
  return 0;
}

int GenmapCreateComm(GenmapComm *c, GenmapCommExternal ce) {
  GenmapMalloc(1, c);
  comm_init(&(*c)->gsComm, ce);
  (*c)->verticesHandle = NULL;
  (*c)->laplacianWeights = NULL;
  buffer_init(&(*c)->buf, 1024);
  return 0;
}

int GenmapDestroyComm(GenmapComm c) {
  if(&c->buf)
    buffer_free(&c->buf);
  if(c->verticesHandle)
    gs_free(c->verticesHandle);
  if(c->laplacianWeights)
    GenmapFree(c->laplacianWeights);
  if(&c->gsComm)
    comm_free(&c->gsComm);
  GenmapFree(c);

  return 0;
}

int GenmapNp(GenmapComm c) {
  return (int) c->gsComm.np;
}

int GenmapId(GenmapComm c) {
  return (int) c->gsComm.id;
}

GenmapElements GenmapGetElements(GenmapHandle h) {
  return (GenmapElements) h->elementArray.ptr;
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

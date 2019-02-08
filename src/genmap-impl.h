#ifndef _GENMAP_IMPL_H_
#define _GENMAP_IMPL_H_

#include <stddef.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#ifdef GENMAP_DEBUG
#include <stdio.h>
#endif

#include "genmap.h"
//
// GenmapComm
//
struct GenmapComm_private {
  struct comm gsComm;
  struct gs_data *verticesHandle;
  GenmapScalar *laplacianWeights;
  buffer buf;
};
//
// GenmapElements
//
struct GenmapElement_private {
  GenmapScalar fiedler;
  GenmapLong globalId;
  GenmapLong vertices[8];
  GenmapInt proc;
};
//
// GenmapElements: Create, Destroy
//
int GenmapCreateElements(GenmapElements *e);
int GenmapDestroyElements(GenmapElements e);
GenmapElements GenmapGetElements_default(GenmapHandle h);
//
// Genmap_Handle
//
struct GenmapHandle_private {
  GenmapComm global;
  GenmapComm local;

  GenmapLong nel;
  GenmapLong npts;
  int nv;
  GenmapInt lelt;
  GenmapLong start;
  GenmapLong Nnodes;

  struct array elementArray;

  int (*Create)(GenmapHandle h);

  GenmapInt dbgLevel;
  GenmapInt printStat;
};
//
// GenmapHandle
//
int GenmapCreateHandle(GenmapHandle h);
int GenmapDestroyHandle(GenmapHandle h);
//
// GenmapVector
//
struct GenmapVector_private {
  GenmapInt size;
  GenmapScalar *data;
};
//
// Memory management routines
//
#define GenmapMalloc(n, p) GenmapMallocArray ((n), sizeof(**(p)), p)
#define GenmapCalloc(n, p) GenmapCallocArray ((n), sizeof(**(p)), p)
#define GenmapRealloc(n, p) GenmapReallocArray((n), sizeof(**(p)), p)

#endif

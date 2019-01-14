#ifndef _GENMAP_IMPL_H_
#define _GENMAP_IMPL_H_

#include <genmap.h>

#include <stddef.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#ifdef GENMAP_DEBUG
#include <stdio.h>
#endif
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
// File I/O
//
// Genmap Input File header
struct GenmapHeader_private {
  GenmapLong nel;
  GenmapLong npts;
  int nv;
  int ndim;
  GenmapInt lelt;
  GenmapLong start;
  GenmapLong Nnodes;
};
//
// GenmapHeader: Create, Destroy
//
int GenmapCreateHeader(GenmapHeader *h);
int GenmapDestroyHeader(GenmapHeader h);
// GenmapElements
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
  int (*Create)(GenmapHandle h);
  int (*Destroy)(GenmapHandle h);

  GenmapComm global;
  GenmapComm local;
  int (*CreateComm)(GenmapComm *c, GenmapCommExternal ce);
  int (*DestroyComm)(GenmapComm c);

  GenmapHeader header;
  int (*CreateHeader)(GenmapHeader *h);
  int (*DestroyHeader)(GenmapHeader h);

  struct array elementArray;
  GenmapElements(*GetElements)(GenmapHandle h);
  int (*CreateElements)(GenmapElements *e);
  int (*DestroyElements)(GenmapElements e);

  GenmapInt(*Np)(GenmapComm c);
  GenmapInt(*Id)(GenmapComm c);

  int (*Ax)(GenmapHandle h, GenmapComm c, GenmapVector u,
            GenmapVector weights, GenmapVector v);
  int (*AxInit)(GenmapHandle h, GenmapComm c, GenmapVector weights);

  int (*Gop)(GenmapComm c, void *v, GenmapInt size, GenmapDataType t,
             GenmapInt op);

  int (*Read)(GenmapHandle h, void *data);

  int (*Write)(GenmapHandle h, char *fileNameBase);
};

// GenmapHandle
int GenmapCreateHandle(GenmapHandle h);
int GenmapDestroyHandle(GenmapHandle h);
//
// Genmap_Vector
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

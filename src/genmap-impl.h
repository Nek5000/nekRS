#ifndef _GENMAP_IMPL_H_
#define _GENMAP_IMPL_H_

#include <stddef.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#ifdef GENMAP_DEBUG
#include <stdio.h>
#endif

#include <genmap.h>
#include <genmap-multigrid-precon.h>

#define GENMAP_FIEDLER 0
#define GENMAP_GLOBALID 1
#define GENMAP_PROC 2
#define GENMAP_ORIGIN 3

#define GENMAP_RCB_ELEMENT 0
#define GENMAP_RSB_ELEMENT 1

#define MAXDIM 3 /* Maximum dimension of the mesh */
#define MAXNV 8 /* Maximum number of vertices per element */

struct GenmapComm_private{
  struct comm gsc;
  /* Un-weighted Laplacian */
  struct gs_data *gsh;
  csr_mat M;
  /* Weighted Laplacian */
  struct gs_data *gsw;
  buffer buf;
  GenmapScalar *b;
};

/* parRCB internals */
struct rcb_element{
  unsigned char type;
  GenmapInt proc;
  GenmapInt origin;
  GenmapInt seq;
  GenmapLong globalId;
  GenmapScalar coord[MAXDIM];
};

void rcb_local(struct array *a,uint start,uint end,int ndim,buffer *buf);
int rcb_level(struct comm *c,struct array *a,int ndim);
int rcb(struct comm *ci,struct array *a,int ndim);

/* parRSB internals */
/* rsb_element should be a superset of rsb_element */
struct rsb_element{
  unsigned char type;
  GenmapInt proc;
  GenmapInt origin;
  GenmapInt seq;
  GenmapLong globalId;
  GenmapScalar coord[MAXDIM];
  GenmapLong vertices[8];
  GenmapInt part;
  GenmapULong globalId0;
  GenmapScalar fiedler;
};

int GenmapCreateElements(GenmapElements *e);
int GenmapDestroyElements(GenmapElements e);
GenmapElements GenmapGetElements_default(genmap_handle h);

struct genmap_handle_private {
  GenmapComm global;
  GenmapComm local;

  GenmapLong nel;
  GenmapLong Nnodes;
  GenmapLong start;
  int nv;

  GenmapScalar *weights;
  struct array *elements;
  struct crystal cr;

  int verbose_level;
  int print_stat;
};

int GenmapCreateHandle(genmap_handle h);
int GenmapDestroyHandle(genmap_handle h);

struct GenmapVector_private {
  GenmapInt size;
  GenmapScalar *data;
};

#define GenmapMalloc(n, p) GenmapMallocArray ((n), sizeof(**(p)), p)
#define GenmapCalloc(n, p) GenmapCallocArray ((n), sizeof(**(p)), p)
#define GenmapRealloc(n, p) GenmapReallocArray((n), sizeof(**(p)), p)

/* Genmap Metrics */
typedef enum{
  ALLREDUCE=0,
  AXISLEN,
  BINN1,
  BINN2,
  BISECT,
  COMMSPLIT,
  CRYSTAL,
  FIEDLER,
  FMG,
  GSOP,
  GSSETUP,
  GRAMMIAN,
  LANCZOS,
  LAPLACIANSETUP0,
  LAPLACIANSETUP1,
  LAPLACIAN,
  LOCALSORT,
  NCONN,
  NNEIGHBORS,
  NFIEDLER,
  NDISCON,
  NPROJECTPF,
  NRQI,
  PAIRWISE,
  PARSORT,
  PRECONSETUP,
  PRECONAX,
  PROJECTPF,
  PROJECT,
  RCB,
  RCBTRANSFER,
  RQI,
  SETPROC,
  UPDATEPROBE,
  VCYCLE
} metric;

void metric_init();
void metric_finalize();
void metric_acc(metric m,double count);
void metric_tic(struct comm *c,metric m);
void metric_toc(struct comm *c,metric m);
double metric_get_value(int level,metric m);
void metric_push_level();
uint metric_get_levels();
void metric_print(struct comm *c);

/* genCon */
typedef struct{
  GenmapULong sequenceId;
  int nNeighbors;
  GenmapULong elementId;
  GenmapULong vertexId;
  uint workProc;
} vertex;

/* Components */
sint is_disconnected(struct comm *c,struct gs_data *gsh,buffer *buf,
  uint nelt,uint nv);

/* Matrix inverse */
void matrix_inverse(int N,double *A);

#endif

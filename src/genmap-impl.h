#ifndef _GENMAP_IMPL_H_
#define _GENMAP_IMPL_H_

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#ifdef GENMAP_DEBUG
#include <stdio.h>
#endif

#include <genmap-multigrid-precon.h>
#include <genmap.h>

#define GENMAP_RCB_ELEMENT 0
#define GENMAP_RSB_ELEMENT 1

#define GENMAP_ALIGN 32

#define GENMAP_SP_TOL 1e-05
#define GENMAP_DP_TOL 1e-12
#define GENMAP_TOL GENMAP_DP_TOL

#define MAXDIM 3 /* Maximum dimension of the mesh */
#define MAXNV 8  /* Maximum number of vertices per element */

/* rcb_element is used for rcb and rib */
struct rcb_element {
  int type;
  GenmapInt proc;
  GenmapInt origin;
  GenmapInt seq;
  GenmapLong globalId;
  GenmapScalar coord[MAXDIM];
  GenmapScalar fiedler;
};

/* rsb_element should be a superset of rcb_element */
struct rsb_element {
  int type;
  GenmapInt proc;
  GenmapInt origin;
  GenmapInt seq;
  GenmapLong globalId;
  GenmapScalar coord[MAXDIM];
  GenmapScalar fiedler;
  GenmapLong vertices[8];
  GenmapInt part;
  GenmapULong globalId0;
};

int rcb(struct comm *ci, struct array *elements, int ndim, buffer *bfr);
int rib(struct comm *ci, struct array *elements, int ndim, buffer *bfr);

struct genmap_handle_private {
  genmap_comm global;
  genmap_comm local;

  GenmapLong nel;
  GenmapLong Nnodes;
  GenmapLong start;
  int nv;

  struct array *elements;
  struct crystal cr;

  /* Weighted Laplacian */
  struct gs_data *gsw;
  GenmapScalar *weights;
  buffer buf;

  /* Un-weighted Laplacian */
  struct gs_data *gsh;
  csr_mat M;
  GenmapScalar *b;

  parRSB_options *options;
};

struct genmap_vector_private {
  GenmapInt size;
  GenmapScalar *data;
};

int GenmapMallocArray(size_t n, size_t unit, void *p);
int GenmapCallocArray(size_t n, size_t unit, void *p);
int GenmapReallocArray(size_t n, size_t unit, void *p);
int GenmapFree(void *p);

#define GenmapMalloc(n, p) GenmapMallocArray((n), sizeof(**(p)), p)
#define GenmapCalloc(n, p) GenmapCallocArray((n), sizeof(**(p)), p)
#define GenmapRealloc(n, p) GenmapReallocArray((n), sizeof(**(p)), p)

/* Genmap Metrics */
typedef enum {
  RCB,
  WEIGHTEDLAPLACIANSETUP,
  FIEDLER,
  NFIEDLER,
  FIEDLERSORT,
  BISECTANDREPAIR,
  LANCZOS,
  NLANCZOS,
  WEIGHTEDLAPLACIAN,
  TQLI,
  LAPLACIANSETUP,
  FINDNBRS,
  CSRMATSETUP,
  CSRTOPSETUP,
  PRECONDSETUP,
  RQI,
  NRQI,
  PROJECT,
  NPROJECT,
  GRAMMIAN,
  LAPLACIAN,
  VCYCLE,
  END
} metric;

void metric_init();
void metric_acc(metric m, double count);
void metric_tic(struct comm *c, metric m);
void metric_toc(struct comm *c, metric m);
double metric_get_value(int level, metric m);
void metric_push_level();
uint metric_get_levels();
void metric_print(struct comm *c);
void metric_finalize();

/* genCon */
typedef struct {
  GenmapULong sequenceId;
  int nNeighbors;
  GenmapULong elementId;
  GenmapULong vertexId;
  uint workProc;
} vertex;

/* Components */
sint get_components(sint *component, struct rsb_element *elements,
                    struct comm *c, buffer *buf, uint nelt, uint nv);

void split_and_repair_partitions(genmap_handle h, struct comm *lc, int level);
/* Matrix inverse */
void matrix_inverse(int N, double *A);

/* Dump data */
int GenmapFiedlerDump(const char *fname, genmap_handle h, slong start,
                      struct comm *c);
int GenmapVectorDump(const char *fname, GenmapScalar *y, uint size,
                     struct comm *c);
int GenmapCentroidDump(const char *fname, genmap_handle h, sint g_id,
                       struct comm *c);
int GenmapElementIdDump(const char *fname, genmap_handle h, struct comm *c);

#endif

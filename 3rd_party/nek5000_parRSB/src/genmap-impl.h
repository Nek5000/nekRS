#ifndef _GENMAP_IMPL_H_
#define _GENMAP_IMPL_H_

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

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
  GenmapULong globalId;
  GenmapScalar coord[MAXDIM];
  GenmapScalar fiedler;
};

/* rsb_element should be a superset of rcb_element */
struct rsb_element {
  int type;
  GenmapInt proc;
  GenmapInt origin;
  GenmapInt seq;
  GenmapULong globalId;
  GenmapScalar coord[MAXDIM];
  GenmapScalar fiedler;
  GenmapLong vertices[8];
  GenmapInt part;
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

  parrsb_options *options;
  size_t elem_size;
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
  FIEDLER,
  FIEDLER_NITER,
  LANCZOS,
  LANCZOS_NITER,
  LANCZOS_TOL_FINAL,
  LANCZOS_TOL_TARGET,
  LAPLACIAN,
  LAPLACIAN_INIT,
  RQI,
  RQI_NITER,
  PROJECT,
  PROJECT_NITER,
  COMPONENTS,
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

/* Repair and balance */
sint get_components(sint *component, struct rsb_element *elements,
                    struct comm *c, buffer *buf, uint nelt, uint nv);

int repair_partitions(genmap_handle h, struct comm *tc, struct comm *lc,
                      int bin, struct comm *gc);
int balance_partitions(genmap_handle h, struct comm *lc, int bin,
                       struct comm *gc);

/* Matrix inverse */
void matrix_inverse(int N, double *A);

/* Dump data */
int GenmapFiedlerDump(const char *fname, genmap_handle h, struct comm *c);
int GenmapVectorDump(const char *fname, GenmapScalar *y, genmap_handle h,
                     struct comm *c);
int GenmapCentroidDump(const char *fname, genmap_handle h, sint g_id,
                       struct comm *c);
int GenmapElementIdDump(const char *fname, genmap_handle h, struct comm *c);

#endif

#ifndef _GENMAP_H_
#define _GENMAP_H_

#include "genmap-gslib.h"
#include "genmap-types.h"
#include <mpi.h>

#include <parRSB.h>

typedef struct comm *genmap_comm;
typedef struct genmap_handle_private *genmap_handle;
typedef struct genmap_vector_private *genmap_vector;

/* genmap_handle */
int genmap_init(genmap_handle *h, comm_ext ce, parRSB_options *options);

void *genmap_get_elements(genmap_handle h);
void genmap_set_elements(genmap_handle h, struct array *elements);

genmap_comm genmap_local_comm(genmap_handle h);
genmap_comm genmap_global_comm(genmap_handle h);

void genmap_set_nvertices(genmap_handle h, int nv);
int genmap_get_nvertices(genmap_handle h);
GenmapULong genmap_get_partition_nel(genmap_handle h);
void genmap_set_partition_nel(genmap_handle h, GenmapULong globalElements);

GenmapLong genmap_get_local_start_index(genmap_handle h);
void genmap_set_local_start_index(genmap_handle h, GenmapLong localStart);

GenmapInt genmap_get_nel(genmap_handle h);

int genmap_finalize(genmap_handle h);

/* genmap_comm */
void genmap_comm_scan(genmap_handle h, genmap_comm c);
int genmap_comm_size(genmap_comm c);
int genmap_comm_rank(genmap_comm c);
void genmap_comm_split(struct comm *old, int bin, int key, struct comm *new_);

/* genmap_vector */
int genmap_vector_create(genmap_vector *x, GenmapInt size);
int genmap_vector_create_zeros(genmap_vector *x, GenmapInt size);
int genmap_vector_copy(genmap_vector x, genmap_vector y);

int genmap_vector_scale(genmap_vector y, genmap_vector x, GenmapScalar alpha);
int genmap_vector_axpby(genmap_vector z, genmap_vector x, GenmapScalar alpha,
                        genmap_vector y, GenmapScalar beta);

GenmapScalar genmap_vector_dot(genmap_vector x, genmap_vector y);

int genmap_vector_ortho_one(struct comm *c, genmap_vector q1, GenmapULong n);

int GenmapDestroyVector(genmap_vector x);

/* Laplacian */
int GenmapInitLaplacian(genmap_handle h, struct comm *c);
int GenmapLaplacian(genmap_handle h, GenmapScalar *u, GenmapScalar *v);

int GenmapInitLaplacianWeighted(genmap_handle h, struct comm *c);
int GenmapLaplacianWeighted(genmap_handle h, GenmapScalar *u, GenmapScalar *v);

/* Eigen */
int GenmapTQLI(genmap_handle h, genmap_vector diag, genmap_vector upper,
               genmap_vector **eVec, genmap_vector *eVal);
int genmap_inverse_power(double *y, int N, double *A, int verbose);
int genmap_power(double *y, int N, double *A, int verbose);

/* Matrix inverse */
void matrix_inverse(int N, double *A);

/* Lanczos */
int GenmapLanczosLegendary(genmap_handle h, struct comm *c, genmap_vector f,
                           GenmapInt niter, genmap_vector **rr,
                           genmap_vector diag, genmap_vector upper);
int GenmapLanczos(genmap_handle h, struct comm *c, genmap_vector init,
                  GenmapInt iter, genmap_vector **q, genmap_vector alpha,
                  genmap_vector beta);

/* Fiedler */
int GenmapFiedlerLanczos(genmap_handle h, struct comm *c, int maxIter,
                         int global);
int GenmapFiedlerRQI(genmap_handle h, struct comm *c, int maxIter, int global);

/* RSB/RCB */
void genmap_load_balance(struct array *eList, uint nel, int nv, double *coord,
                         long long *vtx, struct crystal *cr, buffer *bfr);
int genmap_rsb(genmap_handle h);
int genmap_rcb(genmap_handle h);
int genmap_rib(genmap_handle h);

void genmap_restore_original(int *part, int *seq, struct crystal *cr,
                             struct array *eList, buffer *bfr);

/* Misc */
double GenmapGetMaxRss();
void GenmapPrintStack();

#endif

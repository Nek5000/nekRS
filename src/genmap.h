#ifndef _GENMAP_H_
#define _GENMAP_H_

#include "genmap-types.h"
#include "genmap-gslib.h"
#include <mpi.h>

#define GENMAP_SUM 0
#define GENMAP_MAX 1
#define GENMAP_MIN 2
#define GENMAP_MUL 3

#define GENMAP_ALIGN 32

#define GENMAP_SP_TOL 1e-05
#define GENMAP_DP_TOL 1e-12
#define GENMAP_TOL GENMAP_DP_TOL

#define GENMAP_READER_LEN 256
#define GENMAP_MAX_READERS 32

typedef MPI_Datatype GenmapDataType;
typedef struct GenmapComm_private *GenmapComm;
typedef struct genmap_handle_private *genmap_handle;
typedef struct GenmapVector_private *GenmapVector;
typedef struct rsb_element *GenmapElements;

int genmap_init(genmap_handle *h, comm_ext ce, int *options);
int genmap_finalize(genmap_handle h);

int GenmapMallocArray(size_t n, size_t unit, void *p);
int GenmapCallocArray(size_t n, size_t unit, void *p);
int GenmapReallocArray(size_t n, size_t unit, void *p);
int GenmapFree(void *p);

GenmapElements GenmapGetElements(genmap_handle h);

GenmapComm GenmapGetLocalComm(genmap_handle h);
void GenmapSetLocalComm(genmap_handle h, GenmapComm c);

GenmapComm GenmapGetGlobalComm(genmap_handle h);
void GenmapSetGlobalComm(genmap_handle h, GenmapComm c);

GenmapInt GenmapGetNLocalElements(genmap_handle h);
void GenmapSetArrayElements(genmap_handle h, struct array *localElements);

GenmapLong GenmapGetNGlobalElements(genmap_handle h);
void GenmapSetNGlobalElements(genmap_handle h, GenmapLong globalElements);

GenmapLong GenmapGetLocalStartIndex(genmap_handle h);
void GenmapSetLocalStartIndex(genmap_handle h, GenmapLong localStart);

int GenmapGetNVertices(genmap_handle h);
void GenmapSetNVertices(genmap_handle h, int nVertices);

void GenmapScan(genmap_handle h, GenmapComm c);

int GenmapCreateComm(GenmapComm *c, comm_ext ce);
int GenmapCommSize(GenmapComm c);
int GenmapCommRank(GenmapComm c);

int GenmapGop(GenmapComm c, void *v, GenmapInt size, GenmapDataType type,
              GenmapInt op);
int GenmapReduce(GenmapComm c, void *out, void *in, GenmapInt size,
                 GenmapDataType type, GenmapInt op);
int GenmapBcast(GenmapComm c, void *in, GenmapInt count, GenmapDataType type);

int GenmapDestroyComm(GenmapComm c);
void GenmapSplitComm(genmap_handle h, GenmapComm *c, int bin);
int GenmapCrystalInit(genmap_handle h, GenmapComm c);
int GenmapCrystalTransfer(genmap_handle h, int field);
int GenmapCrystalFinalize(genmap_handle h);

int GenmapRead(genmap_handle h, void *data);

int GenmapCreateVector(GenmapVector *x, GenmapInt size);
int GenmapSetVector(GenmapVector x, GenmapScalar *array);
int GenmapGetVector(GenmapVector x, GenmapScalar *array);

int GenmapCreateRandomVector(GenmapVector *x, GenmapInt size, GenmapInt seed);
int GenmapCreateOnesVector(GenmapVector *x, GenmapInt size);
int GenmapCreateZerosVector(GenmapVector *x, GenmapInt size);

int GenmapScaleVector(GenmapVector y, GenmapVector x, GenmapScalar alpha);
int GenmapAxpbyVector(GenmapVector z, GenmapVector x, GenmapScalar alpha,
                      GenmapVector y, GenmapScalar beta);

int GenmapVectorsEqual(GenmapVector x, GenmapVector y, GenmapScalar tol);
int GenmapCopyVector(GenmapVector x, GenmapVector y);
GenmapScalar GenmapDotVector(GenmapVector x, GenmapVector y);
GenmapScalar GenmapAbsMaxVector(GenmapVector x);
GenmapScalar GenmapMaxVector(GenmapVector x);
GenmapScalar GenmapAbsMinVector(GenmapVector x);
GenmapScalar GenmapMinVector(GenmapVector x);
GenmapScalar GenmapNormVector(GenmapVector x, GenmapInt p);

int GenmapOrthogonalizebyOneVector(GenmapComm c, GenmapVector q1, GenmapLong n);

int GenmapPrintVector(GenmapVector x);
int GenmapDestroyVector(GenmapVector x);

/* Laplacian */
struct array *GenmapFindNeighbors(genmap_handle h,GenmapComm c);

int GenmapInitLaplacian(genmap_handle h, GenmapComm c);
int GenmapLaplacian(genmap_handle h, GenmapComm c, GenmapScalar *u, GenmapScalar *v);
int GenmapInitLaplacianWeighted(genmap_handle h, GenmapComm c);
int GenmapLaplacianWeighted(genmap_handle h, GenmapComm c, GenmapScalar *u, GenmapScalar *v);

/* Eigen */
int GenmapInvPowerIter(GenmapVector eVector, GenmapVector alpha, GenmapVector beta, GenmapVector init, int iter);
int GenmapTQLI(genmap_handle h, GenmapVector diag, GenmapVector upper, GenmapVector **eVec, GenmapVector *eVal);
int genmap_inverse_power(double *y, int N, double *A, int verbose);
int genmap_power(double *y, int N, double *A, int verbose);

/* Lanczos */
int GenmapLanczosLegendary(genmap_handle h, GenmapComm c, GenmapVector f,
                           GenmapInt niter, GenmapVector **rr, GenmapVector diag,
                           GenmapVector upper);
int GenmapLanczos(genmap_handle h, GenmapComm c, GenmapVector init,
                  GenmapInt iter, GenmapVector **q, GenmapVector alpha,
                  GenmapVector beta);

/* Fiedler */
int GenmapFiedlerLanczos(genmap_handle h,GenmapComm c,int maxIter,
  int global);
int GenmapFiedlerRQI(genmap_handle h,GenmapComm c,int maxIter,int global);

/* RSB */
void genmap_rsb(genmap_handle h,int verbose);

double GenmapGetMaxRss();
void GenmapPrintStack();

void matrix_inverse(int N,double *A);

#endif

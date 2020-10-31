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
typedef MPI_Comm GenmapCommExternal;

typedef struct GenmapComm_private *GenmapComm;
typedef struct GenmapHandle_private *GenmapHandle;
typedef struct GenmapVector_private *GenmapVector;
typedef struct GenmapElement_private *GenmapElements;

int GenmapInit(GenmapHandle *h, GenmapCommExternal ce);
int GenmapFinalize(GenmapHandle h);

int GenmapMallocArray(size_t n, size_t unit, void *p);
int GenmapCallocArray(size_t n, size_t unit, void *p);
int GenmapReallocArray(size_t n, size_t unit, void *p);
int GenmapFree(void *p);

GenmapElements GenmapGetElements(GenmapHandle h);
void GenmapSetElements(GenmapHandle h, GenmapElements elements);

GenmapComm GenmapGetLocalComm(GenmapHandle h);
void GenmapSetLocalComm(GenmapHandle h, GenmapComm c);

GenmapComm GenmapGetGlobalComm(GenmapHandle h);
void GenmapSetGlobalComm(GenmapHandle h, GenmapComm c);

GenmapInt GenmapGetNLocalElements(GenmapHandle h);
void GenmapSetNLocalElements(GenmapHandle h, GenmapInt localElements);

GenmapLong GenmapGetNGlobalElements(GenmapHandle h);
void GenmapSetNGlobalElements(GenmapHandle h, GenmapLong globalElements);

GenmapLong GenmapGetLocalStartIndex(GenmapHandle h);
void GenmapSetLocalStartIndex(GenmapHandle h, GenmapLong localStart);

int GenmapGetNVertices(GenmapHandle h);
void GenmapSetNVertices(GenmapHandle, int nVertices);

void GenmapScan(GenmapHandle h, GenmapComm c);

int GenmapCreateComm(GenmapComm *c, GenmapCommExternal ce);
int GenmapCommSize(GenmapComm c);
int GenmapCommRank(GenmapComm c);

int GenmapGop(GenmapComm c, void *v, GenmapInt size, GenmapDataType type,
              GenmapInt op);
int GenmapReduce(GenmapComm c, void *out, void *in, GenmapInt size,
                 GenmapDataType type, GenmapInt op);
int GenmapBcast(GenmapComm c, void *in, GenmapInt count, GenmapDataType type);

int GenmapDestroyComm(GenmapComm c);
void GenmapSplitComm(GenmapHandle h, GenmapComm *c, int bin);
int GenmapCrystalInit(GenmapHandle h, GenmapComm c);
int GenmapCrystalTransfer(GenmapHandle h, int field);
int GenmapCrystalFinalize(GenmapHandle h);

int GenmapRead(GenmapHandle h, void *data);

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

int GenmapPrintVector(GenmapVector x);
int GenmapDestroyVector(GenmapVector x);

/* Laplacian */
struct array *GenmapFindNeighbors(GenmapHandle h,GenmapComm c);

int GenmapInitLaplacianWeighted(GenmapHandle h, GenmapComm c,
                        GenmapVector weights);
int GenmapInitLaplacian(GenmapHandle h, GenmapComm c);

int GenmapLaplacianWeighted(GenmapHandle h, GenmapComm c, GenmapVector u,
                    GenmapVector weights, GenmapVector v);
int GenmapLaplacian(GenmapHandle h, GenmapComm c, GenmapVector u,
                    GenmapVector v);

/* Eigen */
int GenmapInvPowerIter(GenmapVector eVector, GenmapVector alpha,
                       GenmapVector beta, GenmapVector init, int iter);
int GenmapTQLI(GenmapHandle h, GenmapVector diagonal, GenmapVector upper,
               GenmapVector **eVectors, GenmapVector *eValues);

int GenmapOrthogonalizebyOneVector(GenmapHandle h, GenmapComm c,
                                   GenmapVector q1, GenmapLong n);

/* Lanczos */
int GenmapLanczosLegendary(GenmapHandle h, GenmapComm c, GenmapVector f,
                           GenmapInt niter, GenmapVector **rr, GenmapVector diag,
                           GenmapVector upper);
int GenmapLanczos(GenmapHandle h, GenmapComm c, GenmapVector init,
                  GenmapInt iter, GenmapVector **q, GenmapVector alpha,
                  GenmapVector beta);

/* Fiedler */
int GenmapFiedlerLanczos(GenmapHandle h,GenmapComm c,int maxIter,
  int global);
int GenmapFiedlerRQI(GenmapHandle h,GenmapComm c,int maxIter,int global);

/* RSB */
void GenmapRSB(GenmapHandle h,int verbose);

double GenmapGetMaxRss();
void GenmapPrintStack();
#endif

#ifndef _PARRSB_H_
#define _PARRSB_H_

#include <gslib.h>

typedef struct {
  /* General options */
  int global_partitioner; // -1 - None, 0 - RSB, 1 - RCB, 2 - RIB (Default: 0)
  int local_partitioner;  // -1 - None, 0 - RSB, 1 - RCB, 2 - RIB (Default: -1)
  int debug_level;        // 0, 1, 2, .. etc (Default: 0)
  int print_timing_info;  // 0 or 1 (Default: 0)

  /* RSB specific */
  int rsb_algo;         // 0 - Lanczos, 1 - MG (Default: 0)
  int rsb_prepartition; // 0 - None, 1 - RCB , 2 - RIB (Default: 1)
  int rsb_grammian;     // 0 or 1 (Default: 1)
  int rsb_paul;         // 0 or 1 (Default: 1)
} parRSB_options;

extern parRSB_options parrsb_default_options;

#define fparRSB_partMesh FORTRAN_UNPREFIXED(fparrsb_partmesh, FPARRSB_PARTMESH)
void fparRSB_partMesh(int *part, int *seq, long long *vtx, double *coord,
                      int *nel, int *nve, int *options, int *comm, int *err);

int parRSB_partMesh(int *part, int *seq, long long *vtx, double *coord, int nel,
                    int nv, parRSB_options *options, MPI_Comm comm);

#define fparRSB_findConnectivity                                               \
  FORTRAN_UNPREFIXED(fparrsb_findconnectivity, FPARRSB_FINDCONNECTIVITY)
void fparRSB_findConnectivity(long long *vertexId, double *coord, int *nel,
                              int *nDim, long long *periodicInfo,
                              int *nPeriodicFaces, double *tol, MPI_Fint *fcomm,
                              int *verbose, int *err);

int parRSB_findConnectivity(long long *vertexid, double *coord, int nel,
                            int nDim, long long *periodicInfo,
                            int nPeriodicFaces, double tol, MPI_Comm comm,
                            int verbose);

#endif

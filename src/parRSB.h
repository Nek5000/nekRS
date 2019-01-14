// API

#include "genmap-gslib.h"

#define fparRSB_partMesh FORTRAN_UNPREFIXED(fparrsb_partmesh,FPARRSB_PARTMESH)
void fparRSB_partMesh(long long *egl, long long *vl, int *negl,
                      long long *eglcon, long long *vlcon, int *neglcon,
                      int *nve, int *comm, int *err);

int parRSB_partMesh(long long *egl, long long *vl, int *negl,
                    long long *eglcon, long long *vlcon, int neglcon,
                    int nve, MPI_Comm comm);

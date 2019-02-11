// API
#include "genmap-gslib.h"

#define fparRSB_partMesh FORTRAN_UNPREFIXED(fparrsb_partmesh,FPARRSB_PARTMESH)
void fparRSB_partMesh(int *part, long long *vtx, int *nel,
                      int *nve, int *options, int *comm, int *err);

int parRSB_partMesh(int *part, long long *vtx, int nel, int nve,
                     int *options, MPI_Comm comm);

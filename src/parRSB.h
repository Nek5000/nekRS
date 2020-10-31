#ifndef _PARRSB_H_
#define _PARRSB_H_

#define fparRSB_partMesh \
  FORTRAN_UNPREFIXED(fparrsb_partmesh,FPARRSB_PARTMESH)
void fparRSB_partMesh(int *part,long long *vtx,int *nel,
  int *nve,int *options,int *comm,int *err);

int parRSB_partMesh(int *part,long long *vtx,int nel,int nve,
  int *options,MPI_Comm comm);

#define fparRCB_partMesh \
  FORTRAN_UNPREFIXED(fparrcb_partmesh,FPARRCB_PARTMESH)
void fparRCB_partMesh(int *part,int *seq,double *vtx,int *nel,int *nv,
  int *options,int *comm,int *err);

int parRCB_partMesh  (int *part,int *seq,double *vtx,int  nel,int  nv,
  int *options,MPI_Comm comm);

#define fparRSB_findConnectivity\
  FORTRAN_UNPREFIXED(fparrsb_findconnectivity,FPARRSB_FINDCONNECTIVITY)
void fparRSB_findConnectivity(long long *vertexId,double *coord,int *nel,
  int *nDim,long long *periodicInfo,int *nPeriodicFaces,double *tol,
  MPI_Fint *fcomm,int *verbose,int *err);

int parRSB_findConnectivity(long long *vertexid,double *coord,int nel,
  int nDim,long long *periodicInfo,int nPeriodicFaces,double tol,
  MPI_Comm comm,int verbose);
#endif

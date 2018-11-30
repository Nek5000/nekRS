// API

#define fparRSB_partMesh FORTRAN_NAME(fparrsb_partmesh,FPARRSB_PARTMESH)
void fparRSB_partMesh(long long *egl   , long long *vl   , int *negl,
                      long long *eglcon, long long *vlcon, int *neglcon,
                      int *nve, int *comm, int *err);

int parRSB_partMesh(GenmapLong *egl, GenmapLong *vl, GenmapInt *negl,
                    GenmapLong *eglcon, GenmapLong *vlcon, GenmapInt *neglcon,
                    GenmapInt  *nve, GenmapCommExternal comm);

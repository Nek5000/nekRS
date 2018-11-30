#ifndef _GENMAP_GMSH_H_
# define _GENMAP_GMSH_H_

#include <genmap-impl.h>

#include <stddef.h>
#include <stdlib.h>
//
// GenmapHandle
//
int GenmapCreateHandle_gmsh(GenmapHandle h);
//
// File I/O
//
int GenmapRead_gmsh(GenmapHandle h, void *data);
int GenmapWrite_gmsh(GenmapHandle h, char *fileNameBase);
int GenmapWriteGmsh(GenmapScalar *VX, GenmapScalar *VY,
                    GenmapScalar *VZ,
                    GenmapInt *EToV, GenmapLong Nnodes, GenmapLong Nelements,
                    char *fileName);
#endif

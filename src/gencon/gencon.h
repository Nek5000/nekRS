#ifndef _GENMAP_GENCON_H_
#define _GENMAP_GENCON_H_

#include <genmap-gslib.h>
#include <genmap-types.h>

/* Upper bound for number of dimensions */
#define GC_MAX_DIM 3

/* Boundary condition types */
#define GC_PERIODIC "P  "

/* Upper bounds for elements */
#define GC_MAX_FACES 6
#define GC_MAX_VERTICES 8
#define GC_MAX_NEIGHBORS 3

/* Upper bounds for faces */
#define GC_MAX_FACE_VERTICES 4

/* Header lengths */
#define GC_RE2_HEADER_LEN 80
#define GC_CO2_HEADER_LEN 132

typedef struct Point_private *Point;
typedef struct Element_private *Element;
typedef struct Boundary_private *BoundaryFace;
typedef struct Mesh_private *Mesh;

extern int NEIGHBOR_MAP[GC_MAX_VERTICES][GC_MAX_NEIGHBORS];
extern int PRE_TO_SYM_VERTEX[GC_MAX_VERTICES];
extern int PRE_TO_SYM_FACE[GC_MAX_FACES];

/* Mesh */
int mesh_init(Mesh *m, int nel, int nDim);
Element MeshGetElements(Mesh m);
void get_vertex_ids(long long **ids, Mesh m);
void get_vertex_coordinates(double **coords, Mesh m);
int get_mesh_dim(Mesh m);
int get_mesh_nel(Mesh m);
int mesh_free(Mesh m);

/* Read Nek5000 mesh files */
int read_geometry(Mesh *mesh, char *fname, struct comm *c);
int read_connectivity(Mesh mesh, char *fname, struct comm *c);
int write_connectivity(Mesh mesh, char *fname, struct comm *c);
int read_co2_mesh(Mesh *mesh, char *fname, struct comm *c);

/* Connectivity */
int findMinNeighborDistance(Mesh mesh);
int findSegments(Mesh mesh, struct comm *c, GenmapScalar tol, int verbose, buffer *bfr);
int faceCheck(Mesh mesh, struct comm *c);
int setGlobalID(Mesh mesh, struct comm *c);
int sendBack(Mesh mesh, struct comm *c, buffer *bfr);
int matchPeriodicFaces(Mesh mesh, struct comm *c, buffer *bfr);

#endif

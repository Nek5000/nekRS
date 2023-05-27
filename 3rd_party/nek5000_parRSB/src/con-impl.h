#if !defined(_CON_IMPL_H_)
#define _CON_IMPL_H_

#include "parrsb-impl.h"
#include "sort.h"
#include <stdarg.h>

/*
 Preprocessor Corner notation:      Symmetric Corner notation:

         4+-----+3    ^ s                    3+-----+4    ^ s
         /     /|     |                      /     /|     |
        /     / |     |                     /     / |     |
      8+-----+7 +2    +----> r            7+-----+8 +2    +----> r
       |     | /     /                     |     | /     /
       |     |/     /                      |     |/     /
      5+-----+6    t                      5+-----+6    t



                   i) Preprocessor notation:

                                     +--------+     ^ S
                                    /        /|     |
                                   /    3   / |     |
                             4--> /        /  |     |
                                 +--------+ 2 +     +----> R
                                 |        |  /     /
                                 |    6   | /     /
                                 |        |/     /
                                 +--------+     T
                                     1

                  ii) Symmetric notation:

                                     +--------+     ^ S
                                    /        /|     |
                                   /    4   / |     |
                             1--> /        /  |     |
                                 +--------+ 2 +     +----> R
                                 |        |  /     /
                                 |    6   | /     /
                                 |        |/     /
                                 +--------+     T
                                     3

   EFACE(IFACE)  - Given face number IFACE in symmetric notation,
                   returns preprocessor notation face number.

   EFACE1(IFACE) - Given face number IFACE in preprocessor notation,
                   returns symmetric notation face number.

The following variables all take the symmetric notation of IFACE
as arguments:

   ICFACE(i,IFACE) - Gives the 4 vertices which reside on face IFACE
                     as depicted below, e.g. ICFACE(i,2)=2,4,6,8.

                      3+-----+4    ^ Y
                      /  2  /|     |
   Edge 1 extends    /     / |     |
     from vertex   7+-----+8 +2    +----> X
     1 to 2.        |  4  | /     /
                    |     |/     /
                   5+-----+6    Z
                       3

   IEDGFC(i,IFACE) - Gives the 4 edges which border the face IFACE
                     Edge numbering is as follows:
                        Edge = 1,2,3,4     run in +r direction
                        Edge = 5,6,7,8     run in +s direction
                        Edge = 9,10,11,12  run in +t direction

                     Ordering of each edge is such that a monotonically
                     increasing sequence of vertices is associated with
                     the start point of a corresponding set of
                     monotonically increasing edge numbers, e.g.,

   ICEDG(i,IEDGE)  - Gives 3 variables for determining the stride along
                     a given edge, IEDGE;  i=1 gives the starting vertex
                                           i=2 gives the stopping vertex
                                           i=3 gives the stride size.

*/

// Macros for min and max.
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

// Upper bounds for elements and face quantities.
#define GC_MAX_FACES 6
#define GC_MAX_VERTICES 8
#define GC_MAX_NEIGHBORS 3
#define GC_MAX_FACE_VERTICES 4

extern int faces3D[GC_MAX_FACES][GC_MAX_FACE_VERTICES];
extern int faces2D[GC_MAX_FACES][GC_MAX_FACE_VERTICES];

struct point_t {
  scalar dx, x[3];
  uint proc, origin;
  int ifSegment;
  ulong sequenceId, elementId, globalId, pntid;
};
typedef struct point_t *Point;

struct boundary_t {
  ulong elementId, faceId;
  struct point_t face[4];
  uint proc;
  long bc[2];
};
typedef struct boundary_t *BoundaryFace;

struct mesh_t {
  ulong nelgt;
  uint nelt, ndim, nv, nnbrs;
  struct array elements, boundary;
};
typedef struct mesh_t *Mesh;

double diff_sqr(double x, double y);

int send_back(Mesh mesh, struct comm *c, buffer *bfr);

int find_unique_vertices(Mesh mesh, struct comm *c, scalar tol, int verbose,
                         buffer *bfr);

int matchPeriodicFaces(Mesh mesh, struct comm *c, buffer *bfr);

int elementCheck(Mesh mesh, struct comm *c, buffer *bfr);

int faceCheck(Mesh mesh, struct comm *c, buffer *bfr);

#endif // _CON_IMPL_H_

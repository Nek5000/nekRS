#ifndef _GENMAP_GENCON_IMPL_H_
#define _GENMAP_GENCON_IMPL_H_

#include <gencon.h>

#define sqrDiff(x,y) (((x)-(y))*((x)-(y)))
#define distance2D(a,b)\
  (sqrDiff((a).x[0],(b).x[0])+sqrDiff((a).x[1],(b).x[1]))
#define distance3D(a,b) (distance2D(a,b)+sqrDiff((a).x[2],(b).x[2]))

#define min(a,b) ((a)<(b) ? (a) : (b))
#define max(a,b) ((a)>(b) ? (a) : (b))

#define READ_T(coords,buf,T,nVertex) do{\
  memcpy((coords),buf,sizeof(T)*nVertex);\
} while(0)

#define WRITE_INT(dest,val) do{\
  memcpy(dest,&(val),sizeof(int));\
} while(0)

// TODO: Use rsb_element here
struct Point_private{
  GenmapScalar dx;
  GenmapScalar x[GC_MAX_DIM];
  uint proc;
  uint origin;
  int ifSegment;
  ulong sequenceId;
  ulong elementId;
  slong globalId;
};

extern int faces3D[GC_MAX_FACES][GC_MAX_FACE_VERTICES];
extern int faces2D[GC_MAX_FACES][GC_MAX_FACE_VERTICES];

struct Face_private{
  struct Point_private vertex[GC_MAX_FACE_VERTICES];
};

struct Element_private{
  struct Point_private vertex[GC_MAX_VERTICES];
};

struct Boundary_private{
  ulong elementId;
  ulong faceId;
  struct Face_private face;
  uint proc;
  long bc[2];
};

struct Mesh_private{
  ulong nelgt;
  ulong nelgv;
  uint nelt;
  int nDim;
  int nVertex;
  int nNeighbors;
  struct array elements;
  struct array boundary;
};

int transferBoundaryFaces(Mesh mesh,struct comm *c);

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
#endif

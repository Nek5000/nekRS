#include "parrsb-impl.h"
#include "sort.h"

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

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

// Upper bounds for elements and face quantities
#define GC_MAX_FACES 6
#define GC_MAX_VERTICES 8
#define GC_MAX_NEIGHBORS 3
#define GC_MAX_FACE_VERTICES 4

struct point_t {
  scalar dx, x[3];
  uint proc, origin;
  int ifSegment;
  ulong sequenceId, elementId, globalId;
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
  ulong nelgt, nelgv;
  uint nelt;
  int nDim, nVertex, nNeighbors;
  struct array elements;
  struct array boundary;
};
typedef struct mesh_t *Mesh;

typedef struct {
  ulong sequenceId;
  int nNeighbors;
  ulong elementId;
  ulong vertexId;
  uint workProc;
} vertex;

//==============================================================================
// Mesh struct
//
int mesh_init(Mesh *m_, int nel, int ndim) {
  Mesh m = *m_ = tcalloc(struct mesh_t, 1);
  m->nelt = nel, m->nDim = ndim, m->nNeighbors = ndim;
  m->nVertex = (ndim == 2) ? 4 : 8;

  array_init(struct point_t, &m->elements, 10);
  array_init(struct boundary_t, &m->boundary, 10);

  return 0;
}

int mesh_free(Mesh m) {
  array_free(&m->elements);
  array_free(&m->boundary);
  free(m);

  return 0;
}

void get_vertex_ids(long long **vtx_, Mesh mesh) {
  int nelt = mesh->nelt;
  int nv = (mesh->nDim == 3) ? 8 : 4;

  long long *vtx = *vtx_ = tcalloc(long long, nelt *nv);

  struct point_t *ptr = (struct point_t *)mesh->elements.ptr;
  for (uint e = 0, count = 0; e < nelt; e++) {
    for (int v = 0; v < nv; v++)
      vtx[count] = ptr[count++].globalId;
  }
}

void get_vertex_coordinates(double **coords_, Mesh mesh) {
  int nelt = mesh->nelt;
  int ndim = mesh->nDim;
  int nv = (ndim == 3) ? 8 : 4;

  size_t size = nelt;
  size = size * nv * ndim;

  double *coords = *coords_ = tcalloc(double, size);

  struct point_t *ptr = (struct point_t *)mesh->elements.ptr;
  int e, v, d;
  int count = 0;
  for (e = 0; e < nelt; e++)
    for (v = 0; v < nv; v++)
      for (d = 0; d < ndim; d++) {
        coords[count] = ptr[e * nv + v].x[d];
        count++;
      }
}

int get_bcs(unsigned int *nbcs_, long long **bcs_, Mesh m) {
  unsigned int nbcs = *nbcs_ = m->boundary.n;
  long long *bcs = *bcs_ = tcalloc(long long, 4 * nbcs);

  struct boundary_t *ptr = m->boundary.ptr;
  uint i;
  for (i = 0; i < nbcs; i++) {
    bcs[4 * i + 0] = ptr[i].elementId;
    bcs[4 * i + 1] = ptr[i].faceId;
    bcs[4 * i + 2] = ptr[i].bc[0];
    bcs[4 * i + 3] = ptr[i].bc[1];
  }

  return 0;
}

int get_mesh_dim(Mesh mesh) { return mesh->nDim; }

int get_mesh_nel(Mesh mesh) { return mesh->nelt; }

//==============================================================================
// Find the minimum distance between a vertex and its neighbors
//
int PRE_TO_SYM_VERTEX[GC_MAX_VERTICES] = {0, 1, 3, 2, 4, 5, 7, 6};
int PRE_TO_SYM_FACE[GC_MAX_FACES] = {2, 1, 3, 0, 4, 5};
int NEIGHBOR_MAP[GC_MAX_VERTICES][GC_MAX_NEIGHBORS] = {
    {1, 2, 4}, {0, 3, 5}, {0, 3, 6}, {1, 2, 7},
    {0, 5, 6}, {1, 4, 7}, {2, 4, 7}, {3, 5, 6}};

static inline double diff_sqr(double x, double y) { return (x - y) * (x - y); }

static inline double distance_2d(struct point_t *a, struct point_t *b) {
  return diff_sqr(a->x[0], b->x[0]) + diff_sqr(a->x[1], b->x[1]);
}

static inline double distance_3d(struct point_t *a, struct point_t *b) {
  return distance_2d(a, b) + diff_sqr(a->x[2], b->x[2]);
}

int findMinNeighborDistance(Mesh mesh) {
  struct point_t *p = mesh->elements.ptr;
  int nDim = mesh->nDim;
  int nVertex = mesh->nVertex;

  uint i, j, k;
  int neighbor;
  scalar d;

  if (nDim == 3) {
    for (i = 0; i < mesh->elements.n; i += nVertex) {
      for (j = 0; j < nVertex; j++) {
        p[i + j].dx = SCALAR_MAX;
        for (k = 0; k < mesh->nNeighbors; k++) {
          neighbor = NEIGHBOR_MAP[j][k];
          d = distance_3d(&p[i + j], &p[i + neighbor]);
          p[i + j].dx = MIN(p[i + j].dx, d);
        }
      }
    }
  } else if (nDim == 2) {
    for (i = 0; i < mesh->elements.n; i += nVertex) {
      for (j = 0; j < nVertex; j++) {
        p[i + j].dx = SCALAR_MAX;
        for (k = 0; k < mesh->nNeighbors; k++) {
          neighbor = NEIGHBOR_MAP[j][k];
          d = distance_2d(&p[i + j], &p[i + neighbor]);
          p[i + j].dx = MIN(p[i + j].dx, d);
        }
      }
    }
  } else {
    return 1;
  }

  return 0;
}

//==============================================================================
// Identify segements
//
// Tuple sort
static void tuple_sort_(void *ra, uint n, uint usize, uint offset) {
  sint i, ir, j, l;
  void *rra = calloc(1, usize);
  assert(rra != NULL);

#define cpy(rra_, i_, ra_, l_)                                                 \
  {                                                                            \
    char *dst = (char *)rra_ + (i_ - 1) * usize;                               \
    char *src = (char *)ra_ + (l_ - 1) * usize;                                \
    memcpy(dst, src, usize);                                                   \
  }

#define get(ra_, l_) (*((double *)((char *)ra_ + (l_ - 1) * usize + offset)))

  if (n < 2)
    return;
  l = n / 2 + 1;
  ir = n;

  for (;;) {
    if (l > 1) {
      --l;
      assert(l >= 1 && l <= n && "l");
      cpy(rra, 1, ra, l);
    } else {
      cpy(rra, 1, ra, ir);
      cpy(ra, ir, ra, 1);
      if (--ir == 1) {
        cpy(ra, 1, rra, 1);
        break;
      }
    }
    i = l;
    j = l + l;
    while (j <= ir) {
      if (j < ir && get(ra, j) < get(ra, j + 1))
        j++;
      assert(j >= 1 && j <= n && "j2");
      assert(i >= 1 && i <= n && "i");
      if (get(rra, 1) < get(ra, j)) {
        cpy(ra, i, ra, j);
        i = j;
        j = 2 * j;
      } else
        break;
    }
    assert(i >= 1 && i <= n && "i2");
    cpy(ra, i, rra, 1);
  }

#undef cpy
#undef get

  if (rra != NULL)
    free(rra);
}

#define tuple_sort(T, arr, n, index)                                           \
  tuple_sort_((void *)arr, n, sizeof(T), offsetof(T, index))

static void initSegment(Mesh mesh, struct comm *c) {
  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  uint i;
  for (i = 0; i < nPoints; i++) {
    points[i].ifSegment = 0;
    points[i].globalId = 0;
  }

  /* First rank with nPoints > 0 will have ifSegment = 1 */
  sint rank = c->id;
  if (nPoints == 0)
    rank = c->np;

  sint buf[2];
  comm_allreduce(c, gs_int, gs_min, &rank, 1, buf);

  if (c->id == rank)
    points[0].ifSegment = 1;
}

static int sortSegmentsLocal(Mesh mesh, int dim, buffer *bfr) {
  sint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  sint s = 0, e;
  while (s < nPoints) {
    for (e = s + 1; e < nPoints && points[e].ifSegment == 0; e++)
      ;

    if (s < nPoints - 1 && e - s > 1) {
      switch (dim) {
      case 0:
        tuple_sort(struct point_t, &points[s], e - s, x[0]);
        break;
      case 1:
        tuple_sort(struct point_t, &points[s], e - s, x[1]);
        break;
      case 2:
        tuple_sort(struct point_t, &points[s], e - s, x[2]);
        break;
      default:
        break;
      }
    }

    sint i, sum = 0;
    for (i = s; i < e; i++) {
      sum += points[i].ifSegment;
      points[i].ifSegment = 0;
    }

    if (sum > 0)
      points[s].ifSegment = 1;

    s = e;
  }

  return 0;
}

static int sortSegments(Mesh mesh, struct comm *c, int dim, buffer *bfr) {
  if (c->np > 1) {
    // Parallel sort -- we haven't localized the problem yet
    switch (dim) {
    case 0:
      parallel_sort(struct point_t, &mesh->elements, x[0], gs_scalar, bin_sort,
                    1, c, bfr);
      break;
    case 1:
      parallel_sort(struct point_t, &mesh->elements, x[1], gs_scalar, bin_sort,
                    1, c, bfr);
      break;
    case 2:
      parallel_sort(struct point_t, &mesh->elements, x[2], gs_scalar, bin_sort,
                    1, c, bfr);
      break;
    default:
      break;
    }

    initSegment(mesh, c);
  } else {
    // Local sort: Segments are local
    sortSegmentsLocal(mesh, dim, bfr);
  }

  return 0;
}

static int sendLastPoint(struct array *arr, Mesh mesh, struct comm *c) {
  Point pts = mesh->elements.ptr;
  sint npts = mesh->elements.n;

  struct point_t lastp = pts[npts - 1];
  lastp.proc = (c->id + 1) % c->np;

  array_init(struct point_t, arr, 1);
  array_cat(struct point_t, arr, &lastp, 1);

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct point_t, arr, proc, 1, &cr);
  crystal_free(&cr);

  return 0;
}

static int findSegments(Mesh mesh, struct comm *c, int i, scalar tolSquared) {
  Point pts = mesh->elements.ptr;
  sint npts = mesh->elements.n;
  int nDim = mesh->nDim;

  sint j;
  for (j = 1; j < npts; j++) {
    scalar d = diff_sqr(pts[j].x[i], pts[j - 1].x[i]);
    scalar dx = MIN(pts[j].dx, pts[j - 1].dx) * tolSquared;

    if (d > dx)
      pts[j].ifSegment = 1;
  }

  if (c->np > 1) {
    struct array arr;
    sendLastPoint(&arr, mesh, c);

    if (c->id > 0) {
      struct point_t *lastp = arr.ptr;
      scalar d = diff_sqr(lastp->x[i], pts->x[i]);
      scalar dx = MIN(lastp->dx, pts->dx) * tolSquared;
      if (d > dx)
        pts->ifSegment = 1;
    }

    array_free(&arr);
  }

  return 0;
}

static slong countSegments(Mesh mesh, struct comm *c) {
  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  sint count = 0, i;
  for (i = 0; i < nPoints; i++)
    if (points[i].ifSegment > 0)
      count++;

  slong buf[2][1];
  slong in = count;
  comm_allreduce(c, gs_long, gs_add, &in, 1, buf);

  return in;
}

static int setProc(Mesh mesh, sint rankg, uint index, int inc_proc,
                   struct comm *c) {
  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  slong size[2] = {0};
  if (c->id < rankg)
    size[0] = nPoints;
  if (c->id == rankg) {
    size[0] = index;
    size[1] = nPoints - index;
  }
  if (c->id > rankg)
    size[1] = nPoints;

  slong out[2][2], buf[2][2];
  comm_scan(out, c, gs_long, gs_add, size, 2, buf);

  sint np[2] = {0};
  if (c->id < rankg)
    np[0] = 1;
  if (c->id == rankg) {
    np[0] = inc_proc;
    np[1] = 1 - inc_proc;
  }
  if (c->id > rankg)
    np[1] = 1;

  comm_allreduce(c, gs_int, gs_add, np, 2, buf);

  sint low_size = (out[1][0] + np[0] - 1) / np[0];
  sint high_size = (out[1][1] + np[1] - 1) / np[1];

  uint i;
  for (i = 0; i < size[0]; i++) {
    points[i].globalId = out[0][0] + i;
    points[i].proc = (out[0][0] + i) / low_size;
  }

  for (i = size[0]; i < size[0] + size[1]; i++) {
    points[i].globalId = out[0][1] + i - size[0];
    points[i].proc = np[0] + (out[0][1] + i - size[0]) / high_size;
  }

  return 0;
}

static int rearrangeSegments(Mesh mesh, struct comm *seg, buffer *bfr) {
  while (seg->np > 1 && countSegments(mesh, seg) > 1) {
    uint nPoints = mesh->elements.n;
    Point points = mesh->elements.ptr;

    /* comm_scan */
    slong out[2][1], buf[2][1];
    slong in = nPoints;
    comm_scan(out, seg, gs_long, gs_add, &in, 1, buf);
    slong start = out[0][0];
    slong nelg = out[1][0];

    double min = DBL_MAX;
    int inc_proc = 1;
    uint index = (seg->id == 0) ? 1 : 0;

    uint i;
    for (i = index; i < nPoints; i++) {
      if (points[i].ifSegment > 0) {
        double f0 = fabs((start + i + 0.0) / nelg - (seg->id + 0.0) / seg->np);
        if (seg->id == 0)
          f0 = DBL_MAX;
        double f1 = fabs((start + i + 0.0) / nelg - (seg->id + 1.0) / seg->np);
        if (seg->id == seg->np - 1)
          f1 = DBL_MAX;

        if (f0 < min) {
          inc_proc = 0;
          min = f0;
          index = i;
        }
        if (f1 < min) {
          inc_proc = 1;
          min = f1;
          index = i;
        }
      }
    }

    double dbuf[2];
    double ming = min;
    comm_allreduce(seg, gs_double, gs_min, &ming, 1, dbuf);

    sint rankg = -1;
    if (fabs(ming - min) < 1e-15)
      rankg = seg->id;
    comm_allreduce(seg, gs_int, gs_max, &rankg, 1, buf);

    setProc(mesh, rankg, index, inc_proc, seg);

    int bin = 1;
    if (seg->id < rankg)
      bin = 0;
    if (seg->id == rankg && inc_proc == 1)
      bin = 0;

    struct crystal cr;
    crystal_init(&cr, seg);
    sarray_transfer(struct point_t, &mesh->elements, proc, 0, &cr);
    crystal_free(&cr);

    struct comm new;
    comm_split(seg, bin, seg->id, &new);
    comm_free(seg);
    comm_dup(seg, &new);
    comm_free(&new);

    parallel_sort(struct point_t, &mesh->elements, globalId, gs_long, bin_sort,
                  1, seg, bfr);
  }

  return 0;
}

static int findUniqueVertices(Mesh mesh, struct comm *c, scalar tol,
                              int verbose, buffer *bfr) {
  scalar tolSquared = tol * tol;
  int nDim = mesh->nDim;

  initSegment(mesh, c);

  struct comm seg;
  comm_dup(&seg, c);

  int t, d;
  for (t = 0; t < nDim; t++) {
    for (d = 0; d < nDim; d++) {
      sortSegments(mesh, &seg, d, bfr);
      findSegments(mesh, &seg, d, tolSquared);

      slong n_pts = mesh->elements.n;
      slong buf[2];
      comm_allreduce(c, gs_long, gs_add, &n_pts, 1, buf);

      slong n_seg = countSegments(mesh, c);
      if (c->id == 0 && verbose > 1)
        printf("locglob: %d %d %lld %lld\n", t + 1, d + 1, n_seg, n_pts);

      rearrangeSegments(mesh, &seg, bfr);
      parrsb_barrier(c);
    }
  }

  comm_free(&seg);

  return 0;
}

#undef tuple_sort

//==============================================================================
// Global numbering
//
static int setGlobalID(Mesh mesh, struct comm *c) {
  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  sint bin = 1;
  if (nPoints == 0)
    bin = 0;

  comm_ext old = c->c;
  struct comm nonZeroRanks;
#ifdef MPI
  MPI_Comm new;
  MPI_Comm_split(old, bin, c->id, &new);
  comm_init(&nonZeroRanks, new);
  MPI_Comm_free(&new);
#else
  comm_init(&nonZeroRanks, 1);
#endif

  sint rank = nonZeroRanks.id;
  sint size = nonZeroRanks.np;

  if (bin == 1) {
    slong count = 0;
    sint i;
    for (i = 0; i < nPoints; i++)
      if (points[i].ifSegment)
        count++;

    slong out[2][1], buf[2][1];
    slong in = count;
    comm_scan(out, &nonZeroRanks, gs_long, gs_add, &in, 1, buf);
    slong start = out[0][0];

    count = -1;
    for (i = 0; i < nPoints; i++) {
      if (points[i].ifSegment)
        count++;
      points[i].globalId = start + count;
    }
  }

  comm_free(&nonZeroRanks);

  return 0;
}

static int sendBack(Mesh mesh, struct comm *c, buffer *bfr) {
  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct point_t, &mesh->elements, origin, 0, &cr);
  crystal_free(&cr);

  sarray_sort(struct point_t, mesh->elements.ptr, mesh->elements.n, sequenceId,
              1, bfr);

  return 0;
}

//==============================================================================
// Handle periodic BCs
//
int faces3D[GC_MAX_FACES][GC_MAX_FACE_VERTICES] = {{1, 5, 7, 3}, {2, 4, 8, 6},
                                                   {1, 2, 6, 5}, {3, 7, 8, 4},
                                                   {1, 3, 4, 2}, {5, 6, 8, 7}};

int faces2D[GC_MAX_FACES][GC_MAX_FACE_VERTICES] = {{3, 1, 0, 0}, {2, 4, 0, 0},
                                                   {1, 2, 0, 0}, {4, 3, 0, 0},
                                                   {0, 0, 0, 0}, {0, 0, 0, 0}};

#define distance2D(a, b) (diff_sqr(a.x[0], b.x[0]) + diff_sqr(a.x[1], b.x[1]))
#define distance3D(a, b) (distance2D(a, b) + diff_sqr(a.x[2], b.x[2]))

struct mpair_t {
  uint proc;
  ulong orig, min;
};

static int compressPeriodicVertices(Mesh mesh, struct comm *c, buffer *bfr) {
  parallel_sort(struct point_t, &mesh->elements, globalId, gs_long, 0, 0, c,
                bfr);

  Point points = mesh->elements.ptr;
  uint npoints = mesh->elements.n;

  sint i, nunique = 0;
  if (npoints > 0) {
    slong current = points[0].globalId;
    points[0].globalId = nunique;
    for (i = 1; i < npoints; i++)
      if (points[i].globalId == current)
        points[i].globalId = nunique;
      else {
        current = points[i].globalId, ++nunique;
        points[i].globalId = nunique;
      }
  }

  slong out[2][1], buf[2][1], in[1];
  if (npoints > 0)
    in[0] = nunique + 1;
  else
    in[0] = 0;
  comm_scan(out, c, gs_long, gs_add, in, 1, buf);
  slong start = out[0][0];

  for (i = 0; i < npoints; i++)
    points[i].globalId += start;

  return 0;
}

static ulong findMinBelowI(ulong min, uint I, struct array *arr) {
  struct mpair_t *ptr = (struct mpair_t *)arr->ptr;
  for (uint i = 0; i < I; i++)
    if (ptr[i].orig == min)
      return ptr[i].min;
  return min;
}

static int renumberPeriodicVertices(Mesh mesh, struct comm *c,
                                    struct array *matched, buffer *bfr) {
  uint size1 = mesh->elements.n, size2 = matched->n;
  slong *mids = tcalloc(slong, size1 + 2 * size2),
        *mnew = tcalloc(slong, size1 + 2 * size2),
        *mcur = tcalloc(slong, size1);

  struct point_t *pe = (struct point_t *)mesh->elements.ptr;
  for (uint i = 0; i < size1; i++)
    mids[i] = pe[i].globalId;
  struct mpair_t *pm = (struct mpair_t *)matched->ptr;
  for (uint i = 0; i < size2; i++)
    mids[size1 + i] = pm[i].orig;
  struct gs_data *gsh = gs_setup(mids, size1 + size2, c, 0, gs_pairwise, 0);

  for (uint i = 0; i < size1; i++)
    mnew[i] = pe[i].globalId;
  for (uint i = 0; i < size2; i++)
    mnew[size1 + i] = pm[i].min;
  gs(mnew, gs_long, gs_min, 0, gsh, bfr);

  sint changed, wrk;
  do {
    for (uint i = 0; i < size1; i++)
      mcur[i] = mnew[i];
    for (uint i = 0; i < size2; i++)
      mids[size1 + size2 + i] = -mnew[size1 + i];
    struct gs_data *gsh1 =
        gs_setup(mids, size1 + 2 * size2, c, 0, gs_pairwise, 0);

    gs(mnew, gs_long, gs_min, 0, gsh1, bfr);
    gs_free(gsh1);

    for (uint i = 0; i < size2; i++)
      mnew[size1 + i] = mnew[size1 + size2 + i];
    gs(mnew, gs_long, gs_min, 0, gsh, bfr);

    changed = 0;
    for (uint i = 0; i < size1; i++)
      changed += (mnew[i] != mcur[i]);
    comm_allreduce(c, gs_int, gs_max, &changed, 1, &wrk);
  } while (changed);

  for (uint i = 0; i < size1; i++)
    pe[i].globalId = mcur[i];

  gs_free(gsh);
  free(mids), free(mnew), free(mcur);

  return 0;
}

static int findConnectedPeriodicPairs(Mesh mesh, BoundaryFace f_,
                                      BoundaryFace g_, struct array *matched) {
  struct boundary_t f = *f_, g = *g_;

  int nvf = mesh->nVertex / 2;
  int nDim = mesh->nDim;

  int i, j;
  scalar fMax = 0.0, gMax = 0.0;

  for (i = 0; i < nDim; i++) {
    scalar meanF = 0.0, meanG = 0.0;

    for (j = 0; j < nvf; j++) {
      fMax = MAX(fMax, fabs(f.face[j].x[i]));
      gMax = MAX(gMax, fabs(g.face[j].x[i]));
      meanF += f.face[j].x[i];
      meanG += g.face[j].x[i];
    }

    for (j = 0; j < nvf; j++) {
      f.face[j].x[i] -= (meanF / nvf);
      g.face[j].x[i] -= (meanG / nvf);
    }
  }

  int shift = 0, k;
  scalar d2Min = 1.e20, d2;
  for (i = 0; i < nvf; i++) {
    d2 = 0.0;
    for (j = 0; j < nvf; j++) {
      k = (j + i) % nvf;
      k = nvf - 1 - k;
      if (nDim == 3)
        d2 += distance3D(f.face[j], g.face[k]);
      else if (nDim == 2)
        d2 += distance2D(f.face[j], g.face[k]);
    }
    if (d2 < d2Min) {
      d2Min = d2;
      shift = i;
    }
  }
  d2Min = sqrt(d2Min);

  scalar fgMax = MAX(fMax, gMax);
  scalar tol = (1e-3) * fgMax;
  if (d2Min > tol) {
    fprintf(stderr,
            "Faces did not match: (d2Min,tol,face1,face2): "
            "%lf %lf %lld %lld\n",
            d2Min, tol, f.faceId, g.faceId);
    exit(1);
  }

  struct mpair_t m;
  for (i = 0; i < nvf; i++) {
    k = (i + shift) % nvf;
    k = nvf - 1 - k;
    m.min = MIN(f.face[i].globalId, g.face[k].globalId);
    m.orig = MAX(f.face[i].globalId, g.face[k].globalId);
    array_cat(struct mpair_t, matched, &m, 1);
  }
}

static int findConnectedPeriodicFaces(Mesh mesh, struct array *matched) {
  sint bSize = mesh->boundary.n;
  BoundaryFace ptr = mesh->boundary.ptr;
  sint i, j;

  for (i = 0; i < bSize - 1; i++) {
    for (j = i + 1; j < bSize; j++)
      if (ptr[j].bc[0] == ptr[i].elementId && ptr[j].bc[1] == ptr[i].faceId) {
        findConnectedPeriodicPairs(mesh, &ptr[i], &ptr[j], matched);
      }
  }
}

static int gatherMatchingPeriodicFaces(Mesh mesh, struct comm *c) {
  int size = c->np, rank = c->id;

  BoundaryFace bPtr = mesh->boundary.ptr;
  int nFaces = mesh->boundary.n;

  slong nelgt = mesh->nelgt;
  sint nelt = nelgt / size;
  sint nrem = nelgt - nelt * size;
  slong N = (size - nrem) * nelt;

  sint i;
  slong eid;
  for (i = 0; i < nFaces; i++) {
    eid = MAX(bPtr[i].bc[0], bPtr[i].elementId);
    if (eid < N)
      bPtr[i].proc = eid / nelt;
    else
      bPtr[i].proc = (eid - N) / (nelt + 1) + size - nrem;
  }

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct boundary_t, &mesh->boundary, proc, 1, &cr);
  crystal_free(&cr);
}

static int setPeriodicFaceCoordinates(Mesh mesh, struct comm *c, buffer *buf) {
  BoundaryFace bPtr = mesh->boundary.ptr;
  sint bSize = mesh->boundary.n;
  if (bSize == 0)
    return 0;

  Point ePtr = mesh->elements.ptr;
  sint eSize = mesh->elements.n;
  if (eSize == 0)
    return 0;

  /* Need boundary array to be sorted by elementId */
  sarray_sort(struct boundary_t, bPtr, bSize, elementId, 1, buf);

  /* Need element array to be sorted by sequenceId */
  sarray_sort(struct point_t, ePtr, eSize, sequenceId, 1, buf);

  int faces[GC_MAX_FACES][GC_MAX_FACE_VERTICES];
  if (mesh->nDim == 3)
    memcpy(faces, faces3D, GC_MAX_FACES * GC_MAX_FACE_VERTICES * sizeof(int));
  else
    memcpy(faces, faces2D, GC_MAX_FACES * GC_MAX_FACE_VERTICES * sizeof(int));

  sint i = 0, k = 0;
  int nv = mesh->nVertex, nvf = mesh->nVertex / 2, j;
  while (i < bSize) {
    while (k < eSize && ePtr[k].elementId < bPtr[i].elementId)
      k += nv;
    // copy vertices to boundary face
    if (k < eSize && ePtr[k].elementId == bPtr[i].elementId) {
      int faceId = bPtr[i].faceId;
      for (j = 0; j < nvf; j++)
        bPtr[i].face[j] = ePtr[k + faces[faceId][j] - 1];
    }
    i++;
  }
}

static int matchPeriodicFaces(Mesh mesh, struct comm *c, buffer *bfr) {
  setPeriodicFaceCoordinates(mesh, c, bfr);
  gatherMatchingPeriodicFaces(mesh, c);

  struct array matched;
  array_init(struct mpair_t, &matched, 10);
  matched.n = 0;

  findConnectedPeriodicFaces(mesh, &matched);
  renumberPeriodicVertices(mesh, c, &matched, bfr);
  array_free(&matched);

  compressPeriodicVertices(mesh, c, bfr);
  sendBack(mesh, c, bfr);

  return 0;
}

#undef distance2D
#undef distance3D

//==============================================================================
// Various checks
//
typedef struct {
  ulong *elements;
  uint *offsets;
  ulong *globalIds;
  uint size;
} VToEMap;

typedef struct {
  ulong id;
} LongID;

typedef struct {
  ulong procId;
} ProcID;

static VToEMap *getVToEMap(Mesh m, struct comm *c, buffer *bfr) {
  sint nelt = m->nelt;
  sint nv = m->nVertex;

  slong out[2][1], buf[2][1], in = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, buf);
  ulong elemId = out[0][0];
  ulong sequenceId = elemId * nv;

  size_t size = nelt * (size_t)nv;
  struct array vertices;
  array_init(vertex, &vertices, size);

  // Create (globalId, elementId) pairs and send them to globalId % np
  Point ptr = m->elements.ptr;
  sint i, j;
  for (i = 0; i < nelt; i++) {
    for (j = 0; j < nv; j++) {
      ulong globalId = ptr[i * nv + j].globalId + 1;
      vertex t = {.elementId = elemId,
                  .sequenceId = sequenceId,
                  .vertexId = globalId,
                  .workProc = globalId % c->np};
      array_cat(vertex, &vertices, &t, 1);
      sequenceId++;
    }

    elemId++;
  }

  sarray_sort_2(vertex, vertices.ptr, vertices.n, vertexId, 1, elementId, 1,
                bfr);

  struct array vtcsCmpct;
  array_init(vertex, &vtcsCmpct, 10);
  vertex *vPtr = (vertex *)vertices.ptr;

  if (vertices.n > 0) {
    vertex prev = vPtr[0];
    array_cat(vertex, &vtcsCmpct, &prev, 1);

    for (i = 1; i < vertices.n; i++) {
      if ((vPtr[i].elementId != prev.elementId) ||
          (vPtr[i].vertexId != prev.vertexId)) {
        prev = vPtr[i];
        array_cat(vertex, &vtcsCmpct, &prev, 1);
      }
    }
  }
  array_free(&vertices);

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(vertex, &vtcsCmpct, workProc, 1, &cr);

  // Find all the elements which share globalId and send the union
  // back to all the processors which has globalId
  // FIXME: Assumes quads or hexes
  vPtr = vtcsCmpct.ptr;
  sarray_sort_2(vertex, vPtr, vtcsCmpct.n, vertexId, 1, workProc, 0, bfr);

  struct array a;
  array_init(vertex, &a, 10);
  struct array procs;
  array_init(ProcID, &procs, 10);

  vPtr = vtcsCmpct.ptr;
  sint s = 0, e;
  vertex t;
  ProcID p;
  while (s < vtcsCmpct.n) {
    procs.n = 0;

    p.procId = vPtr[s].workProc;
    array_cat(ProcID, &procs, &p, 1);
    for (e = s + 1; e < vtcsCmpct.n && vPtr[s].vertexId == vPtr[e].vertexId;
         e++) {
      if (vPtr[e].workProc != p.procId) {
        p.procId = vPtr[e].workProc;
        array_cat(ProcID, &procs, &p, 1);
      }
    }

    ProcID *pPtr = procs.ptr;
    e = MIN(e, vtcsCmpct.n);
    for (i = 0; i < procs.n; i++) {
      t.workProc = pPtr[i].procId;
      for (j = s; j < e; j++) {
        t.vertexId = vPtr[j].vertexId;
        t.sequenceId = vPtr[j].sequenceId;
        t.elementId = vPtr[j].elementId;
        array_cat(vertex, &a, &t, 1);
      }
    }
    s = e;
  }
  array_free(&vtcsCmpct);
  array_free(&procs);

  sarray_transfer(vertex, &a, workProc, 1, &cr);
  sarray_sort_2(vertex, a.ptr, a.n, vertexId, 1, elementId, 1, bfr);
  crystal_free(&cr);

  // create the map
  if (a.n == 0)
    return NULL;

  VToEMap *map = calloc(1, sizeof(VToEMap));
  map->elements = calloc(a.n, sizeof(ulong));

  uint nGIds = 1, prev = 0;
  vertex *aPtr = (vertex *)a.ptr;
  for (i = 1; i < a.n; i++) {
    if (aPtr[i].vertexId != aPtr[prev].vertexId)
      nGIds++;
    prev = i;
  }

  map->size = nGIds;
  map->globalIds = calloc(nGIds, sizeof(ulong));
  map->offsets = calloc(nGIds + 1, sizeof(ulong));

  map->elements[0] = aPtr[0].elementId;
  map->globalIds[0] = aPtr[0].vertexId;
  map->offsets[0] = 0;

  prev = 0;
  uint nOffsets = 0;
  for (i = 1; i < a.n; i++) {
    if (aPtr[i].vertexId != aPtr[prev].vertexId) {
      nOffsets++;
      map->globalIds[nOffsets] = aPtr[i].vertexId;
      map->offsets[nOffsets] = prev = i;
    }
    map->elements[i] = aPtr[i].elementId;
  }
  map->offsets[++nOffsets] = a.n;
  assert(nOffsets == nGIds);

  array_free(&a);

  return map;
}

// key must be present in globalIds
static int getPosition(VToEMap *map, ulong key) {
  ulong *globalIds = map->globalIds;

  int begin = 0;
  int end = map->size;
  int mid = 0;
  while (begin < end) {
    mid = (begin + end) / 2;

    if (key == globalIds[mid])
      return mid;
    else if (key < globalIds[mid])
      end = mid;
    else
      begin = mid;
  };

  if (globalIds[mid] != key)
    return -1;
  return mid;
}

void freeVToEMap(VToEMap *map) {
  free(map->globalIds);
  free(map->offsets);
  free(map->elements);
  free(map);
}

static int faceCheck(Mesh mesh, struct comm *c, buffer *bfr) {
  VToEMap *map = getVToEMap(mesh, c, bfr);

  sint nelt = mesh->nelt;
  sint ndim = mesh->nDim;

  int faces[GC_MAX_FACES][GC_MAX_FACE_VERTICES];
  if (ndim == 3)
    memcpy(faces, faces3D, GC_MAX_FACES * GC_MAX_FACE_VERTICES * sizeof(int));
  else
    memcpy(faces, faces2D, GC_MAX_FACES * GC_MAX_FACE_VERTICES * sizeof(int));

  Point ptr = mesh->elements.ptr;
  int nf = (ndim == 3) ? 6 : 4;
  int nfv = (ndim == 3) ? 4 : 2;
  int nv = (ndim == 3) ? 8 : 4;

  struct array shared;
  array_init(LongID, &shared, 200);

  int err = 0;

  int i, j, k, l;
  for (i = 0; i < nelt && err == 0; i++) {
    for (j = 0; j < nf && err == 0; j++) {
      shared.n = 0;

      for (k = 0; k < nfv; k++) {
        ulong globalId = ptr[i * nv + faces[j][k] - 1].globalId + 1;
        int indx = getPosition(map, globalId);
        assert(indx >= 0);
        LongID elemId;
        for (l = map->offsets[indx]; l < map->offsets[indx + 1]; l++) {
          elemId.id = map->elements[l];
          array_cat(LongID, &shared, &elemId, 1);
        }
      }

      sarray_sort(LongID, shared.ptr, shared.n, id, 1, bfr);

      ulong prev = 0;
      int ncount = 1;
      LongID *sptr = shared.ptr;
      for (l = 1; l < shared.n; l++) {
        if (sptr[l].id != sptr[prev].id) {
          if (ncount == 3) {
            err = 1;
            break;
          }
          prev = l;
          ncount = 1;
        } else
          ncount++;
      }

      if (ncount == 3) {
        err = 1;
        break;
      }
    }
  }

  array_free(&shared);
  freeVToEMap(map);

  return err;
}

static int elementCheck(Mesh mesh, struct comm *c, buffer *bfr) {
  uint nelt = mesh->nelt;
  uint ndim = mesh->nDim;
  int nv = (ndim == 3) ? 8 : 4;

  LongID globalIds[8];
  Point ptr = mesh->elements.ptr;
  uint i, j;
  int err = 0;
  for (i = 0; i < nelt && err == 0; i++) {
    for (j = 0; j < nv; j++)
      globalIds[j].id = ptr[i * nv + j].globalId + 1;

    sarray_sort(LongID, globalIds, nv, id, 1, bfr);

    for (j = 0; j < nv - 1; j++) {
      if (globalIds[j].id == globalIds[j + 1].id) {
        err = 1;
        break;
      }
    }
  }

  return err;
}

//==============================================================================
// C interface to find_conn
//
#define check_error(err, msg)                                                  \
  do {                                                                         \
    sint buf;                                                                  \
    comm_allreduce(&c, gs_int, gs_max, &err, 1, &buf);                         \
    if (err) {                                                                 \
      if (c.id == 0) {                                                         \
        printf("parCon error in %s.\n", msg);                                  \
        fflush(stdout);                                                        \
      }                                                                        \
      buffer_free(&bfr);                                                       \
      mesh_free(mesh);                                                         \
      comm_free(&c);                                                           \
      return err;                                                              \
    }                                                                          \
  } while (0)

static int transferBoundaryFaces(Mesh mesh, struct comm *c) {
  uint size = c->np;

  struct array *boundary = &mesh->boundary;
  BoundaryFace ptr = boundary->ptr;
  int nFaces = boundary->n;

  slong nelgt = mesh->nelgt;
  sint nelt = nelgt / size;
  sint nrem = nelgt - nelt * size;
  slong N = (size - nrem) * nelt;

  sint i;
  slong eid;
  for (i = 0; i < nFaces; i++) {
    eid = ptr[i].elementId;
    if (eid < N)
      ptr[i].proc = eid / nelt;
    else
      ptr[i].proc = (eid - N) / (nelt + 1) + size - nrem;
  }

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct boundary_t, boundary, proc, 1, &cr);
  crystal_free(&cr);

  return 0;
}

// Input:
//   nelt: Number of elements, nv: Number of vertices in an element
//   coord [nelt, nv, ndim]: Coordinates of elements vertices in preprocessor
//     ordering, nv = 8 if ndim == 3 (Hex) or nv = 4 if ndim = 2 (Quad).
// Output:
//   vtx[nelt, nv]: Global numbering of vertices of elements
int parrsb_conn_mesh(long long *vtx, double *coord, int nelt, int ndim,
                     long long *periodicInfo, int nPeriodicFaces, double tol,
                     MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  int verbose = 0;
  {
    const char *val = getenv("PARRSB_VERBOSE_LEVEL");
    if (val != NULL) verbose= atoi(val);  
  }

  if (c.id == 0 && verbose > 0) {
    printf("Running parCon ... \n");
    fflush(stdout);
  }

  double duration[8] = {0};
  const char *name[8] = {"transferBoundaryFaces", "findMinNbrDistance   ",
                         "findSegments         ", "setGlobalId          ",
                         "elementCheck         ", "faceCheck            ",
                         "matchPeriodicFaces   "};

  parrsb_barrier(&c);
  double tcon = comm_time();

  Mesh mesh;
  mesh_init(&mesh, nelt, ndim);

  slong out[2][1], wrk[2][1], in = nelt;
  comm_scan(out, &c, gs_long, gs_add, &in, 1, wrk);
  ulong start = out[0][0], nelgt = out[1][0];
  mesh->nelgt = mesh->nelgv = nelgt;

  int nvertex = mesh->nVertex;
  struct point_t p = {.origin = c.id};
  uint i, j, k, l;
  for (i = 0; i < nelt; i++) {
    for (k = 0; k < nvertex; k++) {
      j = PRE_TO_SYM_VERTEX[k];
      for (l = 0; l < ndim; l++)
        p.x[l] = coord[i * nvertex * ndim + j * ndim + l];
      p.elementId = start + i;
      p.sequenceId = nvertex * (start + i) + k;
      array_cat(struct point_t, &mesh->elements, &p, 1);
    }
  }

  struct boundary_t b;
  for (i = 0; i < nPeriodicFaces; i++) {
    b.elementId = periodicInfo[4 * i + 0] - 1;
    b.faceId = PRE_TO_SYM_FACE[periodicInfo[4 * i + 1] - 1];
    b.bc[0] = periodicInfo[4 * i + 2] - 1;
    b.bc[1] = PRE_TO_SYM_FACE[periodicInfo[4 * i + 3] - 1];
    array_cat(struct boundary_t, &mesh->boundary, &b, 1);
  }

  buffer bfr;
  buffer_init(&bfr, 1024);

  parrsb_barrier(&c);
  double t = comm_time();
  sint err = transferBoundaryFaces(mesh, &c);
  check_error(err, "transferBoundaryFaces");
  duration[0] = comm_time() - t;

  parrsb_barrier(&c);
  t = comm_time();
  err = findMinNeighborDistance(mesh);
  check_error(err, "findMinNeighborDistance");
  duration[1] = comm_time() - t;

  parrsb_barrier(&c);
  t = comm_time();
  err = findUniqueVertices(mesh, &c, tol, verbose, &bfr);
  check_error(err, "findSegments");
  duration[2] = comm_time() - t;

  parrsb_barrier(&c);
  t = comm_time();
  setGlobalID(mesh, &c);
  sendBack(mesh, &c, &bfr);
  duration[3] = comm_time() - t;

  parrsb_barrier(&c);
  t = comm_time();
  err = elementCheck(mesh, &c, &bfr);
  check_error(err, "elementCheck");
  duration[4] = comm_time() - t;

  parrsb_barrier(&c);
  t = comm_time();
  err = faceCheck(mesh, &c, &bfr);
  check_error(err, "faceCheck");
  duration[5] = comm_time() - t;

  parrsb_barrier(&c);
  t = comm_time();
  err = matchPeriodicFaces(mesh, &c, &bfr);
  check_error(err, "matchPeriodicFaces");
  duration[6] = comm_time() - t;

  // Copy output
  Point ptr = mesh->elements.ptr;
  for (i = 0; i < nelt; i++) {
    for (j = 0; j < nvertex; j++)
      vtx[i * nvertex + j] = ptr[i * nvertex + j].globalId + 1;
  }

  // Report timing info and finish
  parrsb_barrier(&c);
  tcon = comm_time() - tcon;
  if (c.id == 0) {
    printf("parCon (tol = %e) finished in %g s\n", tol, tcon);
    fflush(stdout);
  }

  double gmin[8], gmax[8], buf[8];
  for (unsigned i = 0; i < 8; i++)
    gmax[i] = gmin[i] = duration[i];
  comm_allreduce(&c, gs_double, gs_min, gmin, 8, buf);
  comm_allreduce(&c, gs_double, gs_max, gmax, 8, buf);

  if (c.id == 0 && verbose > 1) {
    for (unsigned i = 0; i < 7; i++)
      printf("%s: %e %e (min max)\n", name[i], gmin[i], gmax[i]);
    fflush(stdout);
  }

  buffer_free(&bfr);
  mesh_free(mesh);
  comm_free(&c);

  return err;
}

//=============================================================================
// Fortran interface
//
void fparrsb_conn_mesh(long long *vtx, double *coord, int *nelt, int *ndim,
                       long long *periodicInfo, int *nPeriodicFaces,
                       double *tol, MPI_Fint *fcomm, int *err) {
  *err = 1;
  MPI_Comm c = MPI_Comm_f2c(*fcomm);
  *err = parrsb_conn_mesh(vtx, coord, *nelt, *ndim, periodicInfo,
                          *nPeriodicFaces, *tol, c);
}

#undef GC_MAX_FACES
#undef GC_MAX_VERTICES
#undef GC_MAX_NEIGHBORS
#undef GC_MAX_FACE_VERTICES
#undef MIN
#undef MAX

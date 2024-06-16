#include "con-impl.h"

int PRE_TO_SYM_VERTEX[GC_MAX_VERTICES] = {0, 1, 3, 2, 4, 5, 7, 6};
int PRE_TO_SYM_FACE[GC_MAX_FACES] = {2, 1, 3, 0, 4, 5};
int NEIGHBOR_MAP[GC_MAX_VERTICES][GC_MAX_NEIGHBORS] = {
    {1, 2, 4}, {0, 3, 5}, {0, 3, 6}, {1, 2, 7},
    {0, 5, 6}, {1, 4, 7}, {2, 4, 7}, {3, 5, 6}};

double diff_sqr(double x, double y) { return (x - y) * (x - y); }

//==============================================================================
// Mesh struct
//
static struct mesh_t *mesh_init(uint nelt, unsigned ndim, double *coord,
                                long long *pinfo, uint npinfo,
                                const struct comm *c) {
  struct mesh_t *m = tcalloc(struct mesh_t, 1);
  m->nelt = nelt, m->ndim = ndim, m->nnbrs = ndim;
  m->nv = (ndim == 2) ? 4 : 8;

  slong out[2][1], wrk[2][1], in = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  ulong start = out[0][0];
  m->nelgt = out[1][0];

  uint nv = m->nv;
  array_init(struct point_t, &m->elements, nelt * nv);
  struct point_t p = {.origin = c->id};
  for (uint i = 0; i < nelt; i++) {
    for (uint k = 0; k < nv; k++) {
      uint j = PRE_TO_SYM_VERTEX[k];
      for (uint l = 0; l < ndim; l++)
        p.x[l] = coord[i * nv * ndim + j * ndim + l];
      p.elementId = start + i, p.sequenceId = nv * (start + i) + k;
      array_cat(struct point_t, &m->elements, &p, 1);
    }
  }

  array_init(struct boundary_t, &m->boundary, npinfo);
  struct boundary_t b;
  for (uint i = 0; i < npinfo; i++) {
    b.elementId = pinfo[4 * i + 0] - 1;
    b.faceId = PRE_TO_SYM_FACE[pinfo[4 * i + 1] - 1];
    b.bc[0] = pinfo[4 * i + 2] - 1;
    b.bc[1] = PRE_TO_SYM_FACE[pinfo[4 * i + 3] - 1];
    array_cat(struct boundary_t, &m->boundary, &b, 1);
  }

  return m;
}

static int mesh_free(struct mesh_t *m) {
  array_free(&m->elements), array_free(&m->boundary), free(m);
  return 0;
}

//==============================================================================
// Find the minimum distance between a vertex and its neighbors
//
static inline double distance_2d(struct point_t *a, struct point_t *b) {
  return diff_sqr(a->x[0], b->x[0]) + diff_sqr(a->x[1], b->x[1]);
}

static inline double distance_3d(struct point_t *a, struct point_t *b) {
  return distance_2d(a, b) + diff_sqr(a->x[2], b->x[2]);
}

int find_min_neighbor_distance(Mesh mesh) {
  struct point_t *p = (struct point_t *)mesh->elements.ptr;
  uint ndim = mesh->ndim;
  uint nv = mesh->nv;

  if (ndim < 2 || ndim > 3) return 1;

  uint i, j, k, neighbor;
  if (ndim == 3) {
    for (i = 0; i < mesh->elements.n; i += nv) {
      for (j = 0; j < nv; j++) {
        p[i + j].dx = SCALAR_MAX;
        for (k = 0; k < mesh->nnbrs; k++) {
          neighbor = NEIGHBOR_MAP[j][k];
          scalar d = distance_3d(&p[i + j], &p[i + neighbor]);
          p[i + j].dx = MIN(p[i + j].dx, d);
        }
      }
    }
  }

  if (ndim == 2) {
    for (i = 0; i < mesh->elements.n; i += nv) {
      for (j = 0; j < nv; j++) {
        p[i + j].dx = SCALAR_MAX;
        for (k = 0; k < mesh->nnbrs; k++) {
          neighbor = NEIGHBOR_MAP[j][k];
          scalar d = distance_2d(&p[i + j], &p[i + neighbor]);
          p[i + j].dx = MIN(p[i + j].dx, d);
        }
      }
    }
  }

  return 0;
}

//==============================================================================
// Global numbering
//
static int set_global_id(Mesh mesh, struct comm *c) {
  uint nPoints = mesh->elements.n;
  Point points = (struct point_t *)mesh->elements.ptr;

  sint bin = (nPoints > 0);
  struct comm nonZeroRanks;
  comm_split(c, bin, c->id, &nonZeroRanks);

  if (bin == 1) {
    slong count = 0;
    for (uint i = 0; i < nPoints; i++)
      if (points[i].ifSegment) count++;

    slong in = count, out[2][1], buf[2][1];
    comm_scan(out, &nonZeroRanks, gs_long, gs_add, &in, 1, buf);
    slong start = out[0][0];

    count = -1;
    for (uint i = 0; i < nPoints; i++) {
      if (points[i].ifSegment) count++;
      assert(start + count >= 0);
      points[i].globalId = start + count;
    }
  }

  comm_free(&nonZeroRanks);

  return 0;
}

int send_back(Mesh mesh, struct comm *c, buffer *bfr) {
  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct point_t, &mesh->elements, origin, 0, &cr);
  crystal_free(&cr);

  sarray_sort(struct point_t, mesh->elements.ptr, mesh->elements.n, sequenceId,
              1, bfr);

  return 0;
}

static int transfer_boundary_faces(Mesh mesh, struct comm *c) {
  uint size = c->np;

  struct array *boundary = &mesh->boundary;
  BoundaryFace ptr = (struct boundary_t *)boundary->ptr;
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

//==============================================================================
// C interface to find_conn
//
// Input:
//   nelt: Number of elements, nv: Number of vertices in an element
//   coord [nelt, nv, ndim]: Coordinates of elements vertices in preprocessor
//     ordering, nv = 8 if ndim == 3 (Hex) or nv = 4 if ndim = 2 (Quad).
// Output:
//   vtx[nelt, nv]: Global numbering of vertices of elements
int parrsb_conn_mesh(long long *vtx, double *coord, uint nelt, unsigned ndim,
                     long long *pinfo, int npinfo, double tol, MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  buffer bfr;
  buffer_init(&bfr, 1024);

  int verbose = 1;
  {
    const char *val = getenv("PARRSB_VERBOSE_LEVEL");
    if (val != NULL) verbose = atoi(val);
  }

  parrsb_print(&c, verbose, "Running parCon ...");

  parrsb_barrier(&c);
  double tall = comm_time(), t;

  double duration[8] = {0};
  const char *name[8] = {
      "transfer_boundary_faces    ", "find_min_neighbor_distance ",
      "find_unique_vertices       ", "set_global_id              ",
      "element_check              ", "face_check                 ",
      "match_periodic_faces       ", "copy_output                "};

  Mesh mesh = mesh_init(nelt, ndim, coord, pinfo, npinfo, &c);

  parrsb_print(&c, verbose - 1, "\t%s ...", name[0]);
  parrsb_barrier(&c), t = comm_time();
  transfer_boundary_faces(mesh, &c);
  duration[0] = comm_time() - t;

  parrsb_print(&c, verbose - 1, "\t%s ...", name[1]);
  parrsb_barrier(&c), t = comm_time();
  find_min_neighbor_distance(mesh);
  duration[1] = comm_time() - t;

  parrsb_print(&c, verbose - 1, "\t%s ...", name[2]);
  parrsb_barrier(&c), t = comm_time();
  find_unique_vertices(mesh, &c, tol, verbose - 1, &bfr);
  duration[2] = comm_time() - t;

  parrsb_print(&c, verbose - 1, "\t%s ...", name[3]);
  parrsb_barrier(&c), t = comm_time();
  set_global_id(mesh, &c);
  send_back(mesh, &c, &bfr);
  duration[3] = comm_time() - t;

#define check_error(call, msg)                                                 \
  {                                                                            \
    sint err = (call), wrk;                                                    \
    comm_allreduce(&c, gs_int, gs_max, &err, 1, &wrk);                         \
    if (err) {                                                                 \
      parrsb_print(&c, 1, msg, __FILE__, __LINE__);                            \
      buffer_free(&bfr), mesh_free(mesh), comm_free(&c);                       \
      return err;                                                              \
    }                                                                          \
  }

  parrsb_print(&c, verbose - 1, "\t%s ...", name[4]);
  parrsb_barrier(&c), t = comm_time();
  check_error(element_check(mesh, &c, &bfr), "\t%s:%d element_check failed.");
  duration[4] = comm_time() - t;

  parrsb_print(&c, verbose - 1, "\t%s ...", name[5]);
  parrsb_barrier(&c), t = comm_time();
  check_error(face_check(mesh, &c, &bfr), "\t%s:%d face_check failed.");
  duration[5] = comm_time() - t;

#undef check_error

  parrsb_print(&c, verbose - 1, "\t%s ...", name[6]);
  parrsb_barrier(&c), t = comm_time();
  match_periodic_faces(mesh, &c, verbose - 1, &bfr);
  duration[6] = comm_time() - t;

  parrsb_print(&c, verbose - 1, "\t%s ...", name[7]);
  parrsb_barrier(&c), t = comm_time();
  Point ptr = mesh->elements.ptr;
  for (uint i = 0; i < nelt; i++) {
    for (uint j = 0; j < mesh->nv; j++)
      vtx[i * mesh->nv + j] = ptr[i * mesh->nv + j].globalId + 1;
  }
  duration[7] = comm_time() - t;

  // Report timing info and finish
  {
    double gmin[8], gmax[8], buf[8];
    for (unsigned i = 0; i < 8; i++) gmax[i] = gmin[i] = duration[i];
    comm_allreduce(&c, gs_double, gs_min, gmin, 8, buf);
    comm_allreduce(&c, gs_double, gs_max, gmax, 8, buf);

    for (unsigned i = 0; i < 7; i++) {
      parrsb_print(&c, verbose - 1, "%s: %e %e (min max)", name[i], gmin[i],
                   gmax[i]);
    }
  }

  parrsb_barrier(&c);
  tall = comm_time() - tall;
  parrsb_print(&c, verbose, "parCon (tol = %e) finished in %g s", tol, tall);

  buffer_free(&bfr), mesh_free(mesh), comm_free(&c);

  return 0;
}

//=============================================================================
// Fortran interface
//
void fparrsb_conn_mesh(long long *vtx, double *coord, int *nelt, int *ndim,
                       long long *pinfo, int *npinfo, double *tol,
                       MPI_Fint *fcomm, int *err) {
  *err = 1;
  MPI_Comm c = MPI_Comm_f2c(*fcomm);
  *err = parrsb_conn_mesh(vtx, coord, *nelt, *ndim, pinfo, *npinfo, *tol, c);
}

#include "parrsb-impl.h"

#include <float.h>
#include <string.h>

static uint get_partition(const struct comm *const gc,
                          const struct comm *const lc) {
  // Find the partition id. A partition is a group of processors sharing the
  // same local communicator.
  sint out[2][1], wrk[2][1], root = (lc->id == 0);
  comm_scan(out, gc, gs_int, gs_add, &root, 1, wrk);
  sint part = out[0][0] * (lc->id == 0);
  comm_allreduce(lc, gs_int, gs_max, &part, 1, wrk);
  return part;
}

uint parrsb_get_neighbors(const struct array *const elems, const unsigned nv,
                          const struct comm *const gc,
                          const struct comm *const lc, buffer *bfr) {
  const uint n = elems->n;
  const uint size = elems->n * nv;

  struct vertex_t {
    ulong v;
    uint p, partition;
  };

  struct array vertices;
  array_init(struct vertex_t, &vertices, size);

  const struct rsb_element *const pe =
      (const struct rsb_element *const)elems->ptr;
  struct vertex_t vt = {.partition = get_partition(gc, lc)};
  for (uint i = 0; i < n; i++) {
    for (uint v = 0; v < nv; v++) {
      vt.v = pe[i].vertices[v], vt.p = vt.v % gc->np;
      array_cat(struct vertex_t, &vertices, &vt, 1);
    }
  }

  struct crystal cr;
  crystal_init(&cr, gc);

  sarray_transfer(struct vertex_t, &vertices, p, 1, &cr);
  sarray_sort(struct vertex_t, vertices.ptr, vertices.n, v, 1, bfr);

  struct array neighbors;
  array_init(struct vertex_t, &neighbors, vertices.n * 27);

  const struct vertex_t *const pv = (const struct vertex_t *const)vertices.ptr;
  uint s = 0;
  while (s < vertices.n) {
    uint e = s + 1;
    while (e < vertices.n && pv[s].v == pv[e].v) e++;
    for (uint i = s; i < e; i++) {
      struct vertex_t vt = pv[i];
      for (uint j = s; j < e; j++) {
        vt.partition = pv[j].partition;
        array_cat(struct vertex_t, &neighbors, &vt, 1);
      }
    }
    s = e;
  }
  array_free(&vertices);

  sarray_transfer(struct vertex_t, &neighbors, p, 0, &cr);
  crystal_free(&cr);
  sarray_sort(struct vertex_t, neighbors.ptr, neighbors.n, partition, 0, bfr);

  // Now, extract out different partition ids found locally into an array.
  struct unique_t {
    uint p, partition;
  };

  struct array unique;
  array_init(struct unique_t, &unique, 27);

  if (neighbors.n > 0) {
    const struct vertex_t *const pn =
        (const struct vertex_t *const)neighbors.ptr;
    struct unique_t ut = {.partition = pn[0].partition,
                          .p = pn[0].partition % lc->np};
    array_cat(struct unique_t, &unique, &ut, 1);
    for (uint i = 1; i < neighbors.n; i++) {
      if (pn[i].partition > pn[i - 1].partition) {
        ut.partition = pn[i].partition, ut.p = ut.partition % lc->np;
        array_cat(struct unique_t, &unique, &ut, 1);
      }
    }
  }
  array_free(&neighbors);

  crystal_init(&cr, lc);
  sarray_transfer(struct unique_t, &unique, p, 0, &cr);
  crystal_free(&cr);

  sarray_sort(struct unique_t, unique.ptr, unique.n, partition, 0, bfr);
  sint un = 0;
  if (unique.n > 0) {
    un = 1;
    struct unique_t *pu = (struct unique_t *)unique.ptr;
    for (uint i = 1; i < unique.n; i++) {
      if (pu[i].partition > pu[un - 1].partition) pu[un] = pu[i], un++;
    }
  }
  array_free(&unique);

  sint wrk;
  comm_allreduce(lc, gs_int, gs_add, &un, 1, &wrk);
  assert(un >= 1);

  return un - 1;
}

static struct array pgeom;
static buffer bfr;
static uint pgeom_initialized = 0;
static uint nv = 0;
static uint level = 0;

struct pgeom_t {
  uint partition, level;
  double centroid[3], min[3], max[3];
  uint p;
};

void parrsb_dump_stats_start(const uint nv_) {
  if (pgeom_initialized) return;

  nv = nv_;
  level = 0;
  array_init(struct pgeom_t, &pgeom, 1024);
  buffer_init(&bfr, 1024);

  pgeom_initialized = 1;
}

void parrsb_dump_stats(const struct comm *const gc, const struct comm *const lc,
                       const struct array *const elems, buffer *bfr) {
  if (!pgeom_initialized) return;

  const struct rsb_element *const pe =
      (const struct rsb_element *const)elems->ptr;

  // Find the centroid and the bounding box of the partition.
  double centroid[3] = {0.0, 0.0, 0.0};
  double max[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
  double min[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
  const uint n = elems->n;
  const unsigned ndim = (nv == 8) ? 3 : 2;
  for (uint e = 0; e < n; e++) {
    for (uint d = 0; d < ndim; d++) {
      double c = pe[e].coord[d];
      centroid[d] += c;
      max[d] = (max[d] < c) ? c : max[d];
      min[d] = (min[d] > c) ? c : min[d];
    }
  }
  for (uint d = 0; d < ndim; d++) centroid[d] /= n;

  double wrk[3];
  comm_allreduce(lc, gs_double, gs_min, min, ndim, wrk);
  comm_allreduce(lc, gs_double, gs_max, max, ndim, wrk);
  comm_allreduce(lc, gs_double, gs_add, centroid, ndim, wrk);
  for (uint d = 0; d < ndim; d++) centroid[d] /= lc->np;

  // Partition root accumulates the partition geometry.
  level++;
  struct pgeom_t pg = {.partition = get_partition(gc, lc),
                       .level = level,
                       .centroid = {centroid[0], centroid[1], centroid[2]},
                       .max = {max[0], max[1], max[2]},
                       .min = {min[0], min[1], min[2]},
                       .p = 0};
  if (lc->id == 0) array_cat(struct pgeom_t, &pgeom, &pg, 1);
}

void parrsb_dump_stats_end(const struct comm *const gc, const char *prefix) {
  if (!pgeom_initialized) return;

  const uint size = strnlen(prefix, 64);
  assert(size < 64 && "Prefix must be less than 64 characters.");

  // Send all the data to global root.
  struct crystal cr;
  crystal_init(&cr, gc);
  sarray_transfer(struct pgeom_t, &pgeom, p, 0, &cr);
  crystal_free(&cr);

  // Sort by level first, then by partition id.
  sarray_sort_2(struct pgeom_t, pgeom.ptr, pgeom.n, level, 0, partition, 0,
                &bfr);

  if (gc->id == 0) {
    const char name[BUFSIZ];
    snprintf((char *)name, BUFSIZ, "%s_partition_geom_p%06d.txt", prefix,
             gc->np);

    FILE *fp = fopen(name, "w");
    if (!fp) {
      fprintf(stderr, "Failed to open %s for writing.\n", name);
      exit(EXIT_FAILURE);
    }

    fprintf(fp, "%zu\n", pgeom.n);
    fprintf(fp, "level partition centroid[0] centroid[1] centroid[2] min[0] "
                "min[1] min[2] max[0] max[1] max[2]\n");
    const struct pgeom_t *const pg = (const struct pgeom_t *const)pgeom.ptr;
    for (uint i = 0; i < pgeom.n; i++) {
      fprintf(fp, "%u %u %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", pg[i].level,
              pg[i].partition, pg[i].centroid[0], pg[i].centroid[1],
              pg[i].centroid[2], pg[i].min[0], pg[i].min[1], pg[i].min[2],
              pg[i].max[0], pg[i].max[1], pg[i].max[2]);
    }
    fclose(fp);
  }

  array_free(&pgeom);
  buffer_free(&bfr);

  pgeom_initialized = nv = level = 0;
}

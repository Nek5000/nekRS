#include "con-impl.h"

static struct comm COMM_NULL = {.id = 0, .np = 0, .c = MPI_COMM_NULL};

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

  if (n < 2) return;
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
      if (j < ir && get(ra, j) < get(ra, j + 1)) j++;
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

  free(rra);
}

#define tuple_sort(T, arr, n, index)                                           \
  tuple_sort_((void *)arr, n, sizeof(T), offsetof(T, index))

static void sort_segments_local(struct array *local, int dim) {
  uint npts = local->n;
  if (npts == 0) return;

  struct point_t *const pts = (struct point_t *const)local->ptr;
  uint s = 0, e;
  while (s < npts) {
    for (e = s + 1; e < npts && pts[e].ifSegment == 0; e++)
      ;

    if (s < npts - 1 && e - s > 1) {
      switch (dim) {
      case 0: tuple_sort(struct point_t, &pts[s], e - s, x[0]); break;
      case 1: tuple_sort(struct point_t, &pts[s], e - s, x[1]); break;
      case 2: tuple_sort(struct point_t, &pts[s], e - s, x[2]); break;
      default: break;
      }
    }

    uint sum = 0;
    for (uint i = s; i < e; i++) sum += pts[i].ifSegment, pts[i].ifSegment = 0;

    if (sum > 0) pts[s].ifSegment = 1;

    s = e;
  }
}

static void sort_segments_shared_aux(struct array *arr, int dim, struct comm *c,
                                     int verbose, buffer *bfr) {
  parrsb_print(c, verbose, "\t\t\t\tsss_aux_parallel_sort: ...\n");
  switch (dim) {
  case 0:
    parallel_sort(struct point_t, arr, x[0], gs_double, 0, 1, c, bfr);
    break;
  case 1:
    parallel_sort(struct point_t, arr, x[1], gs_double, 0, 1, c, bfr);
    break;
  case 2:
    parallel_sort(struct point_t, arr, x[2], gs_double, 0, 1, c, bfr);
    break;
  default: break;
  }
  parrsb_print(c, verbose, "\t\t\t\tsss_aux_parallel_sort: done.\n");

  // Mark the first point of the segment to have ifSegment = 1 and zero out
  // everything else.
  struct point_t *const pts = (struct point_t *const)arr->ptr;
  for (uint i = 0; i < arr->n; i++) pts[i].ifSegment = 0;

  sint wrk;
  sint rank = (arr->n > 0) ? c->id : c->np;
  comm_allreduce(c, gs_int, gs_min, &rank, 1, &wrk);

  if ((sint)c->id == rank) pts[0].ifSegment = 1;

  parrsb_print(c, verbose, "\t\t\t\tsss_aux_mark_first_point: done.");
}

static uint find_bin_scan(const sint sum, const struct comm *c,
                          const int verbose, buffer *bfr) {
  sint out[2][1], wrk[2][1], in = sum;
  comm_scan(out, c, gs_int, gs_add, &in, 1, wrk);
  return out[0][0];
}

static uint find_bin_gs(const slong id, const struct comm *c, const int verbose,
                        buffer *bfr) {
  slong gid = id + 1;
  struct gs_data *gsh = gs_setup(&gid, 1, c, 0, gs_crystal_router, verbose);
  parrsb_print(c, verbose, "\t\t\tsss_gs_setup: done.");
  sint bin = c->id;
  gs(&bin, gs_int, gs_min, 0, gsh, bfr);
  gs_free(gsh);

  return bin;
}

static uint find_bin_cr(const slong id, const struct comm *c, const int verbose,
                        buffer *bfr) {
  struct gid_t {
    ulong id;
    uint proc, procm;
  };

  struct array arr;
  array_init(struct gid_t, &arr, 1);

  struct gid_t gid = {.id = id, .proc = id % c->np, .procm = c->id};
  array_cat(struct gid_t, &arr, &gid, 1);

  struct crystal cr;
  crystal_init(&cr, c);

  sarray_transfer(struct gid_t, &arr, proc, 1, &cr);
  if (arr.n > 0) {
    sarray_sort_2(struct gid_t, arr.ptr, arr.n, id, 1, procm, 0, bfr);
    struct gid_t *pa = (struct gid_t *)arr.ptr;
    uint s = 0;
    while (s < arr.n) {
      uint e = s + 1;
      for (; e < arr.n && pa[s].id == pa[e].id; e++) pa[e].procm = pa[s].procm;
      s = e;
    }
  }
  sarray_transfer(struct gid_t, &arr, proc, 0, &cr);

  crystal_free(&cr);

  assert(arr.n == 1);
  struct gid_t *pa = (struct gid_t *)arr.ptr;
  uint procm = pa[0].procm;
  array_free(&arr);

  return procm;
}

static void sort_segments_shared(struct array *shared, int dim, struct comm *c,
                                 int verbose, buffer *bfr) {
  // Each process can only have at most a single ifSegment = 1 in shared
  // array. Otherwise, we can always move the segments into the local segments
  // array till we end up in such a configuration. Let's first check for this
  // condition and find how many distinct global ids (either 1 or 2 ids
  // respectively depending on if there is 0 or 1 ifSegment = 1 points) are
  // residing on this MPI process. While doing so, we will also split shared
  // array to at most two segments each corresponding to a distinct global id.
  // Algorithm listed below discard MPI processes with shared array size 0.
  struct array segments[2];
  array_init(struct point_t, &segments[0], (shared->n + 1) / 2);
  array_init(struct point_t, &segments[1], (shared->n + 1) / 2);
  slong gids[2] = {INT_MIN, INT_MIN};

  sint sum = 0, ngids = 0;
  if (shared->n > 0) {
    struct point_t *pts = (struct point_t *)shared->ptr;
    gids[0] = pts[0].globalId, ngids++;
    sum += pts[0].ifSegment;
    array_cat(struct point_t, &segments[ngids - 1], &pts[0], 1);
    for (uint i = 1; i < shared->n; i++) {
      if (pts[i].ifSegment > 0) gids[1] = pts[i].globalId, ngids++;
      sum += pts[i].ifSegment;
      array_cat(struct point_t, &segments[ngids - 1], &pts[i], 1);
    }
  }
  assert(sum <= 1);
  assert(ngids <= 1 || (ngids == 2 && gids[1] == gids[0] + 1));
  parrsb_print(c, verbose, "\t\t\tsss_local: done.");

  // Algorithm to be used for finding the bin id for segmented shared sort.
  // Default (algo = 0) is the scan. algo = 1 is gs with gs_crystal_router.
  // algo = 2 is a custom crystal router implementation.
  int algo = 0;
  char *val = getenv("PARRSB_FIND_BIN_ALGO");
  if (val) algo = atoi(val);
  assert(algo >= 0 && algo <= 2);

  // We sort the shared segments in two phases. All the segments having an even
  // global id are sorted first and then the segments having an odd global id
  // are sorted. This is done to avoid same process having to work on both the
  // global ids (if ngids = 2) it owns at the same time.
  for (int parity = 0; parity < 2; parity++) {
    parrsb_print(c, verbose, "\t\t\tsss_parity_%d: ...", parity);
    int index = INT_MIN;
    if (gids[0] >= 0 && (gids[0] % 2 == parity))
      index = 0;
    else if (gids[1] >= 0 && (gids[1] % 2 == parity))
      index = 1;

    struct comm active, seg;
    comm_split(c, index >= 0, c->id, &active);
    if (index >= 0) {
      assert(gids[index] >= 0);
      sint bin = -1;
      if (algo == 0) {
        uint off = (ngids == 1 && sum == 1) || (ngids == 2 && index == 1);
        bin = find_bin_scan(sum, &active, verbose - 1, bfr) + off;
      } else if (algo == 1) {
        bin = find_bin_gs(gids[index], &active, verbose - 1, bfr);
      } else if (algo == 2) {
        bin = find_bin_cr(gids[index], &active, verbose - 1, bfr);
      }
      parrsb_print(&active, verbose,
                   "\t\t\tsss_find_bin_algo_%d_parity_%d: done.", algo, parity);
      assert(bin >= 0 && bin <= (sint)active.np);

      // index >= 0 --> gids[index] >= 0 --> segments[index].n > 0
      comm_split(&active, bin, active.id, &seg);
      sort_segments_shared_aux(&segments[index], dim, &seg, verbose - 1, bfr);
      comm_free(&seg);
      parrsb_print(&active, verbose, "\t\t\tsss_aux_%d: done.", parity);
    }
    comm_free(&active);
    parrsb_print(c, verbose, "\t\t\tsss_parity_%d: done.", parity);
  }
  parrsb_print(c, verbose, "\t\t\tsss_shared: done.");

  // Combine the segments after sorting.
  shared->n = 0;
  array_cat(struct point_t, shared, segments[0].ptr, segments[0].n);
  array_cat(struct point_t, shared, segments[1].ptr, segments[1].n);
  array_free(&segments[0]), array_free(&segments[1]);
}

static int talk_to_neighbor(struct point_t *pnt, const struct array *arr,
                            int dir, const struct comm *c) {
  assert(dir == -1 || dir == 1);

  if (c->np <= 1) return 0;

  struct comm active;
  comm_split(c, arr->n > 0, c->id, &active);

  if (arr->n == 0) {
    comm_free(&active);
    return 0;
  }

  struct array tmp;
  array_init(struct point_t, &tmp, 1);

  struct point_t *pts = (struct point_t *)arr->ptr;
  sint dest = (sint)c->id + dir;
  if (dest >= 0 && dest < (sint)c->np) {
    struct point_t p = (dir == 1) ? pts[arr->n - 1] : pts[0];
    p.proc = dest;
    array_cat(struct point_t, &tmp, &p, 1);
  }

  struct crystal cr;
  crystal_init(&cr, &active);
  sarray_transfer(struct point_t, &tmp, proc, 1, &cr);
  crystal_free(&cr);

  if (tmp.n == 0) return 0;

  pts = (struct point_t *)tmp.ptr, pnt[0] = pts[0];
  array_free(&tmp), comm_free(&active);
  return 1;
}

static void find_segments(struct array *arr, int i, scalar tol2,
                          const struct comm *c) {
  struct point_t *pts = (struct point_t *)arr->ptr;
  for (uint j = 1; j < arr->n; j++) {
    scalar d = diff_sqr(pts[j].x[i], pts[j - 1].x[i]);
    scalar dx = MIN(pts[j].dx, pts[j - 1].dx) * tol2;

    if (d > dx) pts[j].ifSegment = 1;
  }

  struct point_t pnt;
  int npts = talk_to_neighbor(&pnt, arr, 1 /* send */, c);
  if (npts > 0) { // npts > 0 --> arr->n > 0
    scalar d = diff_sqr(pnt.x[i], pts[0].x[i]);
    scalar dx = MIN(pnt.dx, pts[0].dx) * tol2;
    if (d > dx) pts[0].ifSegment = 1;
  }
}

static inline void remove_marked(struct array *arr) {
  struct point_t *pts = (struct point_t *)arr->ptr;
  uint count = 0;
  for (uint i = 0; i < arr->n; i++) {
    if (pts[i].ifSegment != -1) pts[count] = pts[i], count++;
  }
  arr->n = count;
}

static void separate_local_segments(struct array *local, struct array *shared,
                                    const struct comm *c) {
  // Find the first non-zero `ifSegment` value.
  struct point_t *pts = (struct point_t *)shared->ptr;
  uint s = 0;
  for (; s < shared->n && pts[s].ifSegment == 0; s++)
    ;
  sint lcheck = 0;
  while (s < shared->n) {
    uint e = s + 1;
    for (; e < shared->n && pts[e].ifSegment == 0; e++)
      ;
    if (e < shared->n) {
      for (uint i = s; i < e; i++) {
        array_cat(struct point_t, local, &pts[i], 1);
        // Mark the point to be removed from shared array.
        pts[i].ifSegment = -1;
      }
    } else if (e == shared->n) {
      // We reached end of the array without getting to a new segment.  This
      // means we have to talk to neighbor process. We will call
      // talk_to_neighbor if at least one of the process has to talk to the
      // neighbor.
      lcheck = 1;
    }
    s = e;
  }

  sint check = lcheck, wrk;
  comm_allreduce(c, gs_int, gs_add, &check, 1, &wrk);
  if (check) {
    // Bring the first point from next process. Check if `ifSegment` value
    // of that point is a 1 or a 0. If it is a 1, add the current range to
    // the local segment.
    struct point_t pnt;
    int npts = talk_to_neighbor(&pnt, shared, -1 /* recv */, c);
    if (lcheck && (npts == 0 || (npts > 0 && pnt.ifSegment == 1))) {
      for (uint i = s; i < shared->n; i++) {
        array_cat(struct point_t, local, &pts[i], 1);
        // Mark the point to be removed from shared array.
        pts[i].ifSegment = -1;
      }
    }
  }
  remove_marked(shared);
}

static slong number_segments(struct array *local, struct array *shared,
                             const struct comm *c) {
  // First number the local segments.
  struct point_t *pts = (struct point_t *)local->ptr;
  uint lcnt = 0;
  for (uint i = 0; i < local->n; i++) {
    if (pts[i].ifSegment) lcnt++;
  }

  slong out[2][1], wrk[2][1], in = lcnt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong ls = out[0][0], lt = out[1][0];

  ls--;
  for (uint i = 0; i < local->n; i++) {
    if (pts[i].ifSegment) ls++;
    assert(ls >= 0);
    pts[i].globalId = ls;
  }

  uint scnt = 0;
  pts = (struct point_t *)shared->ptr;
  for (uint i = 0; i < shared->n; i++) {
    if (pts[i].ifSegment) scnt++;
  }

  in = scnt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong ss = out[0][0], st = out[1][0];

  ss = lt + ss, ss--;
  for (uint i = 0; i < shared->n; i++) {
    if (pts[i].ifSegment) ss++;
    assert(ss >= lt);
    pts[i].globalId = ss;
  }

  return st + lt;
}

static void number_points(struct array *elems, const struct array *local,
                          const struct array *shared, const struct comm *c,
                          buffer *bfr) {
  // First number local points and then number shared points.
  slong out[2][1], wrk[2][1], in = local->n;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong s = out[0][0], nl = out[1][0];

  struct point_t *pts = (struct point_t *)local->ptr;
  for (uint i = 0; i < local->n; i++) pts[i].pntid = s + i;

  in = shared->n;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  s = out[0][0] + nl;

  pts = (struct point_t *)shared->ptr;
  for (uint i = 0; i < shared->n; i++) pts[i].pntid = s + i;

  // Copy everything back to elements array.
  elems->n = 0;
  array_cat(struct point_t, elems, local->ptr, local->n);
  array_cat(struct point_t, elems, shared->ptr, shared->n);

  // Now do a sort to load balance.
  parallel_sort(struct point_t, elems, pntid, gs_long, 0, 1, c, bfr);
}

int find_unique_vertices(Mesh mesh, struct comm *c, scalar tol, int verbose,
                         buffer *bfr) {
  scalar tol2 = tol * tol;
  int ndim = mesh->ndim;

  // Initialize globalId and ifSegment of each point to zero. These will be
  // updated as the gencon algorithm progress. ifSegment = 1 denotes the start
  // of a segment. A segment is a set of points which have the same global id
  // (or a set of points which we can't distinguish yet). Initially, all the
  // points are a single segment.
  struct array *elems = &mesh->elements;
  struct point_t *pts = (struct point_t *)elems->ptr;
  for (uint i = 0; i < elems->n; i++) pts[i].ifSegment = pts[i].globalId = 0;

  slong npts = elems->n, wrk;
  comm_allreduce(c, gs_long, gs_add, &npts, 1, &wrk);

  // Initialize shared and local arrays and then copy all points in `elems`
  // array to shared array first. Shared array contains only the segments which
  // are split across two or more processes. This means that the shared array
  // on a given process can only have at most one ifSegment value set to 1.
  // Initially all points are just one segment. Points will be removed from
  // shared array and added to local array as we progress in the algorithm.
  struct array shared, local;
  array_init(struct point_t, &shared, elems->n);
  array_init(struct point_t, &local, elems->n);
  array_cat(struct point_t, &shared, elems->ptr, elems->n);

  for (int t = 0; t < ndim; t++) {
    for (int d = 0; d < ndim; d++) {
      // Sort both local and shared segments.
      parrsb_print(c, verbose - 1, "\t\tsort_shared_segments ...");
      sort_segments_shared(&shared, d, c, verbose - 1, bfr);
      parrsb_print(c, verbose - 1, "\t\tsort_local_segments ...");
      sort_segments_local(&local, d);

      // Find segments in local and shared segments now.
      parrsb_print(c, verbose - 1, "\t\tfind_shared_segments ...");
      find_segments(&shared, d, tol2, c);
      parrsb_print(c, verbose - 1, "\t\tfind_local_segments ...");
      find_segments(&local, d, tol2, &COMM_NULL);

      // Separate local segments from the shared segments.
      parrsb_print(c, verbose - 1, "\t\tseparate_local_segments ...");
      separate_local_segments(&local, &shared, c);

      // Number the segments.
      parrsb_print(c, verbose - 1, "\t\tnumber_segments ...");
      slong nseg = number_segments(&local, &shared, c);
      parrsb_print(c, verbose, "\tlocglob: %d %d %lld %lld", t + 1, d + 1, nseg,
                   npts);
    }
  }
  // Number points consecutively -- shared points after local and then load
  // balance.
  parrsb_print(c, verbose - 1, "\tnumber_points_and_load_balance ...");
  number_points(elems, &local, &shared, c, bfr);
  array_free(&shared), array_free(&local);

  return 0;
}

#include <math.h>
#include <stdlib.h>

#include <gencon-impl.h>
#include <sort.h>

static void initSegments(Mesh mesh, struct comm *c) {
  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  // Initialize globalId and ifSegment
  slong out[2][1], buf[2][1], in[1];
  in[0] = nPoints;
  comm_scan(out, c, gs_long, gs_add, in, 1, buf);
  slong start = out[0][0];
  uint i;
  for (i = 0; i < nPoints; i++) {
    points[i].ifSegment = 0;
    points[i].globalId = start + i;
  }

  // First rank with nPoints > 0 set ifSegment = 1
  sint rank = c->id;
  if (nPoints == 0)
    rank = c->np;

  comm_allreduce(c, gs_int, gs_min, &rank, 1, buf);

  if (c->id == rank)
    points[0].ifSegment = 1;
}

static int sortSegments(Mesh mesh, struct comm *c, int dim, buffer *bfr) {
  sint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  sint s = 0, e;
  while (s < nPoints) {
    // find the length of the segment
    for (e = s + 1; e < nPoints && points[e].ifSegment == 0; e++)
      ;

    // sort start to end based on dim
    switch (dim) {
    case 0:
      sarray_sort(struct Point_private, &points[s], e - s, x[0], 3, bfr);
      break;
    case 1:
      sarray_sort(struct Point_private, &points[s], e - s, x[1], 3, bfr);
      break;
    case 2:
      sarray_sort(struct Point_private, &points[s], e - s, x[2], 3, bfr);
      break;
    default:
      break;
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

static int findLocalSegments(Mesh mesh, struct comm *c, int i,
                             GenmapScalar tolSquared) {
  Point pts = mesh->elements.ptr;
  sint npts = mesh->elements.n;
  int nDim = mesh->nDim;

  sint j;
  for (j = 1; j < npts; j++) {
    GenmapScalar d = sqrDiff(pts[j].x[i], pts[j - 1].x[i]);

    GenmapScalar dx = min(pts[j].dx, pts[j - 1].dx) * tolSquared;

    if (d > dx)
      pts[j].ifSegment = 1;
  }

  sint rank = c->id;
  sint size = c->np;

  struct Point_private lastp = pts[npts - 1];
  lastp.proc = (rank + 1) % size;

  struct array arr;
  array_init(struct Point_private, &arr, 1);
  array_cat(struct Point_private, &arr, &lastp, 1);

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct Point_private, &arr, proc, 1, &cr);
  crystal_free(&cr);

  uint n = arr.n;
  assert(n == 1);
  lastp = ((struct Point_private *)arr.ptr)[0];

  if (rank > 0) {
    GenmapScalar d = sqrDiff(lastp.x[i], pts->x[i]);
    GenmapScalar dx = min(lastp.dx, pts->dx) * tolSquared;
    if (d > dx)
      pts->ifSegment = 1;
  }

  array_free(&arr);

  return 0;
}

static int mergeSegments(Mesh mesh, struct comm *c, buffer *bfr) {
  uint npoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  int n, ifseg = 0;
  uint sendn = 0;
  for (n = 0; n < npoints; n++)
    if (points[n].ifSegment == 1) {
      ifseg = 1;
      sendn = n;
      break;
    }

  // TODO: comm_scan
  sint out[2][1], buf[2][1], in[1];
  in[0] = ifseg * c->id;
  comm_scan(out, c, gs_int, gs_max, in, 1, buf);
  sint rank = out[0][0];

  // If rank > 0, send i = 0,... n-1 where points[i].ifSegment == 0 to
  // rank with previous ifSegment == 1
  for (n = 0; n < sendn; n++)
    points[n].proc = rank;
  for (; n < npoints; n++)
    points[n].proc = c->id;

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct Point_private, &mesh->elements, proc, 0, &cr);
  crystal_free(&cr);

  sarray_sort(struct Point_private, mesh->elements.ptr, mesh->elements.n,
              globalId, 1, bfr);

  return 0;
}

slong countSegments(Mesh mesh, struct comm *c) {
  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  sint count = 0, i;
  for (i = 0; i < nPoints; i++)
    if (points[i].ifSegment > 0)
      count++;

  slong in, buf[2][1];
  in = count;
  comm_allreduce(c, gs_long, gs_add, &in, 1, buf);
  return in;
}

struct schedule {
  int dim;
  slong segments;
};

#define sort_by_coord(mesh, c, xa, xb, xc, bfr)                                \
  do {                                                                         \
    parallel_sort(struct Point_private, &(mesh->elements), x[xa], gs_scalar,   \
                  bin_sort, 0, c, bfr);                                        \
    uint nPoints = mesh->elements.n;                                           \
    Point points = mesh->elements.ptr;                                         \
                                                                               \
    int nDim = mesh->nDim;                                                     \
    if (nDim == 3)                                                             \
      sarray_sort_3(struct Point_private, points, nPoints, x[xa], 3, x[xb], 3, \
                    x[xc], 3, bfr);                                            \
    else if (nDim == 2)                                                        \
      sarray_sort_2(struct Point_private, points, nPoints, x[xa], 3, x[xb], 3, \
                    bfr);                                                      \
  } while (0)

#define segments_by_coord(cnt, sched, mesh, c, xa, tolSquared, bfr)            \
  do {                                                                         \
    sched[cnt].dim = xa;                                                       \
    initSegments(mesh, c);                                                     \
    findLocalSegments(mesh, c, xa, tolSquared);                                \
    sched[cnt].segments = -countSegments(mesh, c);                             \
  } while (0)

int findScheduleAndSort(struct schedule sched[3], Mesh mesh, struct comm *c,
                        GenmapScalar tolSquared, int verbose, buffer *bfr) {
  sort_by_coord(mesh, c, 0, 1, 2, bfr);
  segments_by_coord(0, sched, mesh, c, 0, tolSquared, bfr);

  int nDim = mesh->nDim;

  if (nDim == 3) {
    sort_by_coord(mesh, c, 1, 2, 0, bfr);
    segments_by_coord(1, sched, mesh, c, 1, tolSquared, bfr);

    sort_by_coord(mesh, c, 2, 0, 1, bfr);
    segments_by_coord(2, sched, mesh, c, 2, tolSquared, bfr);
  } else {
    sort_by_coord(mesh, c, 1, 0, 2, bfr);
    segments_by_coord(0, sched, mesh, c, 1, tolSquared, bfr);
  }

  sarray_sort(struct schedule, sched, nDim, segments, 1, bfr);

  if (nDim == 2) {
    printf("Not implemented.\n");
    exit(1);
  } else {
    switch (sched[0].dim) {
    case 0:
      sort_by_coord(mesh, c, 0, 1, 2, bfr);
      break;
    case 1:
      sort_by_coord(mesh, c, 1, 2, 0, bfr);
      break;
    case 2:
      sort_by_coord(mesh, c, 2, 0, 1, bfr);
      break;
    default:
      break;
    }
  }

  int i;
  for (i = 1; i < nDim; i++)
    sched[i].dim = (sched[0].dim + i) % nDim;

  return 0;
}

int findSegments(Mesh mesh, struct comm *c, GenmapScalar tol, int verbose,
                 buffer *bfr) {
  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;
  int nDim = mesh->nDim;

  int bin = (nPoints > 0);

  struct comm nonZeroRanks, dup;
  genmap_comm_split(c, bin, c->id, &nonZeroRanks);
  comm_dup(&dup, &nonZeroRanks);

  GenmapScalar tolSquared = tol * tol;
  struct schedule sched[3];
  sched[0].dim = 0, sched[1].dim = 1, sched[2].dim = 2;
  sort_by_coord(mesh, c, 0, 1, 2, bfr);
  initSegments(mesh, c);

  int t, d, merge = 1, err = 0;
  for (t = 0; t < nDim && err == 0; t++) {
    for (d = 0; d < nDim; d++) {
      int dim = sched[d].dim;

      if (bin > 0) {
        sortSegments(mesh, &nonZeroRanks, dim, bfr);
        findLocalSegments(mesh, &nonZeroRanks, dim, tolSquared);

        slong segments = countSegments(mesh, &nonZeroRanks);
        int rank = nonZeroRanks.id;
        if (rank == 0 && verbose > 0)
          printf("\tlocglob: %d %d %lld\n", t + 1, dim + 1, segments);

        if (merge > 0) {
          mergeSegments(mesh, &nonZeroRanks, bfr);
          merge = 0;

          nPoints = mesh->elements.n;
          bin = (nPoints > 0);

          comm_free(&nonZeroRanks);
          genmap_comm_split(&dup, bin, nonZeroRanks.id, &nonZeroRanks);
        }
      }
    }
  }

  comm_free(&dup);
  comm_free(&nonZeroRanks);

  return err;
}

#include <stdlib.h>
#include <math.h>

#include <gencon-impl.h>
#include <sort.h>

static int sortSegments(Mesh mesh, struct comm *c, int dim, buffer *bfr) {
  sint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  sint s = 0, e;
  while (s < nPoints) {
    // find the length of the segment
    for (e = s + 1; e < nPoints && points[e].ifSegment == 0; e++);

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

    s = e;
  }

  comm_barrier(c);

  // Set globalId
  slong out[2][1], buf[2][1], in[1];
  in[0] = nPoints;
  comm_scan(out, c, gs_long, gs_add, in, 1, buf);
  slong start = out[0][0];

  for (e = 0; e < nPoints; e++)
    points[e].globalId = start + e;

  return 0;
}

static int findLocalSegments(Mesh mesh, int i, GenmapScalar tolSquared) {
  Point pts = mesh->elements.ptr;
  sint npts = mesh->elements.n;
  int nDim = mesh->nDim;

  sint j;
  for (j = 0; j < npts - 1; j++) {
    GenmapScalar d = sqrDiff(pts[j].x[i], pts[j + 1].x[i]);

    GenmapScalar dx = min(pts[j].dx, pts[j + 1].dx)*tolSquared;

    if (d > dx)
      pts[j + 1].ifSegment = 1;
  }
}

static int mergeSegments(Mesh mesh, struct comm *c, int i, GenmapScalar tolSquared, buffer *bfr) {
  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;
  int nDim = mesh->nDim;

  sint rank = c->id;
  sint size = c->np;

  struct Point_private lastp = points[nPoints-1];
  lastp.proc = (rank + 1)%size;

  struct array arr;
  array_init(struct Point_private, &arr, 1);
  array_cat(struct Point_private, &arr, &lastp, 1);

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct Point_private, &arr, proc, 1, &cr);

  uint n = arr.n;
  assert(n == 1);
  lastp = ((struct Point_private *) arr.ptr)[0];

  if (rank > 0) {
    GenmapScalar d = sqrDiff(lastp.x[i], points->x[i]);

    GenmapScalar dx = min(lastp.dx, points->dx)*tolSquared;

    if (d > dx)
      points->ifSegment = 1;
  }

  array_free(&arr);

  // If rank > 0, send i = 0,... n-1 where points[i].ifSegment == 0 to rank - 1
  n = 0;
  if (rank > 0) {
    for (; n < nPoints && points[n].ifSegment == 0; n++)
      points[n].proc = rank - 1;
  }
  for (; n < nPoints; n++)
    points[n].proc = rank;

  sarray_transfer(struct Point_private, &mesh->elements, proc, 0, &cr);
  crystal_free(&cr);

  sarray_sort(struct Point_private, mesh->elements.ptr, mesh->elements.n, globalId, 1, bfr);

  return 0;
}

int findSegments(Mesh mesh, struct comm *c, GenmapScalar tol, int verbose, buffer *bfr) {
  int nDim = mesh->nDim;
  int nVertex = mesh->nVertex;
  GenmapScalar tolSquared = tol*tol;

  parallel_sort(struct Point_private, &mesh->elements, x[0], genmap_gs_scalar, bin_sort, 0, c);

  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  if (nDim == 3)
    sarray_sort_3(struct Point_private, points, nPoints, x[0], 3, x[1], 3, x[2], 3, bfr);
  else
    sarray_sort_2(struct Point_private, points, nPoints, x[0], 3, x[1], 3, bfr);

  //TODO: load balance

  // Initialize
  uint i;
  for (i = 0; i < nPoints; i++)
    points[i].ifSegment = 0;

  slong out[2][1], buf[2][1], in[1];
  // in[0] = nPoints;
  // comm_scan(out, &nonZeroRanks, gs_long, gs_add, in, 1, buf);
  // slong start = out[0][0];

  // if (verbose > 1 && rank == 0) {
  //   printf("\tsegments: rank=%d npts=%u start=%lld\n", rank, nPoints, start);
  //   fflush(stdout);
  // }
  // comm_barrier(&nonZeroRanks);

  // if (verbose > 1 && rank == 0) {
  //   printf("\tDone initializing\n");
  //   fflush(stdout);
  // }
  // comm_barrier(&nonZeroRanks);

  comm_ext orig;
#ifdef MPI
  MPI_Comm_dup(c->c, &orig);
#endif
  struct comm nonZeroRanks;
  
  int dim, bin, rank;
  for (dim = 0; dim < nDim; dim++) {
    nPoints = mesh->elements.n;

    bin = 1;
    if (nPoints == 0)
      bin = 0;
  
#ifdef MPI
    MPI_Comm new;
    MPI_Comm_split(orig, bin, rank, &new);
    comm_init(&nonZeroRanks, new);
    MPI_Comm_free(&new);
#else
    comm_init(&nonZeroRanks, 1);
#endif

    rank = nonZeroRanks.id;
    for (i = 0; i < nPoints; i++)
      points[i].proc = rank;

    if (bin == 1) {
      if (verbose > 0 && rank == 0) {
        printf("\tsortSegments started.\n");
      	fflush(stdout);
      }
      sortSegments(mesh, &nonZeroRanks, dim, bfr);
      comm_barrier(&nonZeroRanks);
      if (verbose > 0 && rank == 0)
        printf("\tsortSegments finished.\n");

      findLocalSegments(mesh, dim, tolSquared);
      comm_barrier(&nonZeroRanks);
      if (verbose > 0 && rank == 0)
        printf("\tfindLocalSegments finished.\n");

      mergeSegments(mesh, &nonZeroRanks, dim, tolSquared, bfr);
      comm_barrier(&nonZeroRanks);
      if (verbose > 0 && rank == 0)
        printf("\tmergeSegments finished.\n");

      if (verbose > 0) {
        nPoints = mesh->elements.n;
        points = mesh->elements.ptr;
        sint count = 0;
        for (i = 0; i < nPoints; i++)
          if (points[i].ifSegment > 0)
            count++;

        in[0] = count;
        comm_allreduce(&nonZeroRanks, gs_long, gs_add, in, 1, buf);
        if (rank == 0)
          printf("locglob: %d %lld\n", dim + 1, in[0] + 1);
      }
    }

    comm_free(&nonZeroRanks);
  }

#ifdef MPI
  MPI_Comm_free(&orig);
#endif

  return 0;
}

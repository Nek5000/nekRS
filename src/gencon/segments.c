#include <stdlib.h>
#include <math.h>

#include <gencon-impl.h>
#include <sort.h>

int sortSegments(Mesh mesh, struct comm *c, int dim) {
  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  buffer bfr;
  buffer_init(&bfr, 1024);

  uint s = 0, e;
  while (s < nPoints - 1) {
    // find the length of the segment
    for (e = s + 1; e < nPoints && points[e].ifSegment == 0; e++);

    // sort start to end based on dim
    switch (dim) {
      case 0:
        sarray_sort(struct Point_private, points + s, e - s, x[0], 3, &bfr);
        break;
      case 1:
        sarray_sort(struct Point_private, points + s, e - s, x[1], 3, &bfr);
        break;
      case 2:
        sarray_sort(struct Point_private, points + s, e - s, x[2], 3, &bfr);
        break;
      default:
        break;
    }

    s = e;
  }

  buffer_free(&bfr);

  // set sequenceId
  slong out[2][1], buf[2][1], in[1];
  in[0] = nPoints;
  comm_scan(out, c, gs_long, gs_add, in, 1, buf);
  slong start = out[0][0];

  for (e = 0; e < nPoints; e++)
    points[e].globalId = start + e;

  return 0;
}

int findLocalSegments(Mesh mesh, int i, GenmapScalar tolSquared){
  Point pts = mesh->elements.ptr;
  sint npts = mesh->elements.n;
  int nDim = mesh->nDim;

  sint j;
  for (j = 0; j < npts - 1; j++) {
#if 0
    GenmapScalar d;
    if (nDim == 3)
      d = distance3D(pts[j], pts[j + 1]);
    else
      d = distance2D(pts[j], pts[j + 1]);
#else
    GenmapScalar d = sqrDiff(pts[j].x[i], pts[j + 1].x[i]);
#endif
    GenmapScalar dx = min(pts[j].dx, pts[j + 1].dx)*tolSquared;

    if (d > dx)
      pts[j + 1].ifSegment = 1;
  }
}

int mergeSegments(Mesh mesh, struct comm *c, int i, GenmapScalar tolSquared) {
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
#if 0
    GenmapScalar d;
    if (nDim == 3)
      d = distance3D(lastp, points[0]);
    else
      d = distance2D(lastp, points[0]);
#else
    GenmapScalar d = sqrDiff(lastp.x[i], points->x[i]);
#endif

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

  buffer buf;
  buffer_init(&buf, 1024);
  sarray_sort(struct Point_private, mesh->elements.ptr, mesh->elements.n, globalId, 1, &buf);
  buffer_free(&buf);

  return 0;
}

int findSegments(Mesh mesh, struct comm *c, GenmapScalar tol, int verbose) {
  int nDim = mesh->nDim;
  int nVertex = mesh->nVertex;
  GenmapScalar tolSquared = tol*tol;

  parallel_sort(struct Point_private, &mesh->elements, x[0], genmap_gs_scalar, bin_sort, 0, c);

  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  buffer buf;
  buffer_init(&buf, 1024);

  if (nDim == 3)
    sarray_sort_3(struct Point_private, points, nPoints, x[0], 3, x[1], 3, x[2], 3, &buf);
  else
    sarray_sort_2(struct Point_private, points, nPoints, x[0], 3, x[1], 3, &buf);

  buffer_free(&buf);

  //TODO: load balance

  sint bin = 1;
  if (nPoints == 0)
    bin = 0;

  sint rank = c->id;
  sint size = c->np;
  comm_ext old = c->c;
  struct comm nonZeroRanks;
#ifdef MPI
  MPI_Comm new; MPI_Comm_split(old,bin,rank,&new);
  comm_init(&nonZeroRanks,new);
  MPI_Comm_free(&new);
#else
  comm_init(&nonZeroRanks,1);
#endif

  rank = nonZeroRanks.id;
  size = nonZeroRanks.np;

  if (bin > 0) {
    slong out[2][1], buff[2][1], in[1];
    in[0] = nPoints;
    comm_scan(out, &nonZeroRanks, gs_long, gs_add, in, 1, buff);
    slong start = out[0][0];

    if (verbose > 0) {
      printf("segments: rank=%d npts=%u start=%lld\n", rank, nPoints, start);
      fflush(stdout);
    }

    sint i;
    for (i = 0; i < nPoints; i++) {
      points[i].ifSegment = 0;
      points[i].proc = rank;
    }

    int dim;
    for (dim = 0; dim < nDim; dim++) {
      sortSegments(mesh, &nonZeroRanks, dim); // FIXME: parallel is not working
      findLocalSegments(mesh, dim, tolSquared);
      mergeSegments(mesh, &nonZeroRanks, dim, tolSquared);

      if (verbose > 0) {
        sint count = 0;
        for (i = 0; i < nPoints; i++)
          if (points[i].ifSegment > 0)
            count++;

        in[0] = count;
        comm_allreduce(&nonZeroRanks, gs_long, gs_add, in, 1, buff);
        if (rank == 0)
          printf("locglob: %d %lld\n", dim+1, in[0]+1);
      }
    }
  }

  comm_free(&nonZeroRanks);

  return 0;
}

#include "gencon-impl.h"

int setGlobalID(Mesh mesh, struct comm *c) {
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

    slong out[2][1], buf[2][1], in[1];
    in[0] = count;
    comm_scan(out, &nonZeroRanks, gs_long, gs_add, in, 1, buf);
    slong start = out[0][0];

    assert(points[0].ifSegment == 1);
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

int sendBack(Mesh mesh, struct comm *c, buffer *bfr) {
  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct Point_private, &mesh->elements, origin, 0, &cr);
  crystal_free(&cr);

  sarray_sort(struct Point_private, mesh->elements.ptr, mesh->elements.n,
              sequenceId, 1, bfr);

  return 0;
}

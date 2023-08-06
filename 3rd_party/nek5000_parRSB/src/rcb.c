#include "parrsb-impl.h"
#include "sort.h"

static void get_axis_len(double *length, size_t unit_size, char *elems,
                         uint nel, int ndim, struct comm *c) {
  double min[3] = {DBL_MAX, DBL_MAX, DBL_MAX},
         max[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};

  struct rcb_element *ei;
  for (uint i = 0; i < nel; i++) {
    ei = (struct rcb_element *)(elems + i * unit_size);
    if (ei->coord[0] < min[0])
      min[0] = ei->coord[0];
    if (ei->coord[0] > max[0])
      max[0] = ei->coord[0];

    if (ei->coord[1] < min[1])
      min[1] = ei->coord[1];
    if (ei->coord[1] > max[1])
      max[1] = ei->coord[1];
  }

  if (ndim == 3) {
    for (uint i = 0; i < nel; i++) {
      ei = (struct rcb_element *)(elems + i * unit_size);
      if (ei->coord[2] < min[2])
        min[2] = ei->coord[2];
      if (ei->coord[2] > max[2])
        max[2] = ei->coord[2];
    }
  }

  if (c != NULL) {
    double wrk[3];
    comm_allreduce(c, gs_double, gs_min, min, 3, wrk);
    comm_allreduce(c, gs_double, gs_max, max, 3, wrk);
  }

  for (uint i = 0; i < ndim; i++)
    length[i] = max[i] - min[i];
}

void rcb_local(struct array *a, size_t unit_size, uint start, uint end,
               int ndim, buffer *buf) {
  sint size = end - start;
  if (size <= 1)
    return;

  double length[3];
  char *st = (char *)a->ptr + unit_size * start;
  get_axis_len(length, unit_size, st, size, ndim, NULL);

  int axis = 0;
  if (length[1] > length[0])
    axis = 1;
  if (ndim == 3)
    if (length[2] > length[axis])
      axis = 2;

  if (unit_size == sizeof(struct rcb_element)) {
    switch (axis) {
    case 0:
      sarray_sort(struct rcb_element, st, size, coord[0], 3, buf);
      break;
    case 1:
      sarray_sort(struct rcb_element, st, size, coord[1], 3, buf);
      break;
    case 2:
      sarray_sort(struct rcb_element, st, size, coord[2], 3, buf);
      break;
    default:
      break;
    }
  } else if (unit_size == sizeof(struct rsb_element)) {
    switch (axis) {
    case 0:
      sarray_sort(struct rsb_element, st, size, coord[0], 3, buf);
      break;
    case 1:
      sarray_sort(struct rsb_element, st, size, coord[1], 3, buf);
      break;
    case 2:
      sarray_sort(struct rsb_element, st, size, coord[2], 3, buf);
      break;
    default:
      break;
    }
  }

  uint mid = (start + end) / 2;
  rcb_local(a, unit_size, start, mid, ndim, buf);
  rcb_local(a, unit_size, mid, end, ndim, buf);
}

static int rcb_level(struct array *a, size_t unit_size, int ndim,
                     struct comm *c, buffer *bfr) {
  if (c->np == 1)
    return 0;

  double length[3];
  get_axis_len(length, unit_size, (char *)a->ptr, a->n, ndim, c);

  int axis = 0, d;
  for (d = 1; d < ndim; d++)
    if (length[d] > length[axis])
      axis = d;

  if (unit_size == sizeof(struct rcb_element)) {
    switch (axis) {
    case 0:
      parallel_sort(struct rcb_element, a, coord[0], gs_double, 0, 1, c, bfr);
      break;
    case 1:
      parallel_sort(struct rcb_element, a, coord[1], gs_double, 0, 1, c, bfr);
      break;
    case 2:
      parallel_sort(struct rcb_element, a, coord[2], gs_double, 0, 1, c, bfr);
      break;
    default:
      break;
    }
  } else if (unit_size == sizeof(struct rsb_element)) {
    switch (axis) {
    case 0:
      parallel_sort(struct rsb_element, a, coord[0], gs_double, 0, 1, c, bfr);
      break;
    case 1:
      parallel_sort(struct rsb_element, a, coord[1], gs_double, 0, 1, c, bfr);
      break;
    case 2:
      parallel_sort(struct rsb_element, a, coord[2], gs_double, 0, 1, c, bfr);
      break;
    default:
      break;
    }
  }

  return 0;
}

int rcb(struct array *elements, size_t unit_size, int ndim, struct comm *ci,
        buffer *bfr) {
  struct comm c, t;
  comm_dup(&c, ci);

  int size = c.np;
  int rank = c.id;

  while (size > 1) {
    rcb_level(elements, unit_size, ndim, &c, bfr);

    int bin = 1;
    if (rank < (size + 1) / 2)
      bin = 0;

    comm_split(&c, bin, rank, &t);
    comm_free(&c);
    comm_dup(&c, &t);
    comm_free(&t);

    size = c.np, rank = c.id;
  }
  comm_free(&c);

  rcb_local(elements, unit_size, 0, elements->n, ndim, bfr);

  return 0;
}

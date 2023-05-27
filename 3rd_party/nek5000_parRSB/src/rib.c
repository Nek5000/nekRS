#include "parrsb-impl.h"
#include "sort.h"

extern int power_serial(double *y, int N, double *A, int verbose);

static void get_rib_axis(char *elems, uint nel, size_t unit_size, int ndim,
                         struct comm *c) {
  double avg[3];
  avg[0] = avg[1] = avg[2] = 0.0;
  struct rcb_element *elem;
  uint i;
  for (i = 0; i < nel; i++) {
    elem = (struct rcb_element *)(elems + i * unit_size);
    avg[0] += elem->coord[0];
    avg[1] += elem->coord[1];
    avg[2] += elem->coord[2];
  }
  slong nelg = nel;

  double buf[9];
  if (c != NULL) {
    comm_allreduce(c, gs_double, gs_add, avg, 3, buf);
    comm_allreduce(c, gs_long, gs_add, &nelg, 1, buf);
  }

  avg[0] /= nelg;
  avg[1] /= nelg;
  avg[2] /= nelg;

  double I[3][3];
  for (i = 0; i < 3; i++)
    I[i][0] = I[i][1] = I[i][2] = 0.0;

  double x, y, z;
  for (i = 0; i < nel; i++) {
    elem = (struct rcb_element *)(elems + i * unit_size);
    x = elem->coord[0] - avg[0];
    y = elem->coord[1] - avg[1];
    z = elem->coord[2] - avg[2];
    I[0][0] += x * x, I[0][1] += x * y, I[0][2] += x * z;
    I[1][0] += y * x, I[1][1] += y * y, I[1][2] += y * z;
    I[2][0] += z * x, I[2][1] += z * y, I[2][2] += z * z;
  }

  if (c != NULL)
    comm_allreduce(c, gs_double, gs_add, I, 9, buf);

  double ev[3];                           // ev[2] = 0 if 2D
  power_serial(ev, ndim, (double *)I, 0); // FIXME: 2D does not work

  for (i = 0; i < nel; i++) {
    elem = (struct rcb_element *)(elems + i * unit_size);
    x = elem->coord[0] - avg[0];
    y = elem->coord[1] - avg[1];
    z = elem->coord[2] - avg[2];
    elem->fiedler = x * ev[0] + y * ev[1] + z * ev[2];
  }
}

void rib_local(struct array *a, size_t unit_size, uint start, uint end,
               int ndim, buffer *buf) {
  sint size = end - start;
  if (size <= 1)
    return;

  char *st = (char *)a->ptr + unit_size * start;
  get_rib_axis(st, size, unit_size, ndim, NULL);

  if (unit_size == sizeof(struct rcb_element))
    sarray_sort(struct rcb_element, st, size, fiedler, 3, buf);
  else if (unit_size == sizeof(struct rsb_element))
    sarray_sort(struct rsb_element, st, size, fiedler, 3, buf);

  uint mid = (start + end) / 2;
  rib_local(a, unit_size, start, mid, ndim, buf);
  rib_local(a, unit_size, mid, end, ndim, buf);
}

static int rib_level(struct array *a, size_t unit_size, int ndim,
                     struct comm *c, buffer *bfr) {
  if (c->np == 1)
    return 0;

  get_rib_axis((char *)a->ptr, a->n, unit_size, ndim, c);

  if (unit_size == sizeof(struct rcb_element))
    parallel_sort(struct rcb_element, a, fiedler, gs_double, 0, 1, c, bfr);
  else if (unit_size == sizeof(struct rsb_element))
    parallel_sort(struct rsb_element, a, fiedler, gs_double, 0, 1, c, bfr);

  return 0;
}

int rib(struct array *elements, size_t unit_size, int ndim, struct comm *ci,
        buffer *bfr) {
  struct comm c;
  comm_dup(&c, ci);

  int size = c.np;
  int rank = c.id;

  while (size > 1) {
    rib_level(elements, unit_size, ndim, &c, bfr);

    int p = (size + 1) / 2;
    int bin = (rank >= p);

    MPI_Comm comm_rib;
    MPI_Comm_split(c.c, bin, rank, &comm_rib);
    comm_free(&c);
    comm_init(&c, comm_rib);
    MPI_Comm_free(&comm_rib);

    size = c.np;
    rank = c.id;
  }
  comm_free(&c);

  rib_local(elements, unit_size, 0, elements->n, ndim, bfr);

  return 0;
}

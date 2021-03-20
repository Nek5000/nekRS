#include <float.h>
#include <genmap-impl.h>
#include <sort.h>

void get_rib_axis_local(void *elems, uint nel, int ndim) {
  // TODO: Get rid of this
  size_t unit_size;
  unsigned char *type = elems;
  if (*type == GENMAP_RCB_ELEMENT) {
    unit_size = sizeof(struct rcb_element);
  } else if (*type == GENMAP_RSB_ELEMENT) {
    unit_size = sizeof(struct rsb_element);
  }

  struct rcb_element *elem;
  uint i;

  double avg[3];
  avg[0] = avg[1] = avg[2] = 0.0;
  for (i = 0; i < nel; i++) {
    elem = (struct rcb_element *)((char *)elems + i * unit_size);
    avg[0] += elem->coord[0];
    avg[1] += elem->coord[1];
    avg[2] += elem->coord[2];
  }

  avg[0] /= nel;
  avg[1] /= nel;
  avg[2] /= nel;

  double I[3][3];
  for (i = 0; i < 3; i++)
    I[i][0] = I[i][1] = I[i][2] = 0.0;

  double x, y, z;
  for (i = 0; i < nel; i++) {
    elem = (struct rcb_element *)((char *)elems + i * unit_size);
    x = elem->coord[0] - avg[0];
    y = elem->coord[1] - avg[1];
    z = elem->coord[2] - avg[2];
    I[0][0] += x * x, I[0][1] += x * y, I[0][2] += x * z;
    I[1][0] += y * x, I[1][1] += y * y, I[1][2] += y * z;
    I[2][0] += z * x, I[2][1] += z * y, I[2][2] += z * z;
  }

  double ev[3];                           // ev[2] = 0 if 2D
  genmap_power(ev, ndim, (double *)I, 0); // FIXME: 2D does not work

  for (i = 0; i < nel; i++) {
    elem = (struct rcb_element *)((char *)elems + i * unit_size);
    x = elem->coord[0] - avg[0];
    y = elem->coord[1] - avg[1];
    z = elem->coord[2] - avg[2];
    elem->fiedler = x * ev[0] + y * ev[1] + z * ev[2];
  }
}

void get_rib_axis(struct array *a, struct comm *c, int ndim) {
  void *elems = a->ptr;

  // TODO: Get rid of this
  size_t unit_size;
  unsigned char *type = elems;
  if (*type == GENMAP_RCB_ELEMENT) {
    unit_size = sizeof(struct rcb_element);
  } else if (*type == GENMAP_RSB_ELEMENT) {
    unit_size = sizeof(struct rsb_element);
  }

  uint i;
  uint nel = a->n;
  struct rcb_element *elem;

  double avg[4];
  avg[0] = avg[1] = avg[2] = 0.0;
  for (i = 0; i < nel; i++) {
    elem = (struct rcb_element *)((char *)elems + i * unit_size);
    avg[0] += elem->coord[0];
    avg[1] += elem->coord[1];
    avg[2] += elem->coord[2];
  }
  avg[3] = nel;

  double buf[9];
  comm_allreduce(c, gs_double, gs_add, avg, 4, buf);

  avg[0] /= avg[3];
  avg[1] /= avg[3];
  avg[2] /= avg[3];

  double I[3][3];
  for (i = 0; i < 3; i++)
    I[i][0] = I[i][1] = I[i][2] = 0.0;

  double x, y, z;
  for (i = 0; i < nel; i++) {
    elem = (struct rcb_element *)((char *)elems + i * unit_size);
    x = elem->coord[0] - avg[0];
    y = elem->coord[1] - avg[1];
    z = elem->coord[2] - avg[2];
    I[0][0] += x * x, I[0][1] += x * y, I[0][2] += x * z;
    I[1][0] += y * x, I[1][1] += y * y, I[1][2] += y * z;
    I[2][0] += z * x, I[2][1] += z * y, I[2][2] += z * z;
  }

  comm_allreduce(c, gs_double, gs_add, I, 9, buf);

  double ev[3];                           // ev[2] = 0 if 2D
  genmap_power(ev, ndim, (double *)I, 0); // FIXME: 2D does not work

  for (i = 0; i < nel; i++) {
    elem = (struct rcb_element *)((char *)elems + i * unit_size);
    x = elem->coord[0] - avg[0];
    y = elem->coord[1] - avg[1];
    z = elem->coord[2] - avg[2];
    elem->fiedler = x * ev[0] + y * ev[1] + z * ev[2];
  }
}

void rib_local(struct array *a, uint start, uint end, int ndim, buffer *buf) {
  sint size = end - start;
  assert(size >= 0);

  if (size <= 2)
    return;

  // TODO: Get rid of this
  size_t unit_size;
  unsigned char *type = a->ptr;
  if (*type == GENMAP_RCB_ELEMENT) {
    unit_size = sizeof(struct rcb_element);
  } else if (*type == GENMAP_RSB_ELEMENT) {
    unit_size = sizeof(struct rsb_element);
  }

  void *st = (void *)a->ptr + unit_size * start;
  get_rib_axis_local(st, size, ndim);

  if (*type == GENMAP_RCB_ELEMENT) {
    sarray_sort(struct rcb_element, st, size, fiedler, 3, buf);
  } else if (*type == GENMAP_RSB_ELEMENT) {
    sarray_sort(struct rsb_element, st, size, fiedler, 3, buf);
  }

  uint mid = (start + end) / 2;
  rib_local(a, start, mid, ndim, buf);
  rib_local(a, mid, end, ndim, buf);
}

int rib_level(struct comm *c, struct array *a, int ndim, buffer *bfr) {
  if (c->np == 1)
    return 0;

  get_rib_axis(a, c, ndim);

  unsigned char *type = a->ptr;
  if (*type == GENMAP_RCB_ELEMENT) {
    parallel_sort(struct rcb_element, a, fiedler, gs_double, 0, 1, c, bfr);
  } else if (*type == GENMAP_RSB_ELEMENT) {
    parallel_sort(struct rsb_element, a, fiedler, gs_double, 0, 1, c, bfr);
  }

  return 0;
}

int rib(struct comm *ci, struct array *elements, int ndim, buffer *bfr) {
  struct comm c;
  comm_dup(&c, ci);

  int size = c.np;
  int rank = c.id;

  while (size > 1) {
    rib_level(&c, elements, ndim, bfr);

    int p = (size + 1) / 2;
    int bin = (rank >= p);

    comm_ext comm_rib;
#ifdef MPI
    MPI_Comm_split(c.c, bin, rank, &comm_rib);
#else
    comm_rib = 1;
#endif

    comm_free(&c);
    comm_init(&c, comm_rib);

#ifdef MPI
    MPI_Comm_free(&comm_rib);
#endif

    size = c.np;
    rank = c.id;
  }

  rib_local(elements, 0, elements->n, ndim, bfr);

  comm_free(&c);

  return 0;
}

int genmap_rib(genmap_handle h) {
  struct comm *lc = h->local;

  int ndim = (h->nv == 8) ? 3 : 2;

  rib(lc, h->elements, ndim, &h->buf);

  struct rcb_element *eptr = h->elements->ptr;
  int e;
  for (e = 0; e < h->elements->n; e++)
    eptr[e].seq = e;

  return 0;
}

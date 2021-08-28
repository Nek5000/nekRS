#include <genmap-impl.h>

/* Load balance input data */
size_t genmap_load_balance(struct array *eList, uint nel, int nv, double *coord,
                           long long *vtx, struct crystal *cr, buffer *bfr) {
  slong in = nel;
  slong out[2][1], buf[2][1];
  comm_scan(out, &cr->comm, gs_long, gs_add, &in, 1, buf);
  slong start = out[0][0];
  slong nelg = out[1][0];

  int size = cr->comm.np;
  uint nstar = nelg / size;
  uint nrem = nelg - nstar * size;
  slong lower = (nstar + 1) * nrem;

  size_t unit_size;
  struct rcb_element *element = NULL;

  if (vtx == NULL) { // RCB
    unit_size = sizeof(struct rcb_element);
    element = calloc(1, sizeof(struct rcb_element));
    element->type = GENMAP_RCB_ELEMENT;
  } else {
    unit_size = sizeof(struct rsb_element);
    element = calloc(1, sizeof(struct rsb_element));
    element->type = GENMAP_RSB_ELEMENT;
  }

  element->origin = cr->comm.id;

  array_init_(eList, nel, unit_size, __FILE__, __LINE__);

  int ndim = (nv == 8) ? 3 : 2;
  int e, n, v;
  for (e = 0; e < nel; ++e) {
    slong eg = element->globalId = start + e + 1;
    if (eg <= lower)
      element->proc = (eg - 1) / (nstar + 1);
    else if (nstar != 0)
      element->proc = (eg - 1 - lower) / nstar + nrem;

    element->coord[0] = element->coord[1] = element->coord[2] = 0.0;
    for (v = 0; v < nv; v++)
      for (n = 0; n < ndim; n++)
        element->coord[n] += coord[e * ndim * nv + v * ndim + n];
    for (n = 0; n < ndim; n++)
      element->coord[n] /= nv;

    array_cat_(unit_size, eList, element, 1, __FILE__, __LINE__);
  }
  assert(eList->n == nel);

  if (vtx != NULL) { // RSB
    struct rsb_element *elements = eList->ptr;
    for (e = 0; e < nel; e++) {
      for (v = 0; v < nv; v++)
        elements[e].vertices[v] = vtx[e * nv + v];
    }
  }

  sarray_transfer_(eList, unit_size, offsetof(struct rcb_element, proc), 1, cr);
  nel = eList->n;
  if (vtx != NULL) // RSB
    sarray_sort(struct rsb_element, eList->ptr, nel, globalId, 1, bfr);
  else
    sarray_sort(struct rcb_element, eList->ptr, nel, globalId, 1, bfr);

  free(element);
  return unit_size;
}

void genmap_restore_original(int *part, int *seq, struct crystal *cr,
                             struct array *eList, buffer *bfr) {
  struct rcb_element *element = eList->ptr;
  size_t unit_size;
  if (element->type == GENMAP_RSB_ELEMENT) // RSB
    unit_size = sizeof(struct rsb_element);
  else
    unit_size = sizeof(struct rcb_element);

  sarray_transfer_(eList, unit_size, offsetof(struct rcb_element, origin), 1,
                   cr);
  uint nel = eList->n;

  if (element->type == GENMAP_RSB_ELEMENT) // RSB
    sarray_sort(struct rsb_element, eList->ptr, nel, globalId, 1, bfr);
  else
    sarray_sort(struct rcb_element, eList->ptr, nel, globalId, 1, bfr);

  int e;
  for (e = 0; e < nel; e++) {
    element = (struct rcb_element *)(eList->ptr + e * unit_size);
    part[e] = element->origin; // element[e].origin;
  }

  if (seq != NULL)
    for (e = 0; e < nel; e++) {
      element = eList->ptr + e * unit_size;
      seq[e] = element->seq; // element[e].seq;
    }
}

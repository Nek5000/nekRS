#include <sort-impl.h>

/* assumes array is locally sorted */
int set_bin(uint **proc_, struct sort *s, uint field, struct comm *c) {
  struct array *a = s->a;
  gs_dom t = s->t[field];
  uint offset = s->offset[field];

  uint size = a->n;
  GenmapCalloc(size, proc_);
  uint *proc = *proc_;

  double extrema[2];
  get_extrema((void *)extrema, s, field, c);
  double range = extrema[1] - extrema[0];

  if (size == 0)
    return 0;

  sint np = c->np;
  uint id = 0;
  uint index = 0;
  do {
    double end = extrema[0] + (range / np) * (id + 1);
    while (index < size) {
      double val = get_scalar(a, index, offset, s->unit_size, t);
      if (val <= end)
        proc[index++] = id;
      else
        break;
    }
    id++;
  } while (id < np && index < size);
  for (; index < size; index++)
    proc[index] = np - 1;
}

int parallel_bin_sort(struct sort *s, struct comm *c) {
  /* Local sort */
  sort_local(s);

  /* Set destination bin */
  uint *proc;
  set_bin(&proc, s, 0, c);

  /* Transfer to destination processor */
  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer_ext_(s->a, s->unit_size, proc, sizeof(uint), &cr);
  crystal_free(&cr);

  GenmapFree(proc);

  /* Locally sort again */
  sort_local(s);

  return 0;
}

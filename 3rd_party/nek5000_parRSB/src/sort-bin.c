#include "sort-impl.h"

static uint *set_proc_from_val(struct sort *s, uint field,
                               const struct comm *c) {
  struct array *a = s->a;
  gs_dom t = s->t[field];
  uint offset = s->offset[field];

  double extrema[2];
  get_extrema((void *)extrema, s, field, c);
  double range = extrema[1] - extrema[0];

  uint size = a->n;
  if (size == 0) return NULL;
  uint *proc = tcalloc(uint, size);

  uint np = c->np;
  assert(np > 0);
  uint id = 0, index = 0;
  do {
    double end = extrema[0] + (range / np) * (id + 1);
    while (index < size) {
      double val = get_scalar(a, index, offset, s->unit_size, t);
      if (val <= end)
        proc[index] = id, index++;
      else
        break;
    }
    id++;
  } while (id < np && index < size);
  for (; index < size; index++) proc[index] = np - 1;

  return proc;
}

void parallel_bin_sort(struct sort *s, const struct comm *c) {
  // Locally sort the array first.
  sort_local(s);

  // Set destination bin based on the field value.
  uint *proc = set_proc_from_val(s, 0, c);

  // Transfer the array in chunks.
  sarray_transfer_chunk(s->a, s->unit_size, proc, c);
  free(proc);

  // Locally sort again to make sure that we have both globally and locally
  // sorted array.
  sort_local(s);
}

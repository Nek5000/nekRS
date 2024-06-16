#include "sort-impl.h"
#include <math.h>

struct hypercube {
  struct sort *data;
  int nprobes;
  double *probes;
  ulong *probe_cnt;
};

static void init_probes(struct hypercube *data, const struct comm *c) {
  // Allocate space for probes and counts.
  int nprobes = data->nprobes = 3;
  if (!data->probes) data->probes = tcalloc(double, nprobes);
  if (!data->probe_cnt) data->probe_cnt = tcalloc(ulong, nprobes);

  double extrema[2];
  get_extrema((void *)extrema, data->data, 0, c);
  double range = extrema[1] - extrema[0];
  double delta = range / (nprobes - 1);

  data->probes[0] = extrema[0];
  data->probes[1] = extrema[0] + delta;
  data->probes[2] = extrema[1];
}

static void update_probe_counts(struct hypercube *data, const struct comm *c) {
  struct sort *input = data->data;
  uint offset = input->offset[0];
  gs_dom t = input->t[0];

  uint nprobes = data->nprobes;
  for (uint i = 0; i < nprobes; i++) data->probe_cnt[i] = 0;

  struct array *a = input->a;
  for (uint e = 0; e < a->n; e++) {
    double val = get_scalar(a, e, offset, input->unit_size, t);
    for (uint i = 0; i < nprobes; i++) {
      if (val < data->probes[i]) data->probe_cnt[i]++;
    }
  }

  slong wrk[6];
  comm_allreduce(c, gs_long, gs_add, data->probe_cnt, nprobes, wrk);
}

static void update_probes(slong nelem, double *probes, ulong *probe_cnt,
                          uint threshold) {
  assert(nelem >= 0);
  slong expected = nelem / 2;
  if (llabs(expected - (slong)probe_cnt[1]) < threshold) return;

  if (probe_cnt[1] < (ulong)expected)
    probes[0] = probes[1];
  else
    probes[2] = probes[1];

  probes[1] = probes[0] + (probes[2] - probes[0]) / 2.0;
}

static void transfer_elem(const struct hypercube *data, const struct comm *c) {
  struct sort *input = data->data;
  uint usize = input->unit_size;
  uint offset = input->offset[0];
  gs_dom t = input->t[0];
  struct array *a = input->a;

  uint size = a->n, lown = 0, uppern = 0;
  for (uint e = 0; e < size; e++) {
    double val = get_scalar(a, e, offset, usize, t);
    if (val < data->probes[1])
      lown++;
    else
      uppern++;
  }

  slong out[2][2], in[2] = {lown, uppern}, wrk[2][2];
  comm_scan(out, c, gs_long, gs_add, in, 2, wrk);
  slong lstart = out[0][0], ustart = out[0][1];
  slong lelem = out[1][0], uelem = out[1][1];

  uint np = c->np, lnp = np / 2;
  uint *proc1 = set_proc_from_idx(lnp, lstart, lown, lelem);
  uint *proc2 = set_proc_from_idx(np - lnp, ustart, uppern, uelem);
  proc1 = trealloc(uint, proc1, size);
  for (uint e = lown; e < size; e++) proc1[e] = proc2[e - lown] + lnp;

  sarray_transfer_chunk(a, usize, proc1, c);
  free(proc1), free(proc2);
}

// TODO: Get rid of this recursive implementation.
static void parallel_hypercube_sort_aux(struct hypercube *data,
                                        const struct comm *c) {
  struct sort *input = data->data;
  struct array *a = input->a;

  // FIXME: Replace comm_scan() by comm_allreduce().
  slong out[2][1], buf[2][1], in = a->n;
  comm_scan(out, c, gs_long, gs_add, &in, 1, buf);
  slong nelem = out[1][0];

  uint threshold = nelem / (10 * c->np);
  if (threshold < 2) threshold = 2;

  sort_local(data->data);

  if (c->np == 1) return;

  init_probes(data, c);
  update_probe_counts(data, c);
  int max_iter = log2((data->probes[2] - data->probes[0]) / 1e-12);
  int iter = 0;
  while (llabs(nelem / 2 - (slong)data->probe_cnt[1]) > threshold &&
         iter++ < max_iter) {
    update_probes(nelem, data->probes, data->probe_cnt, threshold);
    update_probe_counts(data, c);
  }

  transfer_elem(data, c);

  // split the communicator
  struct comm nc;
  sint lower = (c->id < c->np / 2);
  comm_split(c, lower, c->id, &nc);

  // TODO: Keep load balancing after each split
  parallel_hypercube_sort_aux(data, &nc);

  comm_free(&nc);
}

void parallel_hypercube_sort(struct sort *sd, const struct comm *c) {
  struct comm dup;
  comm_dup(&dup, c);

  struct hypercube hdata = {.data = sd, .probes = NULL, .probe_cnt = NULL};
  parallel_hypercube_sort_aux(&hdata, &dup);
  free(hdata.probes), free(hdata.probe_cnt);

  comm_free(&dup);
}

#include <float.h>
#include <sort-impl.h>

int parallel_sort_private(struct sort *data, struct comm *c) {
  struct comm dup;
  comm_dup(&dup, c);

  int balance = data->balance;
  sort_algo algo = data->algo;

  struct array *a = data->a;
  size_t usize = data->unit_size;

  struct hypercube hdata;

  switch (algo) {
  case bin_sort:
    parallel_bin_sort(data, c);
    break;
  case hypercube_sort:
    hdata.data = data;
    hdata.probes = NULL;
    hdata.probe_cnt = NULL;
    parallel_hypercube_sort(&hdata, &dup);
    free(hdata.probes);
    free(hdata.probe_cnt);
    break;
  default:
    break;
  }

  if (balance) {
    struct crystal cr;
    crystal_init(&cr, c);
    load_balance(a, usize, c, &cr);
    crystal_free(&cr);
    sort_local(data);
  }

  comm_free(&dup);

  return 0;
}

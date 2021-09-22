#ifndef _SORT_IMPL_H_
#define _SORT_IMPL_H_

#include "sort.h"
#ifndef _SORT_H_
#error "sort.h not included"
#endif

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

/*
 * local sort
 */

int sort_local(struct sort *s);

/*
 * parallel sort
 */

struct hypercube {
  struct sort *data;
  int nprobes;
  double *probes;
  ulong *probe_cnt;
};

int parallel_bin_sort(struct sort *s, struct comm *c);
int parallel_hypercube_sort(struct hypercube *h, struct comm *c);

/*
 * Auxiliary functions
 */
double get_scalar(struct array *a, uint i, uint off, uint usize, gs_dom type);
void get_extrema(void *extrema, struct sort *s, uint field, struct comm *c);

int set_dest(uint *proc, uint size, sint np, slong start, slong nelem);
int load_balance(struct array *a, size_t size, struct comm *c,
                 struct crystal *cr);
#endif

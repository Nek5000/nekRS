#ifndef _SORT_IMPL_H_
#define _SORT_IMPL_H_

#include "sort.h"
#ifndef _SORT_H_
#error "sort.h not included"
#endif

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

double get_scalar(struct array *a, uint i, uint off, uint usize, gs_dom type);
void get_extrema(void *extrema, struct sort *s, uint field, struct comm *c);
int set_dest(uint *proc, uint np, ulong start, uint size, ulong nelem);
int load_balance(struct array *a, size_t size, struct comm *c,
                 struct crystal *cr);
int sort_local(struct sort *s);

#endif

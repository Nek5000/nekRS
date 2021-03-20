#ifndef _SORT_H_
#define _SORT_H_

#include <genmap-impl.h>

typedef enum { bin_sort = 0, hypercube_sort = 1 } sort_algo;

struct sort {
  int nfields;
  gs_dom t[3];
  uint offset[3];

  struct array *a;
  size_t unit_size, align;

  int balance;
  sort_algo algo;

  buffer *buf;
};
//
// parallel_bin_sort
//
int parallel_sort_private(struct sort *s, struct comm *c);
int parallel_bin_sort(struct sort *s, struct comm *c);
//
// parallel_hypercube_sort
//
struct hypercube {
  struct sort *data;
  int nprobes;
  double *probes;
  ulong *probe_cnt;
};
int parallel_hypercube_sort(struct hypercube *h, struct comm *c);
//
// Uniform parallel sort
//
#define parallel_sort(T, A, field, type, method, loadbalance, c, bufp)         \
  do {                                                                         \
    struct sort sd;                                                            \
    sd.unit_size = sizeof(T);                                                  \
    sd.align = ALIGNOF(T);                                                     \
    sd.nfields = 1;                                                            \
    sd.t[0] = type;                                                            \
    sd.offset[0] = offsetof(T, field);                                         \
    sd.a = A;                                                                  \
    sd.algo = method;                                                          \
    sd.balance = loadbalance;                                                  \
    sd.buf = bufp;                                                             \
    parallel_sort_private(&sd, c);                                             \
  } while (0)

#define parallel_sort_2(T, A, f1, t1, f2, t2, method, loadbalance, c, bufp)    \
  do {                                                                         \
    struct sort sd;                                                            \
    sd.unit_size = sizeof(T);                                                  \
    sd.align = ALIGNOF(T);                                                     \
    sd.nfields = 2;                                                            \
    sd.t[0] = t1;                                                              \
    sd.offset[0] = offsetof(T, f1);                                            \
    sd.t[1] = t2;                                                              \
    sd.offset[1] = offsetof(T, f2);                                            \
    sd.a = A;                                                                  \
    sd.algo = method;                                                          \
    sd.balance = loadbalance;                                                  \
    sd.buf = bufp;                                                             \
    parallel_sort_private(&sd, c);                                             \
  } while (0)

#endif

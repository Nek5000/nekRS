#ifndef _PARRSB_SORT_H_
#define _PARRSB_SORT_H_

#include <gslib.h>

struct sort {
  int balance, algo;

  int nfields;
  gs_dom t[3];
  uint offset[3];

  struct array *a;
  size_t unit_size, align;

  buffer *buf;
};

int sort_local(struct sort *s);
int parallel_sort_private(struct sort *s, const struct comm *c);

// Uniform parallel sort
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

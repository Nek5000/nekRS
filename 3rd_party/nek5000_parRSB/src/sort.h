#ifndef _PARRSB_SORT_H_
#define _PARRSB_SORT_H_

#include <gslib.h>
#include <stdarg.h>

void parallel_sort_(struct array *arr, size_t usize, size_t align,
                    unsigned algo, unsigned balance, const struct comm *c,
                    buffer *bfr, unsigned nfields, ...);

#define parallel_sort(T, A, field, type, algo, balance, c, bfr)                \
  parallel_sort_(A, sizeof(T), ALIGNOF(T), algo, balance, c, bfr, 1, type,     \
                 offsetof(T, field))

#define parallel_sort_2(T, A, f1, t1, f2, t2, algo, balance, c, bfr)           \
  parallel_sort_(A, sizeof(T), ALIGNOF(T), algo, balance, c, bfr, 2, t1,       \
                 offsetof(T, f1), t2, offsetof(T, f2))

#endif

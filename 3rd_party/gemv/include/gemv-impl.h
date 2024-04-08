#if !defined(__GEMV_IMPL_H__)
#define __GEMV_IMPL_H__

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "gemv-defs.h"
#include "gemv.h"

// Dynamic memory allocation functions.
#define gemv_malloc(n) malloc((n))
#define gemv_calloc(T, n) ((T *)calloc(n, sizeof(T)))
#define gemv_realloc(ptr, T, n) (T *)realloc(ptr, (n) * sizeof(T))

// Dynamic memory deallocation function.
GEMV_INTERN void gemv_free_(void **p);
#define gemv_free(p) gemv_free_((void **)p)

GEMV_INTERN void gemv_log(const gemv_verbose_t verbose, const char *fmt, ...);

struct gemv_t {
  int device, backend;
  gemv_precision_t precision;
  unsigned m, n;
  double *A;
};

struct gemv_backend_t {
  char name[32];
  void (*init)(struct gemv_backend_t *backend, const struct gemv_t *gemv);
  void (*malloc)(void **ptr, const size_t size);
  void (*free)(void *ptr);
  void (*copy)(void *dest, const void *src, size_t count,
               gemv_direction_t direction);
  void (*run)(void *, const void *, const struct gemv_t *);
  void (*finalize)(void);
};

GEMV_INTERN void gemv_set_backend_impl(struct gemv_t *gemv,
                                       const char *backend);

GEMV_INTERN void gemv_set_verbose_impl(const gemv_verbose_t level);

GEMV_INTERN void gemv_backend_register(
    const char *name,
    void (*init)(struct gemv_backend_t *backend, const struct gemv_t *gemv));

GEMV_INTERN void gemv_backend_init(const struct gemv_t *gemv);

GEMV_INTERN void gemv_backend_malloc(void **ptr, size_t size);

GEMV_INTERN void gemv_backend_free(void **ptr);

GEMV_INTERN void gemv_backend_run(void *y, const void *x,
                                  const struct gemv_t *gemv);

GEMV_INTERN void gemv_backend_copy(void *dest, const void *src, size_t count,
                                   gemv_direction_t direction);

GEMV_INTERN void gemv_backend_finalize(void);

GEMV_INTERN void gemv_backend_deregister(void);

#define GEMV_BACKEND(name) GEMV_INTERN void gemv_register_##name(void);
#include "backends/gemv-backend-list.h"
#undef GEMV_BACKEND

#endif // __GEMV_IMPL_H__

#if !defined(__GEMV_H__)
#define __GEMV_H__

#define GEMV_VISIBILITY(mode) __attribute__((visibility(#mode)))

#if defined(__cplusplus)
#define GEMV_EXTERN extern "C" GEMV_VISIBILITY(default)
#else
#define GEMV_EXTERN extern GEMV_VISIBILITY(default)
#endif

#if defined(__cplusplus)
#define GEMV_INTERN extern "C" GEMV_VISIBILITY(hidden)
#else
#define GEMV_INTERN extern GEMV_VISIBILITY(hidden)
#endif

#include <stddef.h>

struct gemv_t;
GEMV_EXTERN struct gemv_t *gemv_init(int *argc, char ***argv);

typedef enum {
  GEMV_MUTE = 0,
  GEMV_INFO = 1,
  GEMV_WARNING = 2,
  GEMV_ERROR = 3
} gemv_verbose_t;

GEMV_EXTERN void gemv_set_verbose(gemv_verbose_t verbose);

GEMV_EXTERN void gemv_set_device(struct gemv_t *gemv, int device);

GEMV_EXTERN void gemv_set_backend(struct gemv_t *gemv, const char *backend);

typedef enum { GEMV_FP64 = 0, GEMV_FP32 = 1 } gemv_precision_t;
GEMV_EXTERN void gemv_set_precision(struct gemv_t *gemv,
                                    gemv_precision_t precision);

GEMV_EXTERN void gemv_set_matrix(struct gemv_t *gemv, unsigned m, unsigned n,
                                 const double *A);

GEMV_EXTERN void gemv_init_session(const struct gemv_t *gemv);

GEMV_EXTERN void gemv_device_malloc_(void **ptr, size_t size);
#define gemv_device_malloc(ptr, n)                                             \
  gemv_device_malloc_((void **)(ptr), sizeof(**(ptr)) * (n))

GEMV_EXTERN void gemv_device_free_(void **ptr);
#define gemv_device_free(ptr) gemv_device_free_((void **)(ptr))

typedef enum { GEMV_H2D = 0, GEMV_D2H = 1 } gemv_direction_t;
GEMV_EXTERN void gemv_copy_(void *dest, const void *src, size_t size,
                            gemv_direction_t direction);
#define gemv_copy(dest, src, size, direction)                                  \
  gemv_copy_((void *)(dest), (void *)(src), sizeof(*(src)) * (size), direction)

GEMV_EXTERN void gemv_run_(void *y, const void *x, const struct gemv_t *gemv);
#define gemv_run(y, x, gemv) gemv_run_((void *)y, (void *)x, gemv)

GEMV_EXTERN void gemv_finalize_session(void);

GEMV_EXTERN void gemv_finalize(struct gemv_t **gemv);

#endif // __GEMV_H__

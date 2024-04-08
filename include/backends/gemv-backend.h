#if !defined(__GEMV_BACKEND_H__)
#define __GEMV_BACKEND_H__

#include <string.h>

#include "gemv-impl.h"

static void gemv_convert(void *dest, const double *src, const size_t count,
                         const gemv_precision_t precision) {
  switch (precision) {
  case GEMV_FP64: memcpy(dest, src, sizeof(double) * count); break;
  case GEMV_FP32:
    for (size_t i = 0; i < count; i++) ((float *)dest)[i] = src[i];
    break;
  default: break;
  }
}

static size_t gemv_unit_size(const gemv_precision_t precision) {
  switch (precision) {
  case GEMV_FP64: return sizeof(double); break;
  case GEMV_FP32: return sizeof(float); break;
  default:
    gemv_log(GEMV_ERROR, "gemv_unit_size: Invlaid precision = %d", precision);
    return 0;
    break;
  }
}

static const char *gemv_precision_to_str(const gemv_precision_t precision) {
  switch (precision) {
  case GEMV_FP64: return "double"; break;
  case GEMV_FP32: return "float"; break;
  default:
    gemv_log(GEMV_ERROR, "gemv_precision_to_str: Invlaid precision = %d",
             precision);
    return NULL;
    break;
  }
}
#endif // __GEMV_BACKEND_H__

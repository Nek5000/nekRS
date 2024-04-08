#include "gemv.h"

#define GEMV_BACKEND(name)                                                     \
  GEMV_INTERN void gemv_register_##name(void) __attribute__((weak));           \
  void gemv_register_##name(void) { return; }

#include "backends/gemv-backend-list.h"

#undef GEMV_BACKEND

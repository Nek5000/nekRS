#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "gemv-impl.h"

static struct gemv_backend_t *backend_list = NULL;
static unsigned backend_count = 0, backend_max_count = 0;
static int backend_active = -1;

void gemv_backend_register(const char *name,
                           void (*init)(struct gemv_backend_t *backend,
                                        const struct gemv_t *gemv)) {
  gemv_log(GEMV_INFO, "gemv_backend_register: backend = %s", name);

  if (backend_count == backend_max_count) {
    backend_max_count += backend_max_count / 2 + 1;
    backend_list =
        realloc(backend_list, backend_max_count * sizeof(*backend_list));
  }

  strncpy(backend_list[backend_count].name, name, GEMV_MAX_BACKEND_LENGTH);
  backend_list[backend_count].init = init;
  backend_list[backend_count].malloc = NULL;
  backend_list[backend_count].free = NULL;
  backend_list[backend_count].copy = NULL;
  backend_list[backend_count].run = NULL;
  backend_list[backend_count].finalize = NULL;
  backend_count++;

  gemv_log(GEMV_INFO, "gemv_backend_register: backend_count = %d, done.",
           backend_count);
}

void gemv_backend_init(const struct gemv_t *gemv) {
  gemv_log(GEMV_INFO, "gemv_backend_init: backend = %d", gemv->backend);

  struct gemv_backend_t *backend = &backend_list[gemv->backend];
  backend->init(backend, gemv);
  backend_active = gemv->backend;

  gemv_log(GEMV_INFO, "gemv_backend_init: backend_active = %d, done.",
           backend_active);
}

void gemv_backend_run(void *y, const void *x, const struct gemv_t *gemv) {
  if (backend_active == -1)
    gemv_log(GEMV_ERROR, "gemv_backend_run: No active backend !");

  gemv_log(GEMV_INFO, "gemv_backend_run: y = %p, x = %p, gemv = %p", y, x,
           gemv);
  backend_list[backend_active].run(y, x, gemv);
  gemv_log(GEMV_INFO, "gemv_backend_run: done.");
}

void gemv_backend_malloc(void **ptr, const size_t size) {
  if (backend_active == -1)
    gemv_log(GEMV_ERROR, "gemv_backend_malloc: No active backend !");

  gemv_log(GEMV_INFO, "gemv_backend_malloc: size = %zu", size);
  backend_list[backend_active].malloc(ptr, size);
  gemv_log(GEMV_INFO, "gemv_backend_malloc: ptr = %p, done.", *ptr);
}

void gemv_backend_free(void **ptr) {
  if (backend_active == -1)
    gemv_log(GEMV_ERROR, "gemv_backend_malloc: No active backend !");

  gemv_log(GEMV_INFO, "gemv_backend_free: ptr = %p", *ptr);
  backend_list[backend_active].free(*ptr), *ptr = NULL;
  gemv_log(GEMV_INFO, "gemv_backend_free: done.");
}

void gemv_backend_copy(void *dest, const void *src, const size_t count,
                       const gemv_direction_t direction) {
  if (backend_active == -1)
    gemv_log(GEMV_ERROR, "gemv_backend_copy: No active backend !");

  gemv_log(
      GEMV_INFO,
      "gemv_backend_copy: src = %p, dest = %p, count = %zu, direction = %d",
      src, dest, count, direction);

  backend_list[backend_active].copy(dest, src, count, direction);

  gemv_log(GEMV_INFO, "gemv_backend_copy: done.");
}

void gemv_backend_finalize(void) {
  gemv_log(GEMV_INFO, "gemv_backend_finalize: backend_active = %d",
           backend_active);

  backend_list[backend_active].finalize();
  backend_active = -1;

  gemv_log(GEMV_INFO, "gemv_backend_finalize: done.", backend_active);
}

void gemv_backend_deregister(void) {
  gemv_log(GEMV_INFO, "gemv_backend_deregister: backend_count = %d",
           backend_count);

  for (unsigned i = 0; i < backend_count; i++)
    if (backend_list[i].finalize) backend_list[i].finalize();
  gemv_free(&backend_list);
  backend_count = backend_max_count = 0;

  gemv_log(GEMV_INFO, "gemv_backend_deregister: done.");
}

void gemv_set_backend_impl(struct gemv_t *gemv, const char *backend) {
  size_t backend_length = strnlen(backend, GEMV_MAX_BACKEND_LENGTH);
  char backend_lower[GEMV_MAX_BACKEND_LENGTH + 1];
  for (unsigned i = 0; i < backend_length; i++)
    backend_lower[i] = tolower(backend[i]);
  backend_lower[backend_length] = '\0';

  gemv_log(GEMV_INFO, "gemv_set_backend_impl: %s", backend_lower);

  gemv->backend = -1;
  for (unsigned i = 0; i < backend_count; i++) {
    if (strncmp(backend_lower, backend_list[i].name, GEMV_MAX_BACKEND_LENGTH))
      continue;
    gemv->backend = i;
    break;
  }

  if (gemv->backend == -1) {
    gemv_log(GEMV_ERROR, "gemv_set_backend_impl: backend \"%s\" not found.",
             backend_lower);
  }
  gemv_log(GEMV_INFO, "gemv_set_backend_impl: done.");
}

void gemv_free_(void **p) { free(*p), *p = NULL; }

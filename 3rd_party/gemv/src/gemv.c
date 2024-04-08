#include <ctype.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gemv-impl.h"

static void print_help(const char *name, int status) {
  FILE *fp = (status == EXIT_SUCCESS) ? stdout : stderr;
  fprintf(fp, "Usage: %s [OPTIONS]\n", name);
  fprintf(fp, "Options:\n");
  fprintf(fp,
          "  --gemv-verbose=<verbose level>, Verbose level (0, 1, 2, ...).\n");
  fprintf(fp, "  --gemv-device=<device id>, Device ID (0, 1, 2, ...).\n");
  fprintf(fp,
          "  --gemv-backend=<backend>, Backend (CUDA, HIP, OpenCL, etc.).\n");
  fprintf(fp, "  --gemv-help, Prints this help message and exit.\n");
  fflush(fp);
  exit(status);
}

static void parse_opts(struct gemv_t *gemv, int *argc, char ***argv_) {
  static struct option long_options[] = {
      {"gemv-verbose", optional_argument, 0, 10},
      {"gemv-device", optional_argument, 0, 20},
      {"gemv-backend", optional_argument, 0, 30},
      {"gemv-precision", optional_argument, 0, 40},
      {"gemv-help", no_argument, 0, 99},
      {0, 0, 0, 0}};

  // Default values for optional arguments.
  int device = GEMV_DEFAULT_DEVICE;
  gemv_verbose_t verbose = GEMV_DEFAULT_VERBOSE;
  gemv_precision_t precision = GEMV_DEFAULT_PRECISION;
  char backend[GEMV_MAX_BACKEND_LENGTH + 1] = GEMV_DEFAULT_BACKEND;

  if (argc == NULL || *argc == 0 || argv_ == NULL) goto set_options;

  char **argv = *argv_;
  for (;;) {
    int c = getopt_long(*argc, argv, "", long_options, NULL);
    if (c == -1) break;

    switch (c) {
    case 10: verbose = atoi(optarg); break;
    case 20: device = atoi(optarg); break;
    case 30: strncpy(backend, optarg, GEMV_MAX_BACKEND_LENGTH); break;
    case 40: precision = atoi(optarg); break;
    case 99: print_help(argv[0], EXIT_SUCCESS); break;
    default: print_help(argv[0], EXIT_FAILURE); break;
    }
  }

  for (int i = optind; i < *argc; i++) argv[i - optind] = argv[i];
  *argc -= optind;

set_options:
  gemv_set_verbose(verbose);
  gemv_set_device(gemv, device);
  gemv_set_backend(gemv, backend);
  gemv_set_precision(gemv, precision);
}

struct gemv_t *gemv_init(int *argc, char ***argv) {
  // Register all the backends.
#define GEMV_BACKEND(name) gemv_register_##name();
#include "backends/gemv-backend-list.h"
#undef GEMV_BACKEND

  // Initialize the gemv_t struct.
  struct gemv_t *gemv = gemv_calloc(struct gemv_t, 1);

  // Parse command line options if present.
  parse_opts(gemv, argc, (char ***)argv);

  return gemv;
}

void gemv_set_verbose(const gemv_verbose_t verbose) {
  gemv_set_verbose_impl(verbose);
  gemv_log(GEMV_INFO, "gemv_set_verbose: verbose = %d, done.", verbose);
}

void gemv_set_device(struct gemv_t *gemv, int device) {
  gemv_log(GEMV_INFO, "gemv_set_device: device = %d", device);
  gemv->device = device;
  gemv_log(GEMV_INFO, "gemv_set_device: done.");
}

void gemv_set_backend(struct gemv_t *gemv, const char *backend) {
  gemv_log(GEMV_INFO, "gemv_set_backend: backend = %s", backend);
  gemv_set_backend_impl(gemv, backend);
  gemv_log(GEMV_INFO, "gemv_set_backend: done.");
}

void gemv_set_matrix(struct gemv_t *gemv, const unsigned m, const unsigned n,
                     const double *A) {
  gemv_log(GEMV_INFO, "gemv_set_matrix: m = %u, n = %u", m, n);
  gemv->m = m, gemv->n = n;
  gemv->A = gemv_realloc(gemv->A, double, m *n);
  memcpy(gemv->A, A, sizeof(double) * m * n);
  gemv_log(GEMV_INFO, "gemv_set_matrix: done.");
}

void gemv_set_precision(struct gemv_t *gemv, const gemv_precision_t precision) {
  gemv_log(GEMV_INFO, "gemv_set_precision: precision = %d", precision);
  gemv->precision = precision;
  gemv_log(GEMV_INFO, "gemv_set_precision: done.");
}

void gemv_init_session(const struct gemv_t *gemv) {
  gemv_log(GEMV_INFO, "gemv_init_session: backend = %d", gemv->backend);
  gemv_backend_init(gemv);
  gemv_log(GEMV_INFO, "gemv_init_session: done.");
}

void gemv_device_malloc_(void **ptr, const size_t size) {
  gemv_log(GEMV_INFO, "gemv_device_malloc: size = %zu", size);
  gemv_backend_malloc(ptr, size);
  gemv_log(GEMV_INFO, "gemv_device_malloc: ptr = %p, done.", *ptr);
}

void gemv_device_free_(void **ptr) {
  gemv_log(GEMV_INFO, "gemv_device_free: ptr = %p", *ptr);
  gemv_backend_free(ptr);
  gemv_log(GEMV_INFO, "gemv_device_free: done.");
}

void gemv_copy_(void *dest, const void *src, size_t size,
                const gemv_direction_t direction) {
  gemv_log(GEMV_INFO,
           "gemv_copy: dest = %p, src = %p, size = %zu, direction = %d", dest,
           src, size, direction);
  gemv_backend_copy(dest, src, size, direction);
  gemv_log(GEMV_INFO, "gemv_copy: done.");
}

void gemv_run_(void *y, const void *x, const struct gemv_t *gemv) {
  gemv_log(GEMV_INFO, "gemv_run: y = %p, x = %p", y, x);
  gemv_backend_run(y, x, gemv);
  gemv_log(GEMV_INFO, "gemv_run: done.");
}

void gemv_finalize_session(void) {
  gemv_log(GEMV_INFO, "gemv_finalize_session: ..");
  gemv_backend_finalize();
  gemv_log(GEMV_INFO, "gemv_finalize_session: done.");
}

void gemv_finalize(struct gemv_t **gemv) {
  gemv_backend_deregister();
  gemv_free(&(*gemv)->A);
  gemv_free(gemv);
}

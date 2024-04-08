#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gemv-impl.h"

static gemv_verbose_t log_level = 0;
static const char *log_level_to_str[] = {"INFO", "WARNING", "ERROR"};

void gemv_set_verbose_impl(const gemv_verbose_t level) {
  if (level < GEMV_MUTE || level > GEMV_ERROR) {
    fprintf(stderr, "[ERROR] gemv_set_verbose: Invalid verbose level: %d\n",
            level);
    exit(EXIT_FAILURE);
  }

  log_level = level;
}

void gemv_log(const gemv_verbose_t level, const char *fmt, ...) {
  if (level <= GEMV_MUTE || log_level < level || level > GEMV_ERROR) return;

  char *fmt1 = gemv_calloc(char, 32 + strnlen(fmt, BUFSIZ));
  snprintf(fmt1, BUFSIZ, "[%s] ", log_level_to_str[level - 1]);
  fmt1 = strncat(fmt1, fmt, BUFSIZ);

  va_list args;
  va_start(args, fmt);
  char buf[BUFSIZ];
  vsnprintf(buf, BUFSIZ, fmt1, args);
  va_end(args);

  gemv_free(&fmt1);

  fprintf(stderr, "%s\n", buf);
  fflush(stderr);

  if (level == GEMV_ERROR) exit(EXIT_FAILURE);
}

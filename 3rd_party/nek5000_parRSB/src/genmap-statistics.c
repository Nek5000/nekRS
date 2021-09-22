#include <limits.h>
#include <time.h>

#include <genmap-impl.h>

#define MAXMETS 100
#define MAXLVLS 30
#define MAXSIZE (MAXMETS * MAXLVLS)

static double metrics[MAXMETS];
static double *stack;
static uint stack_size;

void metric_init() {
  uint i;
  for (i = 0; i < MAXMETS; i++)
    metrics[i] = 0.0;
  GenmapCalloc(MAXSIZE, &stack);
  stack_size = 0;
}

void metric_finalize() {
  if (stack != NULL)
    GenmapFree(stack);
}

void metric_acc(metric m, double count) { metrics[m] += count; }

void metric_tic(struct comm *c, metric m) {
  genmap_barrier(c);
  metrics[m] -= comm_time();
}

void metric_toc(struct comm *c, metric m) {
  metrics[m] += comm_time();
  genmap_barrier(c);
}

double metric_get_value(int level, metric m) {
  if (level == stack_size)
    return metrics[m];
  else if (level < stack_size)
    return stack[level * MAXMETS + m];
  return 0.0;
}

void metric_push_level() {
  assert(stack_size < MAXLVLS && "stack_size >= MAXLVLS");

  uint i;
  for (i = 0; i < MAXMETS; i++) {
    stack[stack_size * MAXMETS + i] = metrics[i];
    metrics[i] = 0.0;
  }
  stack_size++;
}

uint metric_get_levels() { return stack_size; }

void metric_print(struct comm *c) {
  double *min, *max, *sum, *buf;
  GenmapCalloc(MAXSIZE, &min);
  GenmapCalloc(MAXSIZE, &max);
  GenmapCalloc(MAXSIZE, &sum);
  GenmapCalloc(MAXSIZE, &buf);

  uint max_size = stack_size * MAXMETS;
  assert(max_size <= MAXSIZE);

  uint i;
  for (i = 0; i < max_size; i++)
    min[i] = max[i] = sum[i] = stack[i];

  comm_allreduce(c, gs_double, gs_min, min, MAXSIZE, buf); // min
  comm_allreduce(c, gs_double, gs_max, max, MAXSIZE, buf); // max
  comm_allreduce(c, gs_double, gs_add, sum, MAXSIZE, buf); // sum
  for (i = 0; i < max_size; i++)
    sum[i] /= c->np;

#define SUMMARY(i, m)                                                          \
  sum[i * MAXMETS + m], min[i * MAXMETS + m], max[i * MAXMETS + m]

  int j;
  for (i = 0; i < stack_size; i++) {
    if (c->id == 0) {
      printf("level=%02d\n", i);
      printf("  RCB                    : %g/%g/%g\n", SUMMARY(i, RCB));
      printf("  FIEDLER                : %g/%g/%g\n", SUMMARY(i, FIEDLER));
      printf("  FIEDLER_NITER          : %g/%g/%g\n",
             SUMMARY(i, FIEDLER_NITER));
      printf("    LAPLACIAN_INIT       : %g/%g/%g\n",
             SUMMARY(i, LAPLACIAN_INIT));
      printf("    RQI                  : %g/%g/%g\n", SUMMARY(i, RQI));
      printf("    RQI_NITER            : %g/%g/%g\n", SUMMARY(i, RQI_NITER));
      for (j = 0; j < min[i * MAXMETS + RQI_NITER]; j++)
        printf("      rqi=%02d             : %g/%g/%g\n", j,
               SUMMARY(i, END + j));
      printf("      PROJECT            : %g/%g/%g\n", SUMMARY(i, PROJECT));
      printf("      PROJECT_NITER      : %g/%g/%g\n",
             SUMMARY(i, PROJECT_NITER));
      printf("        LAPLACIAN        : %g/%g/%g\n", SUMMARY(i, LAPLACIAN));
      printf("        PROJECT          : %g/%g/%g\n", SUMMARY(i, PROJECT));
    }
  }

  GenmapFree(min);
  GenmapFree(max);
  GenmapFree(sum);
  GenmapFree(buf);

#undef SUMMARY
}

#undef MAXMETS
#undef MAXLVLS
#undef MAXSIZE

#include <limits.h>
#include <math.h>
#include <stdio.h>
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
  comm_barrier(c);
  metrics[m] -= comm_time();
}

void metric_toc(struct comm *c, metric m) {
  metrics[m] += comm_time();
  comm_barrier(c);
}

double metric_get_value(int level, metric m) {
  if (level == stack_size)
    return metrics[m];
  else if (level < stack_size)
    return stack[level * MAXMETS + m];
  else // FIXME: May be an assert failure?
    return -1.0;
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
      printf("  WEIGHTEDLAPLACIANSETUP : %g/%g/%g\n",
             SUMMARY(i, WEIGHTEDLAPLACIANSETUP));
      printf("  FIEDLER                : %g/%g/%g\n", SUMMARY(i, FIEDLER));
      printf("  NFIEDLER               : %g/%g/%g\n", SUMMARY(i, NFIEDLER));
      printf("    LAPLACIANSETUP       : %g/%g/%g\n",
             SUMMARY(i, LAPLACIANSETUP));
      printf("      FINDNBRS           : %g/%g/%g\n", SUMMARY(i, FINDNBRS));
      printf("      CSRMATSETUP        : %g/%g/%g\n", SUMMARY(i, CSRMATSETUP));
      printf("      CSRTOPSETUP        : %g/%g/%g\n", SUMMARY(i, CSRTOPSETUP));
      printf("    PRECONDSETUP         : %g/%g/%g\n", SUMMARY(i, PRECONDSETUP));
      printf("    RQI                  : %g/%g/%g\n", SUMMARY(i, RQI));
      printf("    NRQI                 : %g/%g/%g\n", SUMMARY(i, NRQI));
      for (j = 0; j < min[i * MAXMETS + NRQI]; j++)
        printf("      rqi=%02d             : %g/%g/%g\n", j,
               SUMMARY(i, END + j));
      printf("      PROJECT            : %g/%g/%g\n", SUMMARY(i, PROJECT));
      printf("      NPROJECT           : %g/%g/%g\n", SUMMARY(i, NPROJECT));
      printf("        VCYCLE           : %g/%g/%g\n", SUMMARY(i, VCYCLE));
      printf("        LAPLACIAN        : %g/%g/%g\n", SUMMARY(i, LAPLACIAN));
      printf("        PROJECT          : %g/%g/%g\n", SUMMARY(i, PROJECT));
      printf("      GRAMMIAN           : %g/%g/%g\n", SUMMARY(i, GRAMMIAN));
      printf("  FIEDLERSORT            : %g/%g/%g\n", SUMMARY(i, FIEDLERSORT));
      printf("  BISECTANDREPAIR        : %g/%g/%g\n",
             SUMMARY(i, BISECTANDREPAIR));
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

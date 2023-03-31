#include "metrics.h"
#include "parrsb-impl.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

extern int rsb(struct array *elements, int nv, int check,
               parrsb_options *options, struct comm *gc, buffer *bfr);
extern int rcb(struct array *elements, size_t unit_size, int ndim,
               struct comm *ci, buffer *bfr);

parrsb_options parrsb_default_options = {
    // General options
    .partitioner = 0,
    .verbose_level = 0,
    .profile_level = 0,
    .two_level = 1,
    .repair = 0,
    // RSB common (Lanczos + MG) options
    .rsb_algo = 0,
    .rsb_pre = 1,
    .rsb_max_iter = 50,
    .rsb_max_passes = 50,
    .rsb_tol = 1e-5,
    // RSB MG specific options
    .rsb_mg_grammian = 0,
    .rsb_mg_factor = 2,
    .rsb_mg_sagg = 0};

static char *ALGO[3] = {"RSB", "RCB", "RIB"};

#define UPDATE_OPTION(OPT, STR, IS_INT)                                        \
  do {                                                                         \
    const char *val = getenv(STR);                                             \
    if (val != NULL) {                                                         \
      if (IS_INT)                                                              \
        options->OPT = atoi(val);                                              \
      else                                                                     \
        options->OPT = atof(val);                                              \
    }                                                                          \
  } while (0)

static void update_options(parrsb_options *options) {
  UPDATE_OPTION(partitioner, "PARRSB_PARTITIONER", 1);
  UPDATE_OPTION(verbose_level, "PARRSB_VERBOSE_LEVEL", 1);
  UPDATE_OPTION(profile_level, "PARRSB_PROFILE_LEVEL", 1);
  UPDATE_OPTION(two_level, "PARRSB_TWO_LEVEL", 1);
  UPDATE_OPTION(repair, "PARRSB_REPAIR", 1);
  UPDATE_OPTION(rsb_algo, "PARRSB_RSB_ALGO", 1);
  UPDATE_OPTION(rsb_pre, "PARRSB_RSB_PRE", 1);
  UPDATE_OPTION(rsb_max_iter, "PARRSB_RSB_MAX_ITER", 1);
  UPDATE_OPTION(rsb_max_passes, "PARRSB_RSB_MAX_PASSES", 1);
  UPDATE_OPTION(rsb_tol, "PARRSB_RSB_TOL", 0);
  UPDATE_OPTION(rsb_mg_grammian, "PARRSB_RSB_MG_GRAMMIAN", 1);
  UPDATE_OPTION(rsb_mg_factor, "PARRSB_RSB_MG_FACTOR", 1);
  UPDATE_OPTION(rsb_mg_sagg, "PARRSB_RSB_MG_SMOOTH_AGGREGATION", 1);
  if (options->verbose_level == 0) 
    options->profile_level = 0;
}

#undef UPDATE_OPTION

#define PRINT_OPTION(OPT, STR, FMT) printf("%s = " FMT "\n", STR, options->OPT)

static void print_options(parrsb_options *options) {
  PRINT_OPTION(partitioner, "PARRSB_PARTITIONER", "%d");
  PRINT_OPTION(verbose_level, "PARRSB_VERBOSE_LEVEL", "%d");
  PRINT_OPTION(profile_level, "PARRSB_PROFILE_LEVEL", "%d");
  PRINT_OPTION(two_level, "PARRSB_TWO_LEVEL", "%d");
  PRINT_OPTION(repair, "PARRSB_REPAIR", "%d");
  PRINT_OPTION(rsb_algo, "PARRSB_RSB_ALGO", "%d");
  PRINT_OPTION(rsb_pre, "PARRSB_RSB_PRE", "%d");
  PRINT_OPTION(rsb_max_iter, "PARRSB_RSB_MAX_ITER", "%d");
  PRINT_OPTION(rsb_max_passes, "PARRSB_RSB_MAX_PASSES", "%d");
  PRINT_OPTION(rsb_tol, "PARRSB_RSB_TOL", "%lf");
  PRINT_OPTION(rsb_mg_grammian, "PARRSB_RSB_MG_GRAMMIAN", "%d");
  PRINT_OPTION(rsb_mg_factor, "PARRSB_RSB_MG_FACTOR", "%d");
  PRINT_OPTION(rsb_mg_sagg, "PARRSB_RSB_MG_SMOOTH_AGGREGATION", "%d");
}

#undef PRINT_OPTION

static size_t load_balance(struct array *elist, uint nel, int nv, double *coord,
                           long long *vtx, struct crystal *cr, buffer *bfr) {
  struct comm *c = &cr->comm;
  slong out[2][1], wrk[2][1], in = nel;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0], nelg = out[1][0];

  uint nstar = nelg / c->np, nrem = nelg - nstar * c->np;
  slong lower = (nstar + 1) * nrem;

  size_t unit_size;
  if (vtx == NULL) // RCB
    unit_size = sizeof(struct rcb_element);
  else // RSB
    unit_size = sizeof(struct rsb_element);

  array_init_(elist, nel, unit_size, __FILE__, __LINE__);

  struct rcb_element *pe = (struct rcb_element *)calloc(1, unit_size);
  pe->origin = c->id;

  int ndim = (nv == 8) ? 3 : 2;
  for (uint e = 0; e < nel; ++e) {
    slong eg = pe->globalId = start + e + 1;
    if (eg <= lower)
      pe->proc = (eg - 1) / (nstar + 1);
    else if (nstar != 0)
      pe->proc = (eg - 1 - lower) / nstar + nrem;

    pe->coord[0] = pe->coord[1] = pe->coord[2] = 0.0;
    for (int v = 0; v < nv; v++)
      for (int n = 0; n < ndim; n++)
        pe->coord[n] += coord[e * ndim * nv + v * ndim + n];
    for (int n = 0; n < ndim; n++)
      pe->coord[n] /= nv;

    array_cat_(unit_size, elist, pe, 1, __FILE__, __LINE__);
  }

  if (vtx != NULL) { // RSB
    struct rsb_element *pr = (struct rsb_element *)elist->ptr;
    for (uint e = 0; e < nel; e++) {
      for (int v = 0; v < nv; v++)
        pr[e].vertices[v] = vtx[e * nv + v];
    }
  }

  sarray_transfer_(elist, unit_size, offsetof(struct rcb_element, proc), 1, cr);
  if (vtx == NULL) // RCB
    sarray_sort(struct rcb_element, elist->ptr, elist->n, globalId, 1, bfr);
  else // RSB
    sarray_sort(struct rsb_element, elist->ptr, elist->n, globalId, 1, bfr);

  free(pe);
  return unit_size;
}

static void restore_original(int *part, int *seq, struct crystal *cr,
                             struct array *elist, size_t usize, buffer *bfr) {
  sarray_transfer_(elist, usize, offsetof(struct rcb_element, origin), 1, cr);
  uint nel = elist->n;

  if (usize == sizeof(struct rsb_element)) // RSB
    sarray_sort(struct rsb_element, elist->ptr, nel, globalId, 1, bfr);
  else if (usize == sizeof(struct rcb_element)) // RCB
    sarray_sort(struct rcb_element, elist->ptr, nel, globalId, 1, bfr);

  struct rcb_element *element;
  uint e;
  for (e = 0; e < nel; e++) {
    element = (struct rcb_element *)((char *)elist->ptr + e * usize);
    part[e] = element->origin; // element[e].origin;
  }

  if (seq != NULL) {
    for (e = 0; e < nel; e++) {
      element = (struct rcb_element *)((char *)elist->ptr + e * usize);
      seq[e] = element->seq; // element[e].seq;
    }
  }
}

int parrsb_part_mesh(int *part, int *seq, long long *vtx, double *coord,
                     int nel, int nv, parrsb_options options, MPI_Comm comm) {
  update_options(&options);

  struct comm c;
  comm_init(&c, comm);
  if (c.id == 0 && options.verbose_level > 0) {
    printf("Running parRSB ...\n");
    print_options(&options);
    fflush(stdout);
  }

  parrsb_barrier(&c);
  double t = comm_time();

  struct crystal cr;
  crystal_init(&cr, &c);

  buffer bfr;
  buffer_init(&bfr, nel * sizeof(struct rsb_element));

  // Load balance input data
  struct array elist;
  size_t esize = load_balance(&elist, nel, nv, coord, vtx, &cr, &bfr);

  struct comm ca;
  comm_split(&c, nel > 0, c.id, &ca);

  // Run RSB now
  metric_init();
  if (nel > 0) {
    slong out[2][1], wrk[2][1], in = nel;
    comm_scan(out, &ca, gs_long, gs_add, &in, 1, wrk);
    slong nelg = out[1][0];

    if (ca.np > nelg) {
      if (ca.id == 0)
        printf("Total number of elements is smaller than the "
               "number of processors.\n"
               "Run with smaller number of processors.\n");
      return 1;
    }

    int ndim = (nv == 8) ? 3 : 2;
    switch (options.partitioner) {
    case 0:
      rsb(&elist, nv, 1, &options, &ca, &bfr);
      break;
    case 1:
      rcb(&elist, esize, ndim, &ca, &bfr);
      break;
    case 2:
      rib(&elist, esize, ndim, &ca, &bfr);
      break;
    default:
      break;
    }

    metric_rsb_print(&ca, options.profile_level);
  }
  metric_finalize();
  comm_free(&ca);

  restore_original(part, seq, &cr, &elist, esize, &bfr);

  // Report time and finish
  parrsb_barrier(&c);
  if (c.id == 0) {
    printf("par%s finished in %g s\n", ALGO[options.partitioner],
           comm_time() - t);
    fflush(stdout);
  }

  array_free(&elist), buffer_free(&bfr);
  crystal_free(&cr), comm_free(&c);

  return 0;
}

void fparrsb_part_mesh(int *part, int *seq, long long *vtx, double *coord,
                       int *nel, int *nv, int *options, int *comm, int *err) {
  *err = 1;
  comm_ext c = MPI_Comm_f2c(*comm);
  parrsb_options opt = parrsb_default_options;
  *err = parrsb_part_mesh(part, seq, vtx, coord, *nel, *nv, opt, c);
}

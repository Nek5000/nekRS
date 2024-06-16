#include "metrics.h"
#include "parrsb-impl.h"

#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))

void parrsb_print(const struct comm *c, int verbose, const char *fmt, ...) {
  comm_barrier(c);

  va_list vargs;
  if (c->id == 0 && verbose > 0) {
    va_start(vargs, fmt);
    vprintf(fmt, vargs);
    va_end(vargs);
    printf("\n");
    fflush(stdout);
  }
}

parrsb_options parrsb_default_options = {
    // General options
    .partitioner = 0,
    .tagged = 0,
    .levels = 2,
    .find_disconnected_comps = 1,
    .repair = 0,
    .verbose_level = 1,
    .profile_level = 0,
    // RSB common (Lanczos and MG) options
    .rsb_algo = 0,
    .rsb_pre = 1,
    .rsb_max_iter = 50,
    .rsb_max_passes = 50,
    .rsb_tol = 1e-5,
    .rsb_dump_stats = 0,
    // RSB MG specific options
    .rsb_mg_grammian = 0,
    .rsb_mg_factor = 2};

static char *ALGO[3] = {"RSB", "RCB", "RIB"};

static void update_options(parrsb_options *const options) {
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

  UPDATE_OPTION(partitioner, "PARRSB_PARTITIONER", 1);
  UPDATE_OPTION(tagged, "PARRSB_TAGGED", 1);
  UPDATE_OPTION(levels, "PARRSB_LEVELS", 1);
  UPDATE_OPTION(find_disconnected_comps, "PARRSB_FIND_DISCONNECTED_COMPONENTS",
                1);
  UPDATE_OPTION(repair, "PARRSB_REPAIR", 1);
  UPDATE_OPTION(verbose_level, "PARRSB_VERBOSE_LEVEL", 1);
  UPDATE_OPTION(profile_level, "PARRSB_PROFILE_LEVEL", 1);
  UPDATE_OPTION(rsb_algo, "PARRSB_RSB_ALGO", 1);
  UPDATE_OPTION(rsb_pre, "PARRSB_RSB_PRE", 1);
  UPDATE_OPTION(rsb_max_iter, "PARRSB_RSB_MAX_ITER", 1);
  UPDATE_OPTION(rsb_max_passes, "PARRSB_RSB_MAX_PASSES", 1);
  UPDATE_OPTION(rsb_tol, "PARRSB_RSB_TOL", 0);
  UPDATE_OPTION(rsb_mg_grammian, "PARRSB_RSB_MG_GRAMMIAN", 1);
  UPDATE_OPTION(rsb_mg_factor, "PARRSB_RSB_MG_FACTOR", 1);

#undef UPDATE_OPTION
}

static void print_options(const struct comm *c,
                          const parrsb_options *const options) {
#define PRINT_OPTION(OPT, STR, FMT)                                            \
  parrsb_print(c, options->verbose_level, "%s = " FMT "", STR, options->OPT)

  PRINT_OPTION(partitioner, "PARRSB_PARTITIONER", "%d");
  PRINT_OPTION(tagged, "PARRSB_TAGGED", "%d");
  PRINT_OPTION(levels, "PARRSB_LEVELS", "%d");
  PRINT_OPTION(find_disconnected_comps, "PARRSB_FIND_DISCONNECTED_COMPONENTS",
               "%d");
  PRINT_OPTION(repair, "PARRSB_REPAIR", "%d");
  PRINT_OPTION(verbose_level, "PARRSB_VERBOSE_LEVEL", "%d");
  PRINT_OPTION(profile_level, "PARRSB_PROFILE_LEVEL", "%d");
  PRINT_OPTION(rsb_algo, "PARRSB_RSB_ALGO", "%d");
  PRINT_OPTION(rsb_pre, "PARRSB_RSB_PRE", "%d");
  PRINT_OPTION(rsb_max_iter, "PARRSB_RSB_MAX_ITER", "%d");
  PRINT_OPTION(rsb_max_passes, "PARRSB_RSB_MAX_PASSES", "%d");
  PRINT_OPTION(rsb_tol, "PARRSB_RSB_TOL", "%lf");
  PRINT_OPTION(rsb_mg_grammian, "PARRSB_RSB_MG_GRAMMIAN", "%d");
  PRINT_OPTION(rsb_mg_factor, "PARRSB_RSB_MG_FACTOR", "%d");

#undef PRINT_OPTION
}

static size_t load_balance(struct array *elist, uint nel, int nv,
                           const double *const xyz, const long long *const vtx,
                           int verbose, struct crystal *cr, buffer *bfr) {
  struct comm *c = &cr->comm;
  slong out[2][1], wrk[2][1], in = nel;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0], nelg = out[1][0];
  parrsb_print(c, verbose, "load_balance: nelg = %lld", nelg);

  uint nstar = nelg / c->np, nrem = nelg - nstar * c->np;
  slong lower = (nstar + 1) * nrem;

  size_t unit_size;
  if (vtx == NULL) // RCB
    unit_size = sizeof(struct rcb_element);
  else             // RSB
    unit_size = sizeof(struct rsb_element);
  parrsb_print(
      c, verbose, "load_balance: unit_size = %zu (rsb = %zu, rcb = %zu)",
      unit_size, sizeof(struct rsb_element), sizeof(struct rcb_element));

  array_init_(elist, nel, unit_size, __FILE__, __LINE__);

  struct rcb_element *pe = (struct rcb_element *)calloc(1, unit_size);
  pe->origin = c->id;

  int ndim = (nv == 8) ? 3 : 2;
  for (uint e = 0; e < nel; ++e) {
    slong eg = pe->globalId = start + e + 1;
    if (nstar == 0)
      pe->proc = eg - 1;
    else if (eg <= lower)
      pe->proc = (eg - 1) / (nstar + 1);
    else
      pe->proc = (eg - 1 - lower) / nstar + nrem;

    pe->coord[0] = pe->coord[1] = pe->coord[2] = 0.0;
    if (xyz != NULL) {
      for (int v = 0; v < nv; v++)
        for (int n = 0; n < ndim; n++)
          pe->coord[n] += xyz[e * ndim * nv + v * ndim + n];
      for (int n = 0; n < ndim; n++) pe->coord[n] /= nv;
    }

    array_cat_(unit_size, elist, pe, 1, __FILE__, __LINE__);
  }

  if (vtx != NULL) { // RSB
    struct rsb_element *pr = (struct rsb_element *)elist->ptr;
    for (uint e = 0; e < nel; e++) {
      for (int v = 0; v < nv; v++) pr[e].vertices[v] = vtx[e * nv + v];
    }
  }

  sarray_transfer_(elist, unit_size, offsetof(struct rcb_element, proc), 1, cr);
  if (vtx == NULL) // RCB
    sarray_sort(struct rcb_element, elist->ptr, elist->n, globalId, 1, bfr);
  else             // RSB
    sarray_sort(struct rsb_element, elist->ptr, elist->n, globalId, 1, bfr);

  free(pe);

  return unit_size;
}

static void restore_original(int *part, struct crystal *cr, struct array *elist,
                             size_t usize, buffer *bfr) {
  sarray_transfer_(elist, usize, offsetof(struct rcb_element, origin), 1, cr);
  uint nel = elist->n;

  if (usize == sizeof(struct rsb_element))      // RSB
    sarray_sort(struct rsb_element, elist->ptr, nel, globalId, 1, bfr);
  else if (usize == sizeof(struct rcb_element)) // RCB
    sarray_sort(struct rcb_element, elist->ptr, nel, globalId, 1, bfr);

  struct rcb_element *element;
  uint e;
  for (e = 0; e < nel; e++) {
    element = (struct rcb_element *)((char *)elist->ptr + e * usize);
    part[e] = element->origin; // element[e].origin;
  }
}

static void initialize_node_aux(struct comm *c, const struct comm *const gc) {
#ifdef MPI
  MPI_Comm node;
  MPI_Comm_split_type(gc->c, MPI_COMM_TYPE_SHARED, gc->id, MPI_INFO_NULL,
                      &node);
  comm_init(c, node);
  MPI_Comm_free(&node);
#else
  comm_init(1, 1);
#endif
}

static void initialize_levels(struct comm *const comms, int *const levels_in,
                              const struct comm *const c, const int verbose) {
  // Level 1 communicator is the global communicator.
  comm_dup(&comms[0], c);
  // Node level communicator is the last level communicator.
  struct comm nc;
  initialize_node_aux(&nc, c);

  // Find the number of nodes under the global communicator and number of MPI
  // ranks in the node level communicator.
  uint num_nodes, nranks_per_node;
  {
    sint in = (nc.id == 0), wrk;
    comm_allreduce(c, gs_int, gs_add, &in, 1, &wrk);
    num_nodes = in;

    nranks_per_node = nc.np;
    // Check invariant: nranks_per_node should be the same across all the nodes.
    sint nranks_max = nranks_per_node, nranks_min = nranks_per_node;
    comm_allreduce(&comms[0], gs_int, gs_max, &nranks_max, 1, &wrk);
    comm_allreduce(&comms[0], gs_int, gs_min, &nranks_min, 1, &wrk);
    assert(nranks_max == nranks_min);
    // Check invariant: nranks_per_node must be larger than 0.
    assert(nranks_per_node > 0);
    parrsb_print(c, verbose,
                 "initialize_levels: num_nodes = %u, num_ranks_per_node = %u",
                 num_nodes, nranks_per_node);
  }

  // Hardcode the maximum number of levels to two for now.
  sint levels = 2;
  *levels_in = levels = MIN(levels, *levels_in);
  if (levels > 1) comm_dup(&comms[levels - 1], &nc);
  comm_free(&nc);

  parrsb_print(c, verbose, "initialize_levels: levels = %u", levels);
}

static void parrsb_part_mesh_v0(int *part, const long long *const vtx,
                                const double *const xyz, const uint nel,
                                const unsigned nv,
                                parrsb_options *const options,
                                const struct comm *const c,
                                struct crystal *const cr, buffer *const bfr) {
  const int verbose = options->verbose_level;

  if (vtx == NULL && xyz == NULL) {
    parrsb_print(
        c, verbose,
        "parrsb_part_mesh_v0: Both vertices and coordinates can't be NULL");
    MPI_Abort(c->c, EXIT_FAILURE);
  }
  if (xyz == NULL) options->rsb_pre = 0;

  struct array elist;
  size_t esize = load_balance(&elist, nel, nv, xyz, vtx, verbose, cr, bfr);

  struct comm ca;
  comm_split(c, elist.n > 0, c->id, &ca);

  // Setup communicators for each level of the partitioning.
  struct comm comms[9];
  {
    // Check invariant: levels > 0 and levels <= sizeof(comms) /
    // sizeof(comms[0]).
    const uint levels = options->levels;
    assert(levels <= sizeof(comms) / sizeof(comms[0]));
    initialize_levels(comms, &options->levels, &ca, verbose);
    parrsb_print(
        c, verbose,
        "parrsb_part_mesh_v0: levels requested = %d, levels enabled = %d",
        levels, options->levels);
  }

  parrsb_print(c, verbose, "parrsb_part_mesh_v0: running partitioner ...");
  if (elist.n > 0) {
    int ndim = (nv == 8) ? 3 : 2;
    switch (options->partitioner) {
    case 0: rsb(&elist, nv, options, comms, bfr); break;
    case 1: rcb(&elist, esize, ndim, &ca, bfr); break;
    case 2: rib(&elist, esize, ndim, &ca, bfr); break;
    default: break;
    }
  }
  comm_free(&ca);

  for (uint l = 0; l < (uint)options->levels; l++) comm_free(&comms[l]);

  parrsb_print(c, verbose, "parrsb_part_mesh_v0: restore original input");
  restore_original(part, cr, &elist, esize, bfr);

  array_free(&elist);
}

void parrsb_check_tagged_partitions(const long long *const eids,
                                    const long long *const vtx, const uint nel,
                                    const unsigned nv, const uint ntags,
                                    const struct comm *const c,
                                    const int verbose) {
  parrsb_print(c, verbose, "Check if the input elements are sorted locally.");
  {
    sint sorted = 1;
    for (uint i = 1; i < nel; i++) {
      if (eids[i] < eids[i - 1]) {
        sorted = 0;
        break;
      }
    }

    sint wrk;
    comm_allreduce(c, gs_int, gs_min, &sorted, 1, &wrk);
    if (!sorted) {
      if (c->id == 0) {
        fprintf(stderr, "Input elements are not sorted.\n");
        fflush(stderr);
      }
      exit(EXIT_FAILURE);
    }
  }

  // Number the elements within the each tag id and setup a gs handle based on
  // 2D element id.
  parrsb_print(c, verbose, "Number elements within each layer.");
  const uint tag_id = c->id / ntags;
  struct comm lc;
  struct gs_data *gse = NULL;
  {
    comm_split(c, tag_id, c->id, &lc);

    slong out[2][1], wrk[2][1], in = nel;
    comm_scan(out, &lc, gs_long, gs_add, &in, 1, wrk);
    slong start = out[0][0];

    slong *lids = tcalloc(slong, nel);
    for (uint i = 0; i < nel; i++) lids[i] = start + i;

    gse = gs_setup(lids, nel, c, 0, gs_pairwise, 0);
    free(lids);
  }

  // Setup a local gs handle based on the original gs vertex ids.
  parrsb_print(c, verbose, "Setup multiplicity.");
  const size_t size = nel * nv;
  buffer bfr;
  buffer_init(&bfr, size);
  sint *mul = tcalloc(sint, size);
  {
    struct gs_data *gsl = gs_setup(vtx, size, &lc, 0, gs_pairwise, 0);
    for (uint i = 0; i < size; i++) mul[i] = 1;
    gs(mul, gs_int, gs_add, 0, gsl, &bfr);
    gs_free(gsl);
  }

  // Now let's compare the multiplicity across the layers.
  parrsb_print(c, verbose, "Check multiplicity across the layers.");
  {
    sint *lmin = tcalloc(sint, nel);
    sint *lmax = tcalloc(sint, nel);
    for (uint v = 0; v < nv; v++) {
      for (uint e = 0; e < nel; e++) {
        lmin[e] = mul[e * nv + v];
        lmax[e] = mul[e * nv + v];
      }

      gs(lmin, gs_int, gs_min, 0, gse, &bfr);
      gs(lmax, gs_int, gs_max, 0, gse, &bfr);

      for (uint e = 0; e < nel; e++) assert(lmin[e] == lmax[e]);
    }

    free(lmin), free(lmax);
  }

  free(mul);
  buffer_free(&bfr);
  gs_free(gse);
  comm_free(&lc);

  return;
}

static void parrsb_part_mesh_v1(int *part, const long long *const vtx,
                                const double *const xyz, const int *const tag,
                                const uint nel, const unsigned nv,
                                parrsb_options *const options,
                                const struct comm *const c,
                                struct crystal *const cr, buffer *const bfr) {
  const int verbose = options->verbose_level;
  parrsb_print(c, verbose, "Find number of tags in the mesh ...");

  struct tag_t {
    uint p, tag, seq, tagn;
  };

  struct array tags;
  array_init(struct tag_t, &tags, nel);

  {
    struct tag_t tt;
    for (uint i = 0; i < nel; i++) {
      tt.seq = i, tt.tag = tag[i], tt.p = tt.tag % c->np;
      array_cat(struct tag_t, &tags, &tt, 1);
    }
    sarray_sort(struct tag_t, tags.ptr, tags.n, tag, 0, bfr);
  }

  struct array unique;
  array_init(struct tag_t, &unique, 1024);

  if (tags.n > 0) {
    const struct tag_t *const pt = (const struct tag_t *const)tags.ptr;
    array_cat(struct tag_t, &unique, &pt[0], 1);
    for (uint i = 1; i < tags.n; i++) {
      if (pt[i].tag > pt[i - 1].tag)
        array_cat(struct tag_t, &unique, &pt[i], 1);
    }
  }

  sint out[2][1];
  {
    sarray_transfer(struct tag_t, &unique, p, 1, cr);
    sarray_sort(struct tag_t, unique.ptr, unique.n, tag, 0, bfr);

    const struct tag_t *const pu = (const struct tag_t *const)unique.ptr;
    sint in = 0;
    if (unique.n > 0) {
      in = 1;
      for (uint i = 1; i < unique.n; i++) {
        if (pu[i].tag > pu[i - 1].tag) in++;
      }
    }

    sint wrk[2][1];
    comm_scan(out, c, gs_int, gs_add, &in, 1, wrk);
  }
  const uint num_tags = out[1][0], tag_start = out[0][0];

  parrsb_print(c, verbose, "Num tags: %d", num_tags);
  if (c->np % num_tags != 0) {
    if (c->id == 0) {
      fprintf(stderr,
              "Number of processes must be a multiple of number of tags: "
              "processes = %d, tags = %d.\n",
              c->np, num_tags);
    }
    exit(EXIT_FAILURE);
  }

  {
    struct tag_t *const pu = (struct tag_t *const)unique.ptr;
    uint start = tag_start;
    if (unique.n > 0) {
      pu[0].tagn = start;
      for (uint i = 1; i < unique.n; i++) {
        if (pu[i].tag > pu[i - 1].tag) start++;
        pu[i].tagn = start;
      }
    }

    sarray_transfer(struct tag_t, &unique, p, 0, cr);
    sarray_sort(struct tag_t, unique.ptr, unique.n, tag, 0, bfr);
  }

  const uint chunk_size = c->np / num_tags;
  parrsb_print(c, verbose, "Processes per tag: %d", chunk_size);
  {
    struct tag_t *const pt = (struct tag_t *const)tags.ptr;
    const struct tag_t *const pu = (const struct tag_t *const)unique.ptr;
    for (uint i = 0, s = 0; i < unique.n; i++) {
      uint e = s + 1;
      assert(pt[s].tag == pu[i].tag);
      while (e < tags.n && pt[e].tag == pu[i].tag) e++;
      for (uint j = s; j < e; j++)
        pt[j].p = chunk_size * pu[i].tagn + pt[i].seq % chunk_size;
      s = e;
    }

    sarray_sort(struct tag_t, tags.ptr, tags.n, seq, 0, bfr);
  }
  array_free(&unique);

  struct element_t {
    uint proc, part, seq;
    scalar xyz[MAXDIM * MAXNV];
    slong vertices[MAXNV];
  };

  struct array elements;
  array_init(struct element_t, &elements, nel);

  parrsb_print(c, verbose,
               "Pack element data for transfering. tags.n=%u, nel=%u", tags.n,
               nel);
  const unsigned ndim = (nv == 8) ? 3 : 2;
  {
    assert(tags.n == nel);
    const struct tag_t *const pt = (const struct tag_t *const)tags.ptr;
    struct element_t et;
    for (uint i = 0; i < tags.n; i++) {
      et.proc = pt[i].p, et.seq = i;
      for (uint j = 0; j < nv; j++) {
        et.vertices[j] = vtx[i * nv + j];
        for (uint k = 0; k < ndim; k++)
          et.xyz[j * ndim + k] = xyz[i * nv * ndim + j * ndim + k];
      }
      array_cat(struct element_t, &elements, &et, 1);
    }

    sarray_transfer(struct element_t, &elements, proc, 1, cr);
  }
  array_free(&tags);

  parrsb_print(c, verbose, "Copy element data for feeding to parRSB.");
  long long *lvtx = tcalloc(long long, (elements.n + 1) * nv);
  double *lxyz = tcalloc(double, (elements.n + 1) * nv * ndim);
  {
    const struct element_t *const pe =
        (const struct element_t *const)elements.ptr;
    for (uint e = 0; e < elements.n; e++) {
      for (uint j = 0; j < nv; j++) {
        lvtx[e * nv + j] = pe[e].vertices[j];
        for (uint k = 0; k < ndim; k++)
          lxyz[e * nv * ndim + j * ndim + k] = pe[e].xyz[j * ndim + k];
      }
    }
  }

  parrsb_print(c, verbose, "Run parRSB locally within a tag now.");
  {
    int *lpart = tcalloc(int, elements.n + 1);

    struct comm lc;
    comm_split(c, c->id / chunk_size, c->id, &lc);

    struct crystal lcr;
    crystal_init(&lcr, &lc);

    options->verbose_level = 0;
    options->profile_level = 0;
    parrsb_part_mesh_v0(lpart, lvtx, lxyz, elements.n, nv, options, &lc, &lcr,
                        bfr);
    crystal_free(&lcr), comm_free(&lc);

    struct element_t *const pe = (struct element_t *const)elements.ptr;
    for (uint e = 0; e < elements.n; e++) {
      pe[e].part = lpart[e] + (c->id / chunk_size) * chunk_size;
      assert(pe[e].part < c->np);
    }
    free(lpart);

    sarray_transfer(struct element_t, &elements, proc, 0, cr);
    assert(nel == elements.n);
  }
  free(lvtx), free(lxyz);

  {
    sarray_sort(struct element_t, elements.ptr, elements.n, seq, 0, bfr);
    const struct element_t *const pe =
        (const struct element_t *const)elements.ptr;
    for (uint i = 0; i < nel; i++) part[i] = pe[i].part;
  }

  array_free(&elements);
}

static void update_frontier(sint *const target, sint *const hop,
                            sint *const frontier, const unsigned nv,
                            const unsigned hid, const struct comm *c,
                            buffer *const bfr) {
  // If target is already set, we don't update either target or hop.
  // We simply update frontier to previous target value and return.
  if (*target >= 0) {
    // Check invariant: *hop < INT_MAX
    assert(*hop < INT_MAX);
    for (uint i = 0; i < nv; i++) frontier[i] = *target;
    return;
  }

  struct dest_t {
    uint target;
  };

  struct array dests;
  array_init(struct dest_t, &dests, nv);
  {
    struct dest_t dt;
    for (uint i = 0; i < nv; i++) {
      if (frontier[i] >= 0) {
        dt.target = frontier[i];
        array_cat(struct dest_t, &dests, &dt, 1);
      }
    }
  }

  if (dests.n > 0) {
    sarray_sort(struct dest_t, dests.ptr, dests.n, target, 0, bfr);

    const struct dest_t *const pd = (const struct dest_t *const)dests.ptr;
    uint current_target = pd[0].target, current_count = 1;
    uint final_target = current_target, final_count = 1;
    for (uint i = 1; i < dests.n; i++) {
      if (pd[i].target == current_target) {
        current_count++;
      } else {
        if (current_count > final_count)
          final_count = current_count, final_target = current_target;
        current_target = pd[i].target, current_count = 1;
      }
    }
    if (current_count > final_count) final_target = current_target;

    // Update frontier, target and hop.
    for (uint j = 0; j < nv; j++) frontier[j] = final_target;
    *target = final_target, *hop = hid + 1;
  }

  array_free(&dests);
}

void parrsb_part_solid(int *part, const long long *const vtx2,
                       const unsigned nel2, const long long *const vtx1,
                       const unsigned nel1, const unsigned nv,
                       const MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);
  parrsb_print(&c, 1, "Running greedy solid ... nel1 = %d nel2 = %d", nel1,
               nel2);

  for (uint i = 0; i < nel2; i++) part[i] = -1;

  buffer bfr;
  buffer_init(&bfr, 1024);

  struct crystal cr;
  crystal_init(&cr, &c);

  // Return if global size is 0.
  const uint nelt = nel1 + nel2;
  slong nelg = nelt;
  {
    slong wrk;
    comm_allreduce(&c, gs_long, gs_add, &nelg, 1, &wrk);
    if (nelg == 0) {
      parrsb_print(&c, 1, "Mesh is empty ...");
      crystal_free(&cr);
      buffer_free(&bfr);
      comm_free(&c);
      return;
    }
  }

  const size_t size1 = nel1 * nv;
  const size_t size2 = nel2 * nv;
  const size_t size = size1 + size2;

  // Setup the gather-scatter handle to find connectivity through BFS.
  parrsb_print(&c, 1, "Setup gather-scatter handle ...");
  struct gs_data *gsh = NULL;
  {
    slong *vtx = tcalloc(slong, size);
    for (size_t i = 0; i < size1; i++) vtx[i] = vtx1[i];
    for (size_t i = 0; i < size2; i++) vtx[size1 + i] = vtx2[i];

    gsh = gs_setup(vtx, size, &c, 0, gs_pairwise, 0);
    free(vtx);
  }

  // Check if the solid + fluid mesh is connected. Otherwise, we cannot use
  // the greedy solid partitioner.
  parrsb_print(&c, 1, "Check if fluid + solid is connected ...");
  {
    slong wrk;
    sint idmin = (c.id + 1) * (size > 0);
    comm_allreduce(&c, gs_int, gs_min, &idmin, 1, &wrk);
    assert(idmin > 0);

    sint *const component = tcalloc(sint, size);
    if (c.id + 1 == (uint)idmin) {
      for (uint i = 0; i < nv; i++) component[i] = 1;
    }

    slong marked0 = 0, marked1 = 1;
    sint epoch = 0;
    while (marked1 > marked0) {
      gs(component, gs_int, gs_max, 0, gsh, &bfr);

      marked0 = marked1, marked1 = 0;
      for (uint i = 0; i < nel1 + nel2; i++) {
        sint v = 0;
        for (uint j = 0; j < nv; j++) v += component[i * nv + j];
        if (v > 0) {
          for (uint j = 0; j < nv; j++) component[i * nv + j] = 1;
          marked1 += 1;
        }
      }

      comm_allreduce(&c, gs_long, gs_add, &marked1, 1, &wrk);
      parrsb_print(&c, 1, "\tepoch = %d marked0 = %lld marked1 = %lld", epoch,
                   marked0, marked1);
      epoch++;
    }
    free(component);

    if (marked1 != nelg) {
      if (c.id == 0) {
        fprintf(stderr, "Fluid + Solid mesh is not connected.\n");
        fflush(stderr);
      }
      exit(EXIT_FAILURE);
    }
  }

  // Calculate the global number of elements in solid mesh and expected number
  // of elements in each partition.
  parrsb_print(&c, 1, "Calculate expected number of elements ...");
  slong nelgt2 = nel2;
  uint nexp2;
  {
    slong wrk;
    comm_allreduce(&c, gs_long, gs_add, &nelgt2, 1, &wrk);
    nexp2 = nelgt2 / c.np;
    nexp2 += (c.id < (nelgt2 - nexp2 * c.np));
    // Check for invariant: (min(nexp2) -  max(nexp2)) <= 1.
    slong nexp2_min = nexp2, nexp2_max = nexp2;
    comm_allreduce(&c, gs_long, gs_min, &nexp2_min, 1, &wrk);
    comm_allreduce(&c, gs_long, gs_max, &nexp2_max, 1, &wrk);
    assert(nexp2_max - nexp2_min <= 1);
    // Check for invariant: (sum(nexp2) == nelgt2).
    slong nexp2_sum = nexp2;
    comm_allreduce(&c, gs_long, gs_add, &nexp2_sum, 1, &wrk);
    assert(nexp2_sum == nelgt2);
  }

  // Initialize array of elements to be sent to each partition.
  struct elem_t {
    sint part;
    uint target, hop, sequence;
  };

  struct array arr;
  array_init(struct elem_t, &arr, nel2);

  // Allocate space for work arrays: frontier, target, and hop.
  sint *const frontier = tcalloc(sint, size);
  sint *const target = tcalloc(sint, nelt);
  sint *const hop = tcalloc(sint, nelt);

  uint nrecv2 = 0;
  slong nrem2 = nelgt2;
  while (nrem2 > 0) {
    parrsb_print(&c, 1, "nrem2 = %lld", nrem2);

    // Check for invariant: nrecv2 <= nexp2.
    assert(nrecv2 <= nexp2);

    // If the partition does not have enough elements, we keep it under
    // consideration for accepting new solid elements. If the partition
    // already has enough elements, we take that partition out of
    // consideration (by setting the frontier to -1). We always initialize solid
    // elements as unassigned (-1) although they may be already assigned. We
    // check for that later when we actually assign the elements to partitions.
    {
      sint id = c.id, hid = 0;
      if (nrecv2 == nexp2) id = -1, hid = INT_MAX;

      // Max id should be >= 0;
      sint wrk, idmax = id;
      comm_allreduce(&c, gs_int, gs_max, &idmax, 1, &wrk);
      assert(idmax >= 0);

      // Initialize frontier, target, and hop.
      for (uint i = 0; i < size1; i++) frontier[i] = id;
      for (uint i = size1; i < size; i++) frontier[i] = -1;
      for (uint i = 0; i < nel1; i++) target[i] = id, hop[i] = hid;
      for (uint i = nel1; i < nelt; i++) target[i] = -1, hop[i] = INT_MAX;
    }

    // Then perform a BFS till we assign all the elements in the solid mesh with
    // a potential partition id.
    parrsb_print(&c, 1, "Assign partition id ...");
    {
      sint assigned = 0;
      slong wrk;
      for (uint hid = 0; !assigned; hid++) {
        gs(frontier, gs_int, gs_max, 0, gsh, &bfr);

        assigned = 1;
        slong unassigned = 0;
        for (uint i = 0; i < nelt; i++) {
          update_frontier(&target[i], &hop[i], &frontier[i * nv], nv, hid, &c,
                          &bfr);
          assigned = assigned && (target[i] >= 0);
          unassigned += (target[i] < 0);
        }

        comm_allreduce(&c, gs_int, gs_min, &assigned, 1, &wrk);
        comm_allreduce(&c, gs_long, gs_add, &unassigned, 1, &wrk);
        parrsb_print(&c, 1, "hid = %d, assigned = %d unassigned = %d", hid,
                     assigned, unassigned);
      }
    }

    // Pack unassigned solid elements and send them to the target partition.
    arr.n = 0;
    {
      struct elem_t et = {.part = -1};
      for (uint i = 0; i < nel2; i++) {
        if (part[i] >= 0) continue;
        et.sequence = i, et.target = target[nel1 + i], et.hop = hop[nel1 + i];
        array_cat(struct elem_t, &arr, &et, 1);
      }

      parrsb_print(&c, 1, "Send elemenets to the target partition ...");
      sarray_transfer(struct elem_t, &arr, target, 1, &cr);
    }

    // Assign elements if the partition still doesn't have enough elements.
    if (nrecv2 < nexp2) {
      // We sort by hop value. Elements with lower hop value are assigned first
      // since they are technically closer to the partition.
      sarray_sort(struct elem_t, arr.ptr, arr.n, hop, 1, &bfr);
      struct elem_t *const pa = (struct elem_t *const)arr.ptr;
      uint keep = MIN(nexp2 - nrecv2, arr.n);
      for (uint i = 0; i < keep; i++) pa[i].part = c.id;
      nrecv2 += keep;
      // Check for invariant: nrecv2 <= nexp2.
      assert(nrecv2 <= nexp2);
    }

    // Send everything back with updated partition id and update the part array.
    {
      parrsb_print(&c, 1, "Send everything back ...");
      sarray_transfer(struct elem_t, &arr, target, 0, &cr);

      const struct elem_t *const pa = (const struct elem_t *const)arr.ptr;
      for (uint j = 0; j < arr.n; j++) part[pa[j].sequence] = pa[j].part;
      arr.n = 0;
    }

    // Update the number of elements remaining.
    {
      slong wrk;
      nrem2 = nexp2 - nrecv2;
      comm_allreduce(&c, gs_long, gs_add, &nrem2, 1, &wrk);
    }
  }

  gs_free(gsh);
  free(frontier), free(target), free(hop);
  array_free(&arr);
  crystal_free(&cr);
  buffer_free(&bfr);
  comm_free(&c);
}

int parrsb_part_mesh(int *part, const long long *const vtx,
                     const double *const xyz, const int *const tag,
                     const int nel, const int nv, parrsb_options *const options,
                     MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  update_options(options);

  // Check verboity and print a message.
  const int verbose = options->verbose_level;
  {
    slong nelg = nel, wrk;
    comm_allreduce(&c, gs_long, gs_add, &nelg, 1, &wrk);
    parrsb_print(&c, verbose, "Running parRSB ..., nv = %d, nelg = %lld", nv,
                 nelg);
  }

  print_options(&c, options);

  if (options->tagged == 1 && !tag) {
    parrsb_print(&c, verbose,
                 "Tagged partitioning requested but tag array is NULL..");
    return 1;
  }

  buffer bfr;
  buffer_init(&bfr, (nel + 1) * 72);

  struct crystal cr;
  crystal_init(&cr, &c);

  metric_init();

  parrsb_barrier(&c);
  const double t = comm_time();

  if (options->tagged == 1)
    parrsb_part_mesh_v1(part, vtx, xyz, tag, nel, nv, options, &c, &cr, &bfr);
  if (options->tagged == 0)
    parrsb_part_mesh_v0(part, vtx, xyz, nel, nv, options, &c, &cr, &bfr);

  parrsb_print(&c, verbose, "par%s finished in %g seconds.",
               ALGO[options->partitioner], comm_time() - t);

  metric_rsb_print(&c, options->profile_level);
  metric_finalize();

  crystal_free(&cr);
  buffer_free(&bfr);
  comm_free(&c);

  return 0;
}

#undef MIN

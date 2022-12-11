#include "metrics.h"
#include "parrsb-impl.h"
#include "sort.h"

static unsigned disconnected = 0;

extern int fiedler(struct array *elements, int nv, parrsb_options *options,
                   struct comm *gsc, buffer *buf, int verbose);

static void test_component_versions(struct array *elements, struct comm *lc,
                                    unsigned nv, unsigned lvl, buffer *bfr) {
  // Send elements to % P processor to test disconnected components
  struct crystal cr;
  crystal_init(&cr, lc);

  struct rsb_element *pe = (struct rsb_element *)elements->ptr;
  for (unsigned e = 0; e < elements->n; e++)
    pe[e].proc = pe[e].globalId % lc->np;

  sarray_transfer(struct rsb_element, elements, proc, 1, &cr);

  MPI_Comm tmp;
  int color = (lc->id < lc->np / 2);
  MPI_Comm_split(lc->c, color, lc->id, &tmp);

  struct comm tc0;
  comm_init(&tc0, tmp);

  sint nc1 = get_components(NULL, elements, nv, &tc0, bfr, 0);
  sint nc2 = get_components_v2(NULL, elements, nv, &tc0, bfr, 0);
  if (nc1 != nc2) {
    if (tc0.id == 0)
      printf("lvl = %u SS BFS != MS BFS: %d %d\n", lvl, nc1, nc2);
    fflush(stdout);
  }
  if (nc1 > 1) {
    if (tc0.id == 0)
      printf("lvl = %u: %d disconnected componets were present.\n", lvl, nc1);
    fflush(stdout);
  }

  comm_free(&tc0);
  MPI_Comm_free(&tmp);

  sarray_transfer(struct rsb_element, elements, proc, 0, &cr);
  crystal_free(&cr);
}

static void check_rsb_partition(struct comm *gc, parrsb_options *opts) {
  int max_levels = log2ll(gc->np);
  int miter = opts->rsb_max_iter, mpass = opts->rsb_max_passes;

  for (int i = 0; i < max_levels; i++) {
    sint converged = 1;
    int val = (int)metric_get_value(i, RSB_FIEDLER_CALC_NITER);
    if (opts->rsb_algo == 0) {
      if (val == miter * mpass)
        converged = 0;
    } else if (opts->rsb_algo == 1) {
      if (val == mpass)
        converged = 0;
    }

    sint ibfr;
    double dbfr;
    comm_allreduce(gc, gs_int, gs_min, &converged, 1, &ibfr);
    if (converged == 0) {
      if (opts->rsb_algo == 0) {
        double final = metric_get_value(i, TOL_FNL);
        comm_allreduce(gc, gs_double, gs_min, &final, 1, &dbfr);

        double target = metric_get_value(i, TOL_TGT);
        comm_allreduce(gc, gs_double, gs_min, &target, 1, &dbfr);

        if (gc->id == 0) {
          printf("Warning: Lanczos only reached a tolerance of %lf given %lf "
                 "after %d x %d iterations in Level=%d!\n",
                 final, target, mpass, miter, i);
          fflush(stdout);
        }
      } else if (opts->rsb_algo == 1) {
        if (gc->id == 0) {
          printf("Warning: Inverse iteration didn't converge after %d "
                 "iterations in Level = %d\n",
                 mpass, i);
          fflush(stdout);
        }
      }
    }

    sint minc, maxc;
    minc = maxc = (sint)metric_get_value(i, RSB_COMPONENTS);
    comm_allreduce(gc, gs_int, gs_min, &minc, 1, &ibfr);
    comm_allreduce(gc, gs_int, gs_max, &maxc, 1, &ibfr);

    if (maxc > 1 && gc->id == 0) {
      printf("Warning: Partition created %d/%d (min/max) disconnected "
             "components in Level=%d!\n",
             minc, maxc, i);
      fflush(stdout);
    }
  }
}

static int check_bin_val(int bin, struct comm *c) {
  if (bin < 0 || bin > 1) {
    if (c->id == 0) {
      printf("%s:%d bin value out of range: %d\n", __FILE__, __LINE__, bin);
      fflush(stdout);
    }
    return 1;
  }
  return 0;
}

int balance_partitions(struct array *elements, int nv, struct comm *lc,
                       struct comm *gc, int bin, buffer *bfr) {
  assert(check_bin_val(bin, gc) == 0);

  struct ielem_t {
    uint index, orig;
    sint dest;
    scalar fiedler;
  };

  // Calculate expected # of elements per processor
  uint ne = elements->n;
  slong nelgt = ne, nglob = ne, wrk;
  comm_allreduce(lc, gs_long, gs_add, &nelgt, 1, &wrk);
  comm_allreduce(gc, gs_long, gs_add, &nglob, 1, &wrk);

  sint ne_ = nglob / gc->np, nrem = nglob - ne_ * gc->np;
  slong nelgt_exp = ne_ * lc->np + nrem / 2 + (nrem % 2) * (1 - bin);
  slong send_cnt = nelgt - nelgt_exp > 0 ? nelgt - nelgt_exp : 0;

  // Setup gather-scatter
  uint size = ne * nv, e, v;
  slong *ids = tcalloc(slong, size);
  struct rsb_element *elems = (struct rsb_element *)elements->ptr;
  for (e = 0; e < ne; e++) {
    for (v = 0; v < nv; v++)
      ids[e * nv + v] = elems[e].vertices[v];
  }
  struct gs_data *gsh = gs_setup(ids, size, gc, 0, gs_pairwise, 0);

  sint *input = (sint *)ids;
  if (send_cnt > 0)
    for (e = 0; e < size; e++)
      input[e] = 0;
  else
    for (e = 0; e < size; e++)
      input[e] = 1;

  gs(input, gs_int, gs_add, 0, gsh, bfr);

  for (e = 0; e < ne; e++)
    elems[e].proc = gc->id;

  sint sid = (send_cnt == 0) ? gc->id : INT_MAX, balanced = 0;
  comm_allreduce(gc, gs_int, gs_min, &sid, 1, &wrk);

  struct crystal cr;

  if (send_cnt > 0) {
    struct array ielems;
    array_init(struct ielem_t, &ielems, 10);

    struct ielem_t ielem = {
        .index = 0, .orig = lc->id, .dest = -1, .fiedler = 0};
    int mul = (sid == 0) ? 1 : -1;
    for (e = 0; e < ne; e++) {
      for (v = 0; v < nv; v++) {
        if (input[e * nv + v] > 0) {
          ielem.index = e, ielem.fiedler = mul * elems[e].fiedler;
          array_cat(struct ielem_t, &ielems, &ielem, 1);
          break;
        }
      }
    }

    // Sort based on fiedler value and sets `orig` field
    parallel_sort(struct ielem_t, &ielems, fiedler, gs_double, 0, 1, lc, bfr);

    slong out[2][1], bfr[2][1], nielems = ielems.n;
    comm_scan(out, lc, gs_long, gs_add, &nielems, 1, bfr);
    slong start = out[0][0];

    sint P = gc->np - lc->np;
    sint part_size = (send_cnt + P - 1) / P;

    if (out[1][0] >= send_cnt) {
      balanced = 1;
      struct ielem_t *ptr = ielems.ptr;
      for (e = 0; start + e < send_cnt && e < ielems.n; e++)
        ptr[e].dest = sid + (start + e) / part_size;

      crystal_init(&cr, lc);
      sarray_transfer(struct ielem_t, &ielems, orig, 0, &cr);
      crystal_free(&cr);

      ptr = ielems.ptr;
      for (e = 0; e < ielems.n; e++)
        if (ptr[e].dest != -1)
          elems[ptr[e].index].proc = ptr[e].dest;
    }

    array_free(&ielems);
  }

  comm_allreduce(gc, gs_int, gs_max, &balanced, 1, &wrk);
  if (balanced == 1) {
    crystal_init(&cr, gc);
    sarray_transfer(struct rsb_element, elements, proc, 0, &cr);
    crystal_free(&cr);

    // Do a load balanced sort in each partition
    parallel_sort(struct rsb_element, elements, fiedler, gs_double, 0, 1, lc,
                  bfr);
  } else {
    // Forget about disconnected components, just do a load balanced partition
    // TODO: Need to change how parallel_sort load balance
    parallel_sort(struct rsb_element, elements, fiedler, gs_double, 0, 1, gc,
                  bfr);
  }

  free(ids);
  gs_free(gsh);
}

int repair_partitions_v2(struct array *elems, unsigned nv, struct comm *tc,
                         struct comm *lc, unsigned bin, unsigned algo,
                         buffer *bfr) {
  assert(check_bin_val(bin, lc) == 0);

  sint ibuf;
  sint nc = get_components_v2(NULL, elems, nv, tc, bfr, 0);
  comm_allreduce(lc, gs_int, gs_max, &nc, 1, &ibuf);
  if (nc > 1) {
    // If nc > 1, send elements back and do RCBx, RCBy and RCBz
    struct crystal cr;
    crystal_init(&cr, lc);
    sarray_transfer(struct rsb_element, elems, proc, 0, &cr);
    crystal_free(&cr);

    // Do rcb or rib
    unsigned ndim = (nv == 8) ? 3 : 2;
    switch (algo) {
    case 0:
      parallel_sort(struct rsb_element, elems, globalId, gs_long, 0, 1, lc,
                    bfr);
      break;
    case 1:
      rcb(elems, sizeof(struct rsb_element), ndim, lc, bfr);
      break;
    case 2:
      rib(elems, sizeof(struct rsb_element), ndim, lc, bfr);
      break;
    default:
      break;
    }

    // And count number of components again. If nc > 1 still, set
    // isconnected = 1
    nc = get_components_v2(NULL, elems, nv, tc, bfr, 0);
    comm_allreduce(lc, gs_int, gs_max, &nc, 1, &ibuf);
    if (nc > 1)
      disconnected = 1;
  }

  return 0;
}

static void get_part(sint *np, sint *nid, int two_lvl, struct comm *lc,
                     struct comm *nc) {
  if (two_lvl) {
    sint out[2][1], wrk[2][1], in = (nc->id == 0);
    comm_scan(out, lc, gs_int, gs_add, &in, 1, &wrk);
    *nid = (nc->id == 0) * out[0][0], *np = out[1][0];
    comm_allreduce(nc, gs_int, gs_max, nid, 1, wrk);
  } else {
    *np = lc->np, *nid = lc->id;
  }
}

int rsb(struct array *elements, int nv, int check, parrsb_options *options,
        struct comm *gc, buffer *bfr) {
  // `gc` is the global communicator. We make a duplicate of it in `lc` and
  // keep splitting it. `nc` is the communicator for the two level partitioning.
  struct comm lc, nc;

  // Duplicate the global communicator to `lc`
  comm_dup(&lc, gc);

  // Initialize `nc` based on `lc`
  if (options->two_level) {
#ifdef MPI
    MPI_Comm node;
    MPI_Comm_split_type(lc.c, MPI_COMM_TYPE_SHARED, lc.id, MPI_INFO_NULL,
                        &node);
    comm_init(&nc, node);
    MPI_Comm_free(&node);
#else
    comm_init(&nc, 1);
#endif
  }

  // Get number of partitions we are going to perform RSB on first level
  sint np, nid;
  get_part(&np, &nid, options->two_level, &lc, &nc);
  if (options->two_level && options->verbose_level) {
    if (gc->id == 0)
      printf("Number of nodes = %d\n", np);
    fflush(stdout);
  }

  unsigned ndim = (nv == 8) ? 3 : 2;
  while (np > 1) {
    // Run the pre-partitioner
    metric_tic(&lc, RSB_PRE);
    switch (options->rsb_pre) {
    case 0: // Sort by global id
      parallel_sort(struct rsb_element, elements, globalId, gs_long, 0, 1, &lc,
                    bfr);
      break;
    case 1: // RCB
      rcb(elements, sizeof(struct rsb_element), ndim, &lc, bfr);
      break;
    case 2: // RIB
      rib(elements, sizeof(struct rsb_element), ndim, &lc, bfr);
      break;
    default:
      break;
    }
    metric_toc(&lc, RSB_PRE);

    // Find the Fiedler vector
    unsigned bin = (nid >= (np + 1) / 2);
    struct comm tc;
    comm_split(&lc, bin, lc.id, &tc);

    struct rsb_element *pe = (struct rsb_element *)elements->ptr;
    for (unsigned i = 0; i < elements->n; i++)
      pe[i].proc = lc.id;

    metric_tic(&lc, RSB_FIEDLER);
    fiedler(elements, nv, options, &lc, bfr, gc->id == 0);
    metric_toc(&lc, RSB_FIEDLER);

    // Sort by Fiedler vector
    metric_tic(&lc, RSB_SORT);
    parallel_sort_2(struct rsb_element, elements, fiedler, gs_double, globalId,
                    gs_long, 0, 1, &lc, bfr);
    metric_toc(&lc, RSB_SORT);

    // Attempt to repair if there are disconnected components
    metric_tic(&lc, RSB_REPAIR);
    if (options->repair)
      repair_partitions_v2(elements, nv, &tc, &lc, bin, options->rsb_pre, bfr);
    metric_toc(&lc, RSB_REPAIR);

    // Bisect and balance
    metric_tic(&lc, RSB_BALANCE);
    balance_partitions(elements, nv, &tc, &lc, bin, bfr);
    metric_toc(&lc, RSB_BALANCE);

    // Split the communicator
    comm_free(&lc);
    comm_dup(&lc, &tc);
    comm_free(&tc);

    get_part(&np, &nid, options->two_level, &lc, &nc);
    metric_push_level();
  }
  comm_free(&lc);

  // Partition within the node
  if (options->two_level) {
    options->two_level = 0;
    rsb(elements, nv, 0, options, &nc, bfr);
    comm_free(&nc);
  }

  if (check)
    check_rsb_partition(gc, options);

  return 0;
}

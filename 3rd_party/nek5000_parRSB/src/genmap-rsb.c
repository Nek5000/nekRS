#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#include <genmap-impl.h>
#include <sort.h>

static int check_convergence(struct comm *gc, int max_pass, int max_iter) {
  int max_levels = log2(gc->np);

  int i;
  for (i = 0; i < max_levels; i++) {
    sint converged = 1;
    int val = (int)metric_get_value(i, NFIEDLER);
    if (val >= max_pass * max_iter) {
      converged = 0;
      break;
    }

    sint ibfr;
    comm_allreduce(gc, gs_int, gs_min, &converged, 1, &ibfr);

    double dbfr;
    double final = (double)metric_get_value(i, LANCZOSTOLFINAL);
    comm_allreduce(gc, gs_double, gs_min, &final, 1, &dbfr);

    double target = (double)metric_get_value(i, LANCZOSTOLTARGET);

    if (converged == 0 && gc->id == 0) {
      printf("\tWarning: Partitioner only reached a tolerance of %lf given %lf "
             "after %d x %d iterations in Level=%d!\n",
             final, target, max_pass, max_iter, i);
      fflush(stdout);
    }
  }
}

int genmap_rsb(genmap_handle h) {
  int verbose = h->options->debug_level > 1;
  int max_iter = 50;
  int max_pass = 50;

  struct comm *lc = h->local;
  struct comm *gc = h->global;

  genmap_comm_scan(h, lc);

  uint nelt = genmap_get_nel(h);
  struct rsb_element *e = genmap_get_elements(h);
  GenmapInt i;
  for (i = 0; i < nelt; i++)
    e[i].globalId0 = genmap_get_local_start_index(h) + i + 1;

  int nv = h->nv;
  int ndim = (nv == 8) ? 3 : 2;

  int np = gc->np;
  int level = 0;

  while ((np = lc->np) > 1) {
    int global = 1;

    /* Run RCB, RIB pre-step or just sort by global id */
    if (h->options->rsb_prepartition == 1) { // RCB
      metric_tic(lc, RCB);
      rcb(lc, h->elements, ndim, &h->buf);
      metric_toc(lc, RCB);
    } else if (h->options->rsb_prepartition == 2) { // RIB
      metric_tic(lc, RCB);
      rib(lc, h->elements, ndim, &h->buf);
      metric_toc(lc, RCB);
    } else {
      parallel_sort(struct rsb_element, h->elements, globalId0, gs_long, 0, 1,
                    lc, &h->buf);
    }

    /* Initialize the laplacian */
    metric_tic(lc, WEIGHTEDLAPLACIANSETUP);
    GenmapInitLaplacianWeighted(h, lc);
    metric_toc(lc, WEIGHTEDLAPLACIANSETUP);

    /* Run fiedler */
    metric_tic(lc, FIEDLER);
    int ipass = 0, iter;
    do {
      if (h->options->rsb_algo == 0)
        iter = GenmapFiedlerLanczos(h, lc, max_iter, global);
      else if (h->options->rsb_algo == 1)
        iter = GenmapFiedlerRQI(h, lc, max_iter, global);
      metric_acc(NFIEDLER, iter);
      global = 0;
    } while (++ipass < max_pass && iter == max_iter);
    metric_toc(lc, FIEDLER);

    /* Sort by Fiedler vector */
    metric_tic(lc, FIEDLERSORT);
    parallel_sort(struct rsb_element, h->elements, fiedler, gs_double, 0, 1, lc,
                  &h->buf);
    metric_toc(lc, FIEDLERSORT);

    /* Bisect */
    double t = comm_time();
    split_and_repair_partitions(h, lc, level, gc);
    t = comm_time() - t;
    metric_acc(BISECTANDREPAIR, t);

    genmap_comm_scan(h, lc);
    metric_push_level();
    level++;
  }

  check_convergence(gc, max_pass, max_iter);

  return 0;
}

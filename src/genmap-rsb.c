#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#include <genmap-impl.h>
#include <sort.h>

static int dump_fiedler_if_discon(genmap_handle h, int level, int max_levels) {
  struct comm *lc = h->local;
  struct comm *gc = h->global;

  int nv = h->nv;
  int ndim = (nv == 8) ? 3 : 2;

  sint bfr[2];

  struct rsb_element *elements = genmap_get_elements(h);

  /* Dump current partition status */
  if (level > 0 && level < max_levels) {
    slong nelt = genmap_get_nel(h);
    slong out[2][1], buf[2][1];
    comm_scan(out, gc, gs_long, gs_add, &nelt, 1, buf); // max
    slong start = out[0][0];

    sint components = get_components(NULL, elements, lc, &h->buf, nelt, nv);
    comm_allreduce(lc, gs_int, gs_max, &components, 1, bfr); // max

    sint g_id = (components > 1) * gc->id;
    comm_allreduce(gc, gs_int, gs_max, &g_id, 1, bfr); // max

    sint l_id = gc->id;
    comm_allreduce(lc, gs_int, gs_max, &l_id, 1, bfr); // max

    if (g_id == l_id && components > 1) {
      if (lc->id == 0)
        printf("\tLevel %02d PRERCB: There are disconnected components!\n",
               level);
      if (components > 1) {
        // Dump the current partition
        char fname[BUFSIZ];
        sprintf(fname, "fiedler_%02d.dump", level);
        GenmapFiedlerDump(fname, h, start, lc);
      }
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
  int level = 0, max_levels = log2(np);

  sint bfr[2];

  while ((np = lc->np) > 1) {
    int global;
    if (h->options->rsb_paul == 1)
      global = 1;
    else
      global = (np == gc->np);

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
    split_and_repair_partitions(h, lc, level);
    t = comm_time() - t;
    metric_acc(BISECTANDREPAIR, t);

    genmap_comm_scan(h, lc);
    metric_push_level();
    level++;
  }

  /* Check if Fidler converged */
  sint converged = 1;
  for (i = 0; i < metric_get_levels(); i++) {
    int val = (int)metric_get_value(i, NFIEDLER);
    if (val >= max_pass * max_iter) {
      converged = 0;
      level = i;
      break;
    }
  }
  comm_allreduce(gc, gs_int, gs_min, &converged, 1, bfr); // min
  if (converged == 0 && gc->id == 0)
    printf("\tWARNING: Failed to converge while partitioning, Level=%d!\n",
           level);

  return 0;
}

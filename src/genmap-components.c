#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#include <genmap-impl.h>
#include <sort.h>

struct unmarked {
  uint index;
};

/* Find the number of disconnected components */
sint get_components(sint *component, struct rsb_element *elements,
                    struct comm *c, buffer *buf, uint nelt, uint nv) {
  slong nelt_ = nelt;
  slong out[2][1], buff[2][1];
  comm_scan(out, c, gs_long, gs_add, &nelt_, 1, buff);
  slong nelg = out[1][0];
  ulong start = out[0][0];

  if (nelg == 0)
    return 0;

  GenmapLong *p, *ids;
  GenmapMalloc(nelt * nv, &p);
  GenmapMalloc(nelt * nv, &ids);

  int null_input = 0;
  if (component == NULL) {
    GenmapMalloc(nelt, &component);
    null_input = 1;
  }

  uint e;
  for (e = 0; e < nelt; e++)
    component[e] = -1;

  struct array arr;
  array_init(struct unmarked, &arr, nelt);

  struct comm cc;

  struct unmarked u;
  uint d;
  slong nnz1, nnzg, nnzg0, nnzb, nmarked = 0;
  sint count = 0;

  do {
    // Count unmarked elements
    arr.n = 0;
    for (e = 0; e < nelt; e++) {
      if (component[e] == -1) {
        u.index = e;
        array_cat(struct unmarked, &arr, &u, 1);
      }
    }

    int bin = 0;
    if (arr.n > 0)
      bin = 1;

    genmap_comm_split(c, bin, c->id, &cc);

    nnz1 = nnzg = nnzg0 = 0;

    if (bin == 1) {
      // Initialize p
      for (e = 0; e < arr.n; e++)
        for (d = 0; d < nv; d++)
          p[e * nv + d] = 0;

      // Mark the first non-marked element as seed
      struct unmarked *ptr = (struct unmarked *)arr.ptr;
      slong first = start + ptr[0].index;
      slong first_ = first;
      comm_allreduce(&cc, gs_long, gs_min, &first_, 1, buff);

      if (first_ == first) {
        for (d = 0; d < nv; d++)
          p[0 * nv + d] = 1;
      }

      // Setup gs
      for (e = 0; e < arr.n; e++)
        for (d = 0; d < nv; d++)
          ids[e * nv + d] = elements[ptr[e].index].vertices[d];

      struct gs_data *gsh = gs_setup(ids, arr.n * nv, &cc, 0, gs_pairwise, 0);

      do {
        gs(p, gs_long, gs_add, 0, gsh, buf);

        nnz1 = 0;
        for (e = 0; e < arr.n; e++) {
          for (d = 0; d < nv; d++) {
            if (p[e * nv + d] > 0) {
              nnz1++;
              component[ptr[e].index] = count;
              break;
            }
          }

          if (d < nv) { // There was one non-zero vertex in the element
            for (d = 0; d < nv; d++)
              p[e * nv + d] = 1;
          }
        }

        nnzg0 = nnzg, nnzg = nnz1;
        comm_allreduce(&cc, gs_long, gs_add, &nnzg, 1, &nnzb);
      } while (nnzg > nnzg0);

      gs_free(gsh);
    }

    comm_free(&cc);

    comm_allreduce(c, gs_long, gs_add, &nnz1, 1, &nnzb);
    nmarked += nnz1;

    count++;
  } while (nmarked < nelg);

  if (null_input == 1)
    GenmapFree(component);

  GenmapFree(p);
  GenmapFree(ids);

  return count;
}

void split_and_repair_partitions(genmap_handle h, struct comm *lc, int level) {
  sint np = lc->np;
  int bin = 1;
  if (lc->id < (np + 1) / 2)
    bin = 0;
  struct comm tc;
  genmap_comm_split(lc, bin, lc->id, &tc);

  struct rsb_element *e = genmap_get_elements(h);
  uint nelt = genmap_get_nel(h);
  int nv = genmap_get_nvertices(h);

  /* Check for disconnected components */
  GenmapInitLaplacianWeighted(h, &tc);

  sint *comp_ids = NULL;
  GenmapMalloc(nelt, &comp_ids);
  sint ncomp_global, ncomp, buf;
  ncomp_global = ncomp = get_components(comp_ids, e, &tc, &h->buf, nelt, nv);
  comm_allreduce(lc, gs_int, gs_max, &ncomp_global, 1, &buf);

  while (ncomp_global > 1) {
    if (lc->id == 0) {
      printf("\tWarning: There are %d disconnected components in level = %d!\n",
             ncomp_global, level);
      fflush(stdout);
    }

    sint *comp_count = NULL;
    GenmapCalloc(2 * ncomp_global, &comp_count);
    uint i;
    for (i = 0; i < nelt; i++)
      comp_count[comp_ids[i]]++;
    comm_allreduce(&tc, gs_int, gs_add, comp_count, ncomp_global,
                   &comp_count[ncomp_global]);

    sint min_count = INT_MAX, min_id = -1;
    for (i = 0; i < ncomp; i++) {
      if (comp_count[i] < min_count) {
        min_count = comp_count[i];
        min_id = i;
      }
    }
    sint min_count_global = min_count;
    comm_allreduce(lc, gs_int, gs_min, &min_count_global, 1, &buf);
    sint id_global = (min_count_global == min_count) ? lc->id : lc->np;
    comm_allreduce(lc, gs_int, gs_min, &id_global, 1, &buf);

    struct crystal cr;
    crystal_init(&cr, lc);

    for (i = 0; i < nelt; i++)
      e[i].proc = lc->id;

    sint low_np = (np + 1) / 2;
    sint high_np = np - low_np;
    sint start = !bin * low_np;
    sint P = bin * low_np + !bin * high_np;
    sint size = (min_count_global + P - 1) / P;

    sint current = 0;
    if (min_count_global == min_count && lc->id == id_global) {
      for (i = 0; i < nelt; i++) {
        if (comp_ids[i] == min_id) {
          e[i].proc = start + current / size;
          current++;
        }
      }
      assert(min_count == current && "min_count != current");
    }

    sarray_transfer(struct rsb_element, h->elements, proc, 1, &cr);
    crystal_free(&cr);

    // do a load balanced sort in each partition
    parallel_sort(struct rsb_element, h->elements, fiedler, gs_double, 0, 1,
                  &tc, &h->buf);

    genmap_comm_scan(h, &tc);

    e = genmap_get_elements(h);
    nelt = genmap_get_nel(h);
    GenmapRealloc(nelt, &comp_ids);
    GenmapInitLaplacianWeighted(h, &tc);
    ncomp_global = ncomp = get_components(comp_ids, e, &tc, &h->buf, nelt, nv);
    comm_allreduce(lc, gs_int, gs_max, &ncomp_global, 1, &buf);

    GenmapFree(comp_count);
  }

  GenmapFree(comp_ids);
  comm_free(lc);
  comm_dup(lc, &tc);
  comm_free(&tc);
}

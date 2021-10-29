#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#include <genmap-impl.h>
#include <sort.h>

struct unmarked {
  uint index;
};

struct interface_element {
  uint index;
  uint orig;
  sint dest;
  GenmapScalar fiedler;
};

/* Check the bin value */
int check_bin_val(int bin, struct comm *gc) {
  if (bin < 0 || bin > 1) {
    if (gc->id == 0) {
      printf("%s:%d bin value out of range: %d\n", __FILE__, __LINE__, bin);
      fflush(stdout);
    }
    return 1;
  }
  return 0;
}

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
  slong nnz1, nnzg, nnzg0, nnzb;
  slong nmarked = 0;
  sint count = 0;

  do {
    /* Count unmarked elements */
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
      /* Initialize p */
      for (e = 0; e < arr.n; e++)
        for (d = 0; d < nv; d++)
          p[e * nv + d] = 0;

      /* Mark the first non-marked element as seed */
      struct unmarked *ptr = (struct unmarked *)arr.ptr;
      slong first = start + ptr[0].index;
      slong first_ = first;
      comm_allreduce(&cc, gs_long, gs_min, &first_, 1, buff);

      if (first_ == first) {
        for (d = 0; d < nv; d++)
          p[0 * nv + d] = 1;
      }

      /* Setup gs */
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

  GenmapFree(p);
  GenmapFree(ids);
  if (null_input == 1)
    GenmapFree(component);

  return count;
}

int balance_partitions(genmap_handle h, struct comm *lc, int bin,
                       struct comm *gc) {
  assert(check_bin_val(bin, gc) == 0);

  uint nelt = genmap_get_nel(h);
  slong buf;
  slong nelgt = nelt;
  comm_allreduce(lc, gs_long, gs_add, &nelgt, 1, &buf);

  /* Calculate expected # of elements per processor */
  slong nglob = nelt;
  comm_allreduce(gc, gs_long, gs_add, &nglob, 1, &buf);

  slong nelt_ = nglob / gc->np;
  slong nelgt_exp = nelt_ * lc->np;
  sint nrem = nglob - nelt_ * gc->np;
  nelgt_exp += nrem / 2 + (nrem % 2) * (1 - bin);

  slong send_cnt = 0;
  if (nelgt - nelgt_exp > 0)
    send_cnt = nelgt - nelgt_exp;

  /* Setup gather-scatter */
  int nv = genmap_get_nvertices(h);
  uint size = nelt * nv;

  slong *ids = NULL;
  GenmapMalloc(size, &ids);

  struct rsb_element *elems = genmap_get_elements(h);
  uint e, v;
  for (e = 0; e < nelt; e++)
    for (v = 0; v < nv; v++)
      ids[e * nv + v] = elems[e].vertices[v];

  struct gs_data *gsh = gs_setup(ids, size, gc, 0, gs_pairwise, 0);

  sint *input = NULL;
  GenmapMalloc(size, &input);

  if (send_cnt > 0)
    for (e = 0; e < size; e++)
      input[e] = 0;
  else
    for (e = 0; e < size; e++)
      input[e] = 1;

  gs(input, gs_int, gs_add, 0, gsh, &h->buf);

  for (e = 0; e < nelt; e++)
    elems[e].proc = gc->id;

  sint start_id = (send_cnt == 0) ? gc->id : INT_MAX;
  comm_allreduce(gc, gs_int, gs_min, &start_id, 1, &buf);

  struct crystal cr;
  sint balanced = 0;

  if (send_cnt > 0) {
    int mul = -1;
    if (start_id == 0) /* we are sending to lower fiedler values */
      mul = 1;

    struct array ielems;
    array_init(struct interface_element, &ielems, 10);

    struct interface_element ielem;
    ielem.dest = -1;

    for (e = 0; e < nelt; e++) {
      for (v = 0; v < nv; v++)
        if (input[e * nv + v] > 0) {
          ielem.index = e;
          ielem.orig = lc->id;
          ielem.fiedler = mul * elems[e].fiedler;
          array_cat(struct interface_element, &ielems, &ielem, 1);
          break;
        }
    }

    parallel_sort(struct interface_element, &ielems, fiedler, gs_double, 0, 1,
                  lc, &h->buf);

    slong ielems_n = ielems.n;
    slong out[2][1], bfr[2][1];
    comm_scan(out, lc, gs_long, gs_add, &ielems_n, 1, bfr);
    slong start = out[0][0];

    sint P = gc->np - lc->np;
    slong part_size = (send_cnt + P - 1) / P;

    if (out[1][0] < send_cnt)
      balanced = 0;
    else {
      struct interface_element *ptr = ielems.ptr;
      for (e = 0; start + e < send_cnt && e < ielems.n; e++)
        ptr[e].dest = start_id + (start + e) / part_size;

      crystal_init(&cr, lc);
      sarray_transfer(struct interface_element, &ielems, orig, 0, &cr);
      crystal_free(&cr);

      ptr = ielems.ptr;
      for (e = 0; e < ielems.n; e++)
        if (ptr[e].dest != -1)
          elems[ptr[e].index].proc = ptr[e].dest;
    }

    array_free(&ielems);
  }

  comm_allreduce(gc, gs_int, gs_min, &balanced, 1, &buf);
  if (balanced == 1) {
    crystal_init(&cr, gc);
    sarray_transfer(struct rsb_element, h->elements, proc, 0, &cr);
    crystal_free(&cr);

    /* Do a load balanced sort in each partition */
    parallel_sort(struct rsb_element, h->elements, fiedler, gs_double, 0, 1, lc,
                  &h->buf);
  } else {
    /* Forget repair, just do a load balanced partition */
    parallel_sort(struct rsb_element, h->elements, fiedler, gs_double, 0, 1, gc,
                  &h->buf);
  }

  genmap_comm_scan(h, lc);

  nelt = genmap_get_nel(h);
  sint ncomp =
      get_components(NULL, genmap_get_elements(h), lc, &h->buf, nelt, nv);
  metric_acc(COMPONENTS, ncomp);

  GenmapFree(input);
  gs_free(gsh);
  GenmapFree(ids);
}

int repair_partitions(genmap_handle h, struct comm *tc, struct comm *lc,
                      int bin, struct comm *gc) {
  assert(check_bin_val(bin, gc) == 0);

  slong buf;
  uint nelt = genmap_get_nel(h);
  slong nelg = nelt;
  comm_allreduce(lc, gs_long, gs_add, &nelg, 1, &buf);

  /* Check for disconnected components */
  GenmapInitLaplacianWeighted(h, tc);
  sint np = lc->np;
  sint id = lc->id;

  struct rsb_element *e = genmap_get_elements(h);
  uint nv = genmap_get_nvertices(h);

  sint *comp_ids = NULL;
  GenmapMalloc(nelt, &comp_ids);

  slong min, max;
  sint ncomp = get_components(comp_ids, e, tc, &h->buf, nelt, nv);
  slong ncompg = ncomp;
  comm_allreduce(lc, gs_long, gs_max, &ncompg, 1, &buf);

  sint root = (lc->id == 0) * gc->id;
  comm_allreduce(lc, gs_int, gs_max, &root, 1, &buf);

  int attempt = 0;
  int nattempts = 1;

  while (ncompg > 1 && attempt < nattempts) {
    slong *comp_count = NULL;
    GenmapCalloc(3 * ncomp, &comp_count);

    uint i;
    for (i = 0; i < nelt; i++)
      comp_count[comp_ids[i]]++;

    for (i = 0; i < ncomp; i++)
      comp_count[ncomp + i] = comp_count[i];

    comm_allreduce(tc, gs_long, gs_add, &comp_count[ncomp], ncomp,
                   &comp_count[2 * ncomp]);

    slong min_count = LONG_MAX;
    sint min_id = -1;
    for (i = 0; i < ncomp; i++) {
      if (comp_count[ncomp + i] < min_count) {
        min_count = comp_count[ncomp + i];
        min_id = i;
      }
    }

    slong min_count_global = min_count;
    comm_allreduce(lc, gs_long, gs_min, &min_count_global, 1, &buf);

    /* bin is the tie breaker */
    sint min_bin = (min_count_global == min_count) ? bin : INT_MAX;
    comm_allreduce(lc, gs_int, gs_min, &min_bin, 1, &buf);

    e = genmap_get_elements(h);
    for (i = 0; i < nelt; i++)
      e[i].proc = id;

    sint low_np = (np + 1) / 2;
    sint high_np = np - low_np;
    sint start = (1 - bin) * low_np;
    sint P = bin * low_np + (1 - bin) * high_np;
    slong size = (min_count_global + P - 1) / P;

    if (min_count_global == min_count && min_bin == bin) {
      slong in = comp_count[min_id];
      slong out[2][1], buff[2][1];
      comm_scan(out, tc, gs_long, gs_add, &in, 1, buff);
      slong off = out[0][0];

      for (i = 0; i < nelt; i++) {
        if (comp_ids[i] == min_id) {
          e[i].proc = start + off / size;
          off++;
        }
      }
    }

    struct crystal cr;
    crystal_init(&cr, lc);
    sarray_transfer(struct rsb_element, h->elements, proc, 0, &cr);
    crystal_free(&cr);

    attempt++;

    /* Do a load balanced sort in each partition */
    parallel_sort(struct rsb_element, h->elements, fiedler, gs_double, 0, 1, tc,
                  &h->buf);
    genmap_comm_scan(h, tc);
    nelt = genmap_get_nel(h);
    GenmapInitLaplacianWeighted(h, tc);

    GenmapRealloc(nelt, &comp_ids);
    ncompg = ncomp =
        get_components(comp_ids, genmap_get_elements(h), tc, &h->buf, nelt, nv);
    comm_allreduce(lc, gs_long, gs_max, &ncompg, 1, &buf);

    GenmapFree(comp_count);
  }

  GenmapFree(comp_ids);

  return 0;
}

#include "coarse-impl.h"
#include "metrics.h"
#include "sort.h"

//------------------------------------------------------------------------------
// Better API for coarse grid system.
//
uint unique_ids(sint *perm, ulong *uid, uint n, const ulong *ids, buffer *bfr) {
  struct id_t {
    ulong id;
    uint idx;
    sint perm;
  };

  struct array arr;
  array_init(struct id_t, &arr, n);

  uint i;
  struct id_t t = {.id = 0, .idx = 0, .perm = -1};
  for (i = 0; i < n; i++) {
    t.id = ids[i], t.idx = i;
    array_cat(struct id_t, &arr, &t, 1);
  }

  sarray_sort(struct id_t, arr.ptr, arr.n, id, 1, bfr);
  struct id_t *pa = (struct id_t *)arr.ptr;

  // Ignore the ids numbered zero
  sint un = 0;
  ulong last = 0;
  for (uint i = 0; i < arr.n; i++) {
    ulong v = pa[i].id;
    if (v != last)
      last = uid[un] = v, un++;
    pa[i].perm = un - 1;
  }

  sarray_sort(struct id_t, pa, n, idx, 0, bfr);
  pa = (struct id_t *)arr.ptr;
  for (i = 0; i < n; i++)
    perm[i] = pa[i].perm;

  array_free(&arr);
  return un;
}

// Number rows, local first then interface. Returns global number of local
// elements.
struct rsb_t {
  uint i, s;
  slong vtx[8];
};

static void number_dofs(slong *nid, struct coarse *crs, const slong *ids,
                        const ulong *uid) {
  uint un = crs->un;
  buffer *bfr = &crs->bfr;
  struct comm *ci = &crs->c;
  sint *u2c = crs->u2c;

  int nnz = (un > 0);
  struct comm c;
  comm_split(ci, nnz, ci->id, &c);

  uint i, j;
  if (nnz) {
    sint *dof = tcalloc(sint, un);
    int level = 1;
    while (c.np > 1) {
      struct gs_data *gsh = gs_setup(ids, un, &c, 0, gs_pairwise, 0);

      int bin = (c.id >= (c.np + 1) / 2);
      for (i = 0; i < un; i++)
        dof[i] = bin;

      gs(dof, gs_int, gs_add, 0, gsh, bfr);

      if (bin == 1) {
        for (i = 0; i < un; i++)
          dof[i] = 0;
      }

      gs(dof, gs_int, gs_add, 0, gsh, bfr);

      for (i = 0; i < un; i++) {
        if (dof[i] > 0 && u2c[i] >= 0 && !nid[u2c[i]])
          nid[u2c[i]] = level;
      }

      gs_free(gsh);

      struct comm t;
      comm_split(&c, bin, c.id, &t);
      comm_free(&c);
      comm_dup(&c, &t);
      comm_free(&t);

      level++;
    }
    free(dof);
  }

  // Calculate unqiue local and interface nodes based on compress ids.
  // Finding unique local ids is easy. To find unique interface ids, we
  // will have to sort in parallel and then manually find the unique ids.
  struct dof_t {
    ulong id, nid;
    uint p, p0, idx;
  };

  struct array arr;
  array_init(struct dof_t, &arr, crs->cn);

  uint ln = 0;
  struct dof_t t = {.id = 0, .nid = 0, .p = 0, .p0 = ci->id, .idx = 0};
  for (i = 0; i < crs->cn; i++) {
    if (!nid[i])
      ln++;
    else
      t.id = uid[i], t.idx = i, array_cat(struct dof_t, &arr, &t, 1);
  }
  crs->n[0] = ln;

  slong cnt[1] = {ln}, out[2][1], wrk[2][1];
  comm_scan(out, ci, gs_long, gs_add, cnt, 1, wrk);
  crs->s[0] = out[0][0] + 1, crs->ng[0] = out[1][0];

  for (i = 0, ln = 0; i < crs->cn; i++) {
    if (!nid[i])
      nid[i] = crs->s[0] + ln, ln++;
  }
  assert(crs->n[0] == ln);

  // parallel_sort and set nid and send back to p0
  parallel_sort(struct dof_t, &arr, id, gs_long, 0, 0, ci, bfr);

  uint in = 0;
  if (arr.n > 0) {
    struct dof_t *pa = (struct dof_t *)arr.ptr;
    for (i = in = 1; i < arr.n; i++)
      in += (pa[i].id != pa[i - 1].id);
  }

  cnt[0] = in;
  comm_scan(out, ci, gs_long, gs_add, cnt, 1, wrk);
  crs->ng[1] = out[1][0];
  slong s = crs->ng[0] + out[0][0] + 1;

  if (in) {
    struct dof_t *pa = (struct dof_t *)arr.ptr;
    i = 0;
    while (i < arr.n) {
      for (j = i + 1; j < arr.n && pa[j].id == pa[i].id; j++)
        ;
      for (; i < j; i++)
        pa[i].nid = s;
      s++;
    }
  }

  struct crystal cr;
  crystal_init(&cr, ci);
  sarray_transfer(struct dof_t, &arr, p0, 0, &cr);
  crystal_free(&cr);

  sarray_sort(struct dof_t, arr.ptr, arr.n, id, 1, bfr);
  struct dof_t *pa = (struct dof_t *)arr.ptr;
  for (i = 0; i < arr.n; i++)
    nid[pa[i].idx] = pa[i].nid;

  array_free(&arr);
  comm_free(&c);
}

// n  = ncr * nelt
// nz = ncr * ncr * nelt
struct coarse *crs_parrsb_setup(uint n, const ulong *id, uint nz,
                                const uint *Ai, const uint *Aj, const scalar *A,
                                unsigned null_space, unsigned type,
                                const struct comm *c) {
  comm_barrier(c);
  double tcrs = comm_time();

  // crs->un is the user vector size.
  struct coarse *crs = tcalloc(struct coarse, 1);
  crs->null_space = null_space, crs->type = type, crs->un = n;
  for (unsigned i = 0; i < 3; i++)
    crs->ng[i] = crs->s[i] = crs->n[i] = 0;

  // Setup the buffer and duplicate the communicator.
  buffer_init(&crs->bfr, 1024);
  comm_dup(&crs->c, c);

  // Let's renumber the ids just in case its the schur solver. Schwarz solver
  // doesn't need re-numbering but we are going to go ahead and do it.
  slong *tid = tcalloc(slong, crs->un);
  for (uint i = 0; i < n; i++)
    tid[i] = id[i];

  // Find the mapping from user ids to unique ids (compressed ids) local to the
  // processor. Compressed vector size is `crs->cn`.
  ulong *uid = tcalloc(ulong, crs->un);
  crs->u2c = tcalloc(sint, crs->un);
  crs->cn = unique_ids(crs->u2c, uid, crs->un, tid, &crs->bfr);
#if 0
  for (uint i = 0; i < crs->un; i++) {
    printf("p = %d i = %u perm[i] = %d\n", c->id, i, crs->u2c[i]);
    fflush(stdout);
  }
#endif

  // Now renumber unique ids based on whether they are internal or on interface.
  slong *nid = tcalloc(slong, crs->cn);
  number_dofs(nid, crs, tid, uid);
  free(tid), free(uid);

  // Now let's setup the coarse system. Create `struct mij` entries and pass
  // them into schur setup. Which processor owns the dof? All the local dofs
  // are owned by those specific preocessors -- interface dofs are owned in
  // a load balanced manner.
  uint nr = crs->ng[1] / c->np, nrem = crs->ng[1] - nr * c->np;
  uint p0 = c->np - nrem;
  ulong s0 = p0 * nr;

  struct array mijs;
  array_init(struct mij, &mijs, n);

  struct mij m = {.r = 0, .c = 0, .idx = 0, .p = 0, .v = 0};
  for (uint k = 0; k < nz; k++) {
    sint i = crs->u2c[Ai[k]], j = crs->u2c[Aj[k]];
    if (i < 0 || j < 0 || A[k] == 0)
      continue;
    m.r = nid[i], m.c = nid[j], m.v = A[k], m.p = c->id;
    if (m.r > crs->ng[0]) {
      if (m.r - crs->ng[0] <= s0)
        m.p = (m.r - crs->ng[0] - 1) / nr;
      else
        m.p = p0 + (m.r - crs->ng[0] - s0 - 1) / (nr + 1);
    }
    array_cat(struct mij, &mijs, &m, 1);
  }

  // Now let's assemble the matrix by sending load balancing the interface rows.
  // Assembled size is `an`.
  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct mij, &mijs, p, 1, &cr);

  nid = trealloc(slong, nid, crs->cn + crs->n[0] + nr + 1);
  for (uint i = 0; i < crs->cn; i++)
    nid[i] = -nid[i];

  crs->an = 0;
  if (mijs.n > 0) {
    sarray_sort_2(struct mij, mijs.ptr, mijs.n, r, 1, c, 1, &crs->bfr);
    struct mij *pm = (struct mij *)mijs.ptr;
    uint i = 0, j;
    while (i < mijs.n) {
      for (j = i + 1; j < mijs.n && pm[j].r == pm[i].r; j++)
        ;
      nid[crs->cn + crs->an] = pm[i].r, crs->an++, i = j;
    }
  }
  crs->n[1] = crs->an - crs->n[0];
  crs->s[1] = nid[crs->cn + crs->n[0]];
  crs->c2a = gs_setup(nid, crs->cn + crs->an, c, 0, gs_pairwise, 0);

  tcrs = comm_time() - tcrs;
  double wrk, min = tcrs, max = tcrs;
  comm_allreduce(c, gs_double, gs_max, &max, 1, &wrk);
  comm_allreduce(c, gs_double, gs_min, &min, 1, &wrk);
  if (c->id == 0) {
    printf("parrsb_crs_setup: %g %g (min max)\n", min, max);
    fflush(stdout);
  }

  comm_barrier(c);
  tcrs = comm_time();

  switch (type) {
  case 0:
    schur_setup(crs, &mijs, &cr, &crs->bfr);
    break;
  default:
    break;
  }

  min = max = comm_time() - tcrs;
  comm_allreduce(c, gs_double, gs_max, &max, 1, &wrk);
  comm_allreduce(c, gs_double, gs_min, &min, 1, &wrk);
  if (c->id == 0) {
    printf("schur_setup: %g %g (min max)\n", min, max);
    fflush(stdout);
  }

  array_free(&mijs), crystal_free(&cr);

  return crs;
}

void crs_parrsb_solve(scalar *x, struct coarse *crs, scalar *b, scalar tol) {
  metric_init();

  scalar *rhs = tcalloc(scalar, crs->cn + crs->an);
  for (uint i = 0; i < crs->un; i++) {
    if (crs->u2c[i] >= 0)
      rhs[crs->u2c[i]] += b[i];
  }

#if 0
  for (uint i = 0; i < crs->cn; i++) {
    printf("p = %d i = %u before b[i] = %lf\n", crs->c.id, i, rhs[i]);
    fflush(stdout);
  }
#endif

  gs(rhs, gs_double, gs_add, 1, crs->c2a, &crs->bfr);

#if 0
  char name[BUFSIZ];
  snprintf(name, BUFSIZ, "rsb_b_np_%d_id_%d_nl_%lld_ni_%lld.txt", crs->c.np,
           crs->c.id, crs->n[0], crs->n[1]);
  FILE *fp = fopen(name, "w");
  if (fp) {
    for (uint i = 0; i < crs->an; i++)
      fprintf(fp, "%lf\n", rhs[crs->cn + i]);
    fclose(fp);
  }
#endif

#if 0
  for (uint i = 0; i < crs->an; i++) {
    printf("p = %d i = %u after b[i] = %lf\n", crs->c.id, i, rhs[crs->cn + i]);
    fflush(stdout);
  }
#endif

  switch (crs->type) {
  case 0:
    schur_solve(rhs + crs->cn, crs, rhs + crs->cn, tol, &crs->bfr);
    break;
  default:
    break;
  }

#if 0
  for (uint i = 0; i < crs->an; i++) {
    printf("p = %d i = %u x[i] = %lf w[i] = %lf\n", crs->c.id, i,
           rhs[crs->cn + i], weights[crs->cn + i]);
    fflush(stdout);
  }
#endif

  gs(rhs, gs_double, gs_add, 0, crs->c2a, &crs->bfr);
  for (uint i = 0; i < crs->un; i++) {
    if (crs->u2c[i] >= 0)
      x[i] = rhs[crs->u2c[i]];
    else
      x[i] = 0;
  }
  free(rhs);

#if 0
  snprintf(name, BUFSIZ, "rsb_x_np_%d_id_%d_un_%u.txt", crs->c.np, crs->c.id,
           crs->un);
  fp = fopen(name, "w");
  if (fp) {
    for (uint i = 0; i < crs->un; i++)
      fprintf(fp, "%lf\n", x[i]);
    fclose(fp);
  }
#endif

  metric_push_level();
  metric_crs_print(&crs->c, 1);
  metric_finalize();
}

void crs_parrsb_free(struct coarse *crs) {
  if (crs != NULL) {
    switch (crs->type) {
    case 0:
      schur_free(crs);
      break;
    default:
      break;
    }
    if (crs->u2c)
      free(crs->u2c);
    gs_free(crs->c2a);
    comm_free(&crs->c), buffer_free(&crs->bfr);
    free(crs), crs = NULL;
  }
}

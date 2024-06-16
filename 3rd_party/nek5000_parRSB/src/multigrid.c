#include "multigrid.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

struct mg_lvl {
  uint npres, nposts;
  scalar over;
  struct gs_data *J; // Interpolation from level l to l + 1

  struct gs_data *Q; // gs handle for matrix vector product
  struct par_mat *M; // Operator
};

struct mg {
  uint nlevels, *level_off;
  struct mg_lvl **levels;
  scalar *buf;
};

//=============================================================================
// MG setup
//
static scalar sigma_cheb(int k, int n, scalar lmin, scalar lmax) {
  k = (k - 1) % n + 1;
  scalar theta = M_PI * (k - 0.5) / n;
  scalar lamk = lmin + 0.5 * (lmax - lmin) * (cos(theta) + 1);
  return 1 / lamk;
}

inline static void set_proc(struct mij *m, uint nelt, uint nrem, uint np) {
  assert(m->r > 0);

  if (nrem == 0) {
    m->p = (m->r - 1) / nelt;
  } else {
    uint s = np - nrem;
    ulong t = nelt * s;
    if (m->r <= t)
      m->p = (m->r - 1) / nelt;
    else
      m->p = s + (m->r - (t + 1)) / (nelt + 1);
  }

  assert(m->p < np);
}

static int sparse_gemm(struct par_mat *WG, const struct par_mat *W,
                       const struct par_mat *G, int diag_wg, struct crystal *cr,
                       buffer *bfr) {
  // W is in CSR, G is in CSC; we multiply rows of W by shifting
  // the columns of G from processor to processor. This is not scalable
  // at all -- need to do a 2D partition of the matrices W and G.
  assert(IS_CSR(W) && !IS_DIAG(W));
  assert(IS_CSC(G));

  // Put G into an array to transfer from processor to processor
  struct array gij, sij;
  array_init(struct mij, &gij, 100);
  array_init(struct mij, &sij, 100);

  struct mij m = {.r = 0, .c = 0, .idx = 0, .p = cr->comm.id, .v = 0};
  uint i, j, je;
  for (i = 0; i < G->cn; i++) {
    m.c = G->cols[i];
    for (j = G->adj_off[i], je = G->adj_off[i + 1]; j != je; j++) {
      m.r = G->rows[G->adj_idx[j]];
      m.v = G->adj_val[j];
      array_cat(struct mij, &gij, &m, 1);
    }
  }
  if (IS_DIAG(G)) {
    for (i = 0; i < G->cn; i++) {
      m.c = m.r = G->cols[i];
      m.v = G->diag_val[i];
      array_cat(struct mij, &gij, &m, 1);
    }
  }

  sarray_sort_2(struct mij, gij.ptr, gij.n, c, 1, r, 1, bfr);
  struct mij *pg = (struct mij *)gij.ptr;
  for (i = 0; i < gij.n; i++) pg[i].idx = i;

  for (uint p = 0; p < cr->comm.np; p++) {
    // Calculate dot product of each row of W with columns of G
    for (i = 0; i < W->rn; i++) {
      m.r = W->rows[i];
      uint s = 0, e = 0;
      while (s < gij.n) {
        m.c = pg[s].c, m.v = 0;
        for (j = W->adj_off[i], je = W->adj_off[i + 1]; j < je; j++) {
          ulong k = W->cols[W->adj_idx[j]];
          while (e < gij.n && pg[s].c == pg[e].c && pg[e].r < k) e++;
          if (e < gij.n && pg[s].c == pg[e].c && pg[e].r == k)
            m.v += W->adj_val[j] * pg[e].v;
        }
        while (e < gij.n && pg[s].c == pg[e].c) e++;
        if (fabs(m.v) > 1e-12) array_cat(struct mij, &sij, &m, 1);
        s = e;
      }
    }

    sint next = (cr->comm.id + 1) % cr->comm.np;
    for (i = 0; i < gij.n; i++) pg[i].p = next;
    sarray_transfer(struct mij, &gij, p, 0, cr);

    sarray_sort(struct mij, gij.ptr, gij.n, idx, 0, bfr);
    pg = gij.ptr;
  }

  par_csr_setup(WG, &sij, diag_wg, bfr);
  array_free(&gij), array_free(&sij);

  return 0;
}

static uint mg_setup_aux(struct mg *d, const int factor, struct crystal *cr,
                         struct array *mijs, buffer *bfr) {
  uint lvl = d->nlevels;
  struct mg_lvl *l = d->levels[lvl - 1];

  struct par_mat *M = l->M;
  uint nnz = ((M->rn > 0) ? (M->adj_off[M->rn] + M->rn) : 0);

  struct mij m = {.r = 0, .c = 0, .idx = 0, .p = 0, .v = 0};
  array_reserve(struct mij, mijs, nnz);

  // Now we interpolate to find the coarse operator Mc = J^T M J
  // Calculate coarse level parameters: ngc, npc, nelt, nrem
  uint size = (M->rn > 0 ? (M->rows[M->rn - 1] - M->rows[0] + 1) : 0);
  slong ng = size, wrk[2][1];
  struct comm *c = &cr->comm;
  comm_allreduce(c, gs_long, gs_add, &ng, 1, wrk);

  // ng > 1 based on while condition in mg_setup(). so ngc >= 1
  ulong ngc = (ng + factor - 1) / factor;
  const uint npc = (ngc < c->np ? ngc : c->np);
  const uint nelt = ngc / npc, nrem = ngc - nelt * npc;

  // Calculate the minimum row id
  slong rs = (M->rn > 0 ? M->rows[0] : LLONG_MAX);
  comm_allreduce(c, gs_long, gs_min, &rs, 1, wrk);
  assert(rs > 0);
  rs -= 1;

  // Reserve enough memory for ids used for interpolation
  slong *ids = (slong *)tcalloc(slong, 2 * M->rn);

  mijs->n = 0;
  uint i, k;
  for (i = 0, k = 0; i < M->rn; i++) {
    m.r = (M->rows[i] - rs + factor - 1) / factor;
    set_proc(&m, nelt, nrem, npc);
    for (uint j = M->adj_off[i], je = M->adj_off[i + 1]; j < je; j++) {
      m.c = (M->cols[M->adj_idx[j]] - rs + factor - 1) / factor;
      m.v = M->adj_val[j];
      array_cat(struct mij, mijs, &m, 1);
    }
    m.c = m.r, m.v = M->diag_val[i];
    ids[k++] = -m.c;
    array_cat(struct mij, mijs, &m, 1);
  }

  sarray_transfer(struct mij, mijs, p, 0, cr);
  sarray_sort_2(struct mij, mijs->ptr, mijs->n, r, 1, c, 1, bfr);

  d->nlevels++;
  d->levels = trealloc(struct mg_lvl *, d->levels, d->nlevels);
  l = d->levels[lvl] = tcalloc(struct mg_lvl, 1);
  M = l->M = par_csr_setup_ext(mijs, 1, bfr);

  // Setup gs ids for coarse level (rhs interpolation )
  ids = (slong *)trealloc(slong, ids, k + M->rn);
  for (i = 0; i < M->rn; i++) ids[k++] = M->rows[i];
  d->levels[lvl - 1]->J = gs_setup(ids, k, c, 0, gs_pairwise, 0);
  free(ids);

  d->levels[lvl]->Q = setup_Q(M, c, bfr);
  return lvl;
}

struct mg *mg_setup(const struct par_mat *M, const int factor,
                    struct crystal *cr, buffer *bfr) {
  assert(IS_CSR(M));
  assert(M->rn == 0 || IS_DIAG(M));

  // Allocate memory for struct mg
  struct mg *d = (struct mg *)tcalloc(struct mg, 1);

  // Setup Level 1, keeps a pointer to input matrix
  d->nlevels = 1;
  d->levels = (struct mg_lvl **)tcalloc(struct mg_lvl *, d->nlevels);
  d->level_off = tcalloc(uint, d->nlevels + 1);

  d->levels[0] = (struct mg_lvl *)tcalloc(struct mg_lvl, 1);
  d->levels[0]->npres = 3;
  d->levels[0]->nposts = 0;
  d->levels[0]->over = 1.33333;
  d->levels[0]->M = (struct par_mat *)M;

  struct comm *c = &cr->comm;
  d->levels[0]->Q = setup_Q(M, c, bfr);

  uint size = (M->rn > 0 ? (M->rows[M->rn - 1] - M->rows[0] + 1) : 0);
  d->level_off[0] = 0, d->level_off[1] = size;

  uint nnz = (M->rn > 0 ? M->adj_off[M->rn] + M->rn : 0);
  struct array mijs;
  array_init(struct mij, &mijs, nnz);

  slong wrk[2], ng = size;
  comm_allreduce(c, gs_long, gs_add, &ng, 1, wrk);
  while (ng > 1) {
    uint l = mg_setup_aux(d, factor, cr, &mijs, bfr);
    struct par_mat *Ml = d->levels[l]->M;
    if (Ml->rn > 0 && Ml->adj_off[Ml->rn] + Ml->rn > nnz)
      nnz = Ml->adj_off[Ml->rn] + Ml->rn;
    d->levels[l]->npres = 3;
    d->levels[l]->nposts = 0;
    d->levels[l]->over = 1.5;

    d->level_off = trealloc(uint, d->level_off, d->nlevels + 1);
    size = (Ml->rn > 0 ? (Ml->rows[Ml->rn - 1] - Ml->rows[0] + 1) : 0);
    d->level_off[l + 1] = d->level_off[l] + size;

    ng = size;
    comm_allreduce(c, gs_long, gs_add, &ng, 1, wrk);
  }

  d->levels[d->nlevels - 1]->J = NULL;
  d->buf = tcalloc(scalar, 5 * d->level_off[d->nlevels] + nnz);

  array_free(&mijs);

  return d;
}

//==============================================================================
// MG V-cycle and related functions
//
void mg_vcycle(scalar *u1, scalar *rhs, struct mg *d, struct comm *c,
               buffer *bfr) {
  if (d->nlevels == 0) return;

  uint *lvl_off = d->level_off, nnz = lvl_off[d->nlevels];
  scalar *r = d->buf;
  for (uint i = 0; i < 4 * nnz; i++) r[i] = 0;
  for (uint i = 0; i < lvl_off[1]; i++) r[i] = rhs[i];

  scalar *s = r + nnz, *Gs = s + nnz, *u = Gs + nnz, *wrk = u + nnz;

  uint i, j, n, off;
  for (uint lvl = 0; lvl < d->nlevels - 1; lvl++) {
    off = lvl_off[lvl];
    n = lvl_off[lvl + 1] - off;

    struct mg_lvl *l = d->levels[lvl];
    struct par_mat *M = l->M;

    // u = sigma * inv(D) * rhs
    scalar sigma = sigma_cheb(1, l->npres + 1, 1, 2);
    for (j = 0; j < n; j++) u[off + j] = sigma * r[off + j] / M->diag_val[j];

    // G*u
    mat_vec_csr(Gs + off, u + off, M, l->Q, wrk, bfr);

    // r = rhs - Gu
    for (j = 0; j < n; j++) r[off + j] = r[off + j] - Gs[off + j];

    for (i = 1; i <= l->npres - 1; i++) {
      sigma = sigma_cheb(i + 1, l->npres + 1, 1, 2);

      // s = sigma * inv(D) * r
      // u = u + s
      for (j = 0; j < n; j++) {
        s[off + j] = sigma * r[off + j] / M->diag_val[j];
        u[off + j] += s[off + j];
      }

      // r = r - Gs
      mat_vec_csr(Gs + off, s + off, M, l->Q, wrk, bfr);
      for (j = 0; j < n; j++) r[off + j] = r[off + j] - Gs[off + j];
    }

    // Interpolate to coarser level
    gs(r + off, gs_double, gs_add, 1, l->J, bfr);
  }

  // Coarsest level
  off = lvl_off[d->nlevels - 1];
  n = lvl_off[d->nlevels] - off;
  if (n == 1) {
    struct mg_lvl *l = d->levels[d->nlevels - 1];
    struct par_mat *M = l->M;
    if (fabs(M->diag_val[0]) > 1e-6)
      u[off] = r[off] / M->diag_val[0];
    else
      u[off] = 0.0;
    r[off] = u[off];
  }

  for (int lvl = (int)d->nlevels - 2; lvl >= 0; lvl--) {
    struct mg_lvl *l = d->levels[lvl];
    off = lvl_off[lvl];
    // J*e
    gs(r + off, gs_double, gs_add, 0, l->J, bfr);

    // u = u + over*S*J*e
    n = lvl_off[lvl + 1] - off;
    for (j = 0; j < n; j++) r[off + j] = l->over * r[off + j] + u[off + j];
  }

  // Avoid this
  for (i = 0; i < lvl_off[1]; i++) u1[i] = r[i];
}

void mg_free(struct mg *d) {
  if (d != NULL) {
    struct mg_lvl **l = d->levels;
    for (uint i = 0; i < d->nlevels; i++) {
      if (i > 0 && l[i]->M != NULL) par_mat_free(l[i]->M), free(l[i]->M);
      if (l[i]->J != NULL) gs_free(l[i]->J), l[i]->J = NULL;
      if (l[i]->Q != NULL) gs_free(l[i]->Q), l[i]->Q = NULL;
      if (l[i] != NULL) free(l[i]), l[i] = NULL;
    }

    if (d->levels != NULL) free(d->levels), d->levels = NULL;
    if (d->level_off != NULL) free(d->level_off), d->level_off = NULL;
    if (d->buf != NULL) free(d->buf), d->buf = NULL;
    free(d);
  }
}

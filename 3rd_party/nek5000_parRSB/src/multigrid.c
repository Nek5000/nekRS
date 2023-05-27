#include "multigrid.h"
#include <math.h>

struct mg_lvl {
  uint npres, nposts;
  scalar over;
  struct gs_data *J; // Interpolation from level l to l + 1

  struct gs_data *Q; // gs handle for matrix vector product
  struct par_mat *M; // Operator

  struct gs_data *Qs, *Qst; // gs handle for matrix vector product
  struct par_mat *S, *St;   // Smooth aggregation
};

struct mg {
  uint sagg, nlevels, *level_off;
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

static void inline set_proc(struct mij *m, uint nelt, uint nrem, uint np) {
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

  assert(m->p >= 0 && m->p < np);
}

extern int sparse_gemm(struct par_mat *WG, const struct par_mat *W,
                       const struct par_mat *G, int diag_wg, struct crystal *cr,
                       buffer *bfr);

static uint mg_setup_aux(struct mg *d, const int factor, const int sagg,
                         struct crystal *cr, struct array *mijs, buffer *bfr) {
  uint lvl = d->nlevels;
  struct mg_lvl *l = d->levels[lvl - 1];

  struct par_mat *Ml = l->M;
  uint nnz = ((Ml->rn > 0) ? (Ml->adj_off[Ml->rn] + Ml->rn) : 0);

  struct mij m = {.r = 0, .c = 0, .idx = 0, .p = 0, .v = 0};
  array_reserve(struct mij, mijs, nnz);

  struct comm *c = &cr->comm;
  const double sigma = 0.65;
  struct par_mat *M;
  // Replace M by the following if smooth aggregation is used:
  // S = (I - sigma * D^{-1} * Ml)
  // M = ST * Ml * S
  if (sagg) {
    // This is very hacky and not optimal at all. Should be rewritten.
    // Create S is in CSR format, with separate diagonal. Then convert
    // to CSC with no separate diagonal in order to do the mat-vec.
    mijs->n = 0;
    for (uint i = 0; i < Ml->rn; i++) {
      m.c = m.r = Ml->rows[i], m.v = 1 - sigma;
      array_cat(struct mij, mijs, &m, 1);
      double di = 1.0 / Ml->diag_val[i];
      for (uint j = Ml->adj_off[i], je = Ml->adj_off[i + 1]; j < je; j++) {
        m.c = Ml->cols[Ml->adj_idx[j]];
        m.v = -sigma * di * Ml->adj_val[j];
        array_cat(struct mij, mijs, &m, 1);
      }
    }
    l->S = tcalloc(struct par_mat, 1);
    par_mat_setup(l->S, mijs, 1, 1, bfr);
    l->Qs = setup_Q(l->S, c, bfr);

    struct par_mat S;
    par_csr_to_csc(&S, l->S, 0, cr, bfr);

    // Create N = M in CSR format, no separate diagonal.
    mijs->n = 0;
    for (uint i = 0; i < Ml->rn; i++) {
      m.c = m.r = Ml->rows[i], m.v = Ml->diag_val[i];
      array_cat(struct mij, mijs, &m, 1);
      for (uint j = Ml->adj_off[i], je = Ml->adj_off[i + 1]; j < je; j++) {
        m.c = Ml->cols[Ml->adj_idx[j]];
        m.v = Ml->adj_val[j];
        array_cat(struct mij, mijs, &m, 1);
      }
    }
    struct par_mat N;
    par_mat_setup(&N, mijs, 1, 0, bfr);

    // T = N * S, CSR format, no separate diagonal.
    struct par_mat T;
    sparse_gemm(&T, &N, &S, 0, cr, bfr);
    par_mat_free(&N), par_mat_free(&S);

    // N = T, CSC format, no separate diagonal.
    par_csr_to_csc(&N, &T, 0, cr, bfr);
    par_mat_free(&T);

    // Setup S^t, CSR format, no separate diagonal.
    mijs->n = 0;
    for (uint i = 0; i < Ml->rn; i++) {
      m.c = m.r = Ml->rows[i], m.v = 1 - sigma;
      array_cat(struct mij, mijs, &m, 1);
      double di = 1.0 / Ml->diag_val[i];
      for (uint j = Ml->adj_off[i], je = Ml->adj_off[i + 1]; j < je; j++) {
        m.r = Ml->cols[Ml->adj_idx[j]];
        m.v = -sigma * di * Ml->adj_val[j];
        array_cat(struct mij, mijs, &m, 1);
      }
    }
    par_mat_setup(&T, mijs, 0, 0, bfr);
    par_csc_to_csr(&S, &T, 0, cr, bfr);
    par_mat_free(&T);

    // M = ST * N
    M = tcalloc(struct par_mat, 1);
    sparse_gemm(M, &S, &N, 1, cr, bfr);
    par_mat_free(&S), par_mat_free(&N);

    // Normalize M by the largest value
    double max = 0;
    for (uint i = 0; i < M->rn; i++) {
      for (uint j = M->adj_off[i], je = M->adj_off[i + 1]; j < je; j++)
        if (fabs(M->adj_val[j]) > max)
          max = fabs(M->adj_val[j]);
      if (fabs(M->diag_val[i]) > max)
        max = fabs(M->diag_val[i]);
    }
    double wrk[2];
    comm_allreduce(c, gs_double, gs_max, &max, 1, wrk);

    for (uint i = 0; i < M->rn; i++) {
      for (uint j = M->adj_off[i], je = M->adj_off[i + 1]; j < je; j++)
        M->adj_val[j] /= max;
      M->diag_val[i] /= max;
    }

    par_mat_setup(&T, mijs, 0, 0, bfr);
    l->St = tcalloc(struct par_mat, 1);
    par_csc_to_csr(l->St, &T, 1, cr, bfr);
    par_mat_free(&T);
    l->Qst = setup_Q(l->St, c, bfr);
  } else {
    l->S = l->St = NULL;
    l->Qs = l->Qst = NULL;
    M = Ml;
  }

  // Now we interpolate to find the coarse operator Mc = J^T M J
  // Calculate coarse level parameters: ngc, npc, nelt, nrem
  uint size = (M->rn > 0 ? (M->rows[M->rn - 1] - M->rows[0] + 1) : 0);
  slong ng = size, wrk[2][1];
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

  if (sagg) {
    par_mat_free(M);
    free(M);
  }

  sarray_transfer(struct mij, mijs, p, 0, cr);
  sarray_sort_2(struct mij, mijs->ptr, mijs->n, r, 1, c, 1, bfr);

  d->nlevels++;
  d->levels = trealloc(struct mg_lvl *, d->levels, d->nlevels);
  l = d->levels[lvl] = tcalloc(struct mg_lvl, 1);
  M = l->M = par_csr_setup_ext(mijs, 1, bfr);

  // Setup gs ids for coarse level (rhs interpolation )
  ids = (slong *)trealloc(slong, ids, k + M->rn);
  for (i = 0; i < M->rn; i++)
    ids[k++] = M->rows[i];
  d->levels[lvl - 1]->J = gs_setup(ids, k, c, 0, gs_pairwise, 0);
  free(ids);

  d->levels[lvl]->Q = setup_Q(M, c, bfr);
  return lvl;
}

struct mg *mg_setup(const struct par_mat *M, const int factor, const int sagg,
                    struct crystal *cr, buffer *bfr) {
  assert(IS_CSR(M));
  assert(M->rn == 0 || IS_DIAG(M));

  // Allocate memory for struct mg
  struct mg *d = (struct mg *)tcalloc(struct mg, 1);
  d->sagg = sagg;

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
    uint l = mg_setup_aux(d, factor, sagg, cr, &mijs, bfr);
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
  if (d->nlevels == 0)
    return;

  uint *lvl_off = d->level_off, nnz = lvl_off[d->nlevels];
  scalar *r = d->buf;
  for (uint i = 0; i < 4 * nnz; i++)
    r[i] = 0;
  for (uint i = 0; i < lvl_off[1]; i++)
    r[i] = rhs[i];

  scalar *s = r + nnz, *Gs = s + nnz, *u = Gs + nnz, *wrk = u + nnz;

  uint i, j, n, off;
  for (int lvl = 0; lvl < d->nlevels - 1; lvl++) {
    off = lvl_off[lvl];
    n = lvl_off[lvl + 1] - off;

    struct mg_lvl *l = d->levels[lvl];
    struct par_mat *M = l->M;

    // u = sigma * inv(D) * rhs
    scalar sigma = sigma_cheb(1, l->npres + 1, 1, 2);
    for (j = 0; j < n; j++)
      u[off + j] = sigma * r[off + j] / M->diag_val[j];

    // G*u
    mat_vec_csr(Gs + off, u + off, M, l->Q, wrk, bfr);

    // r = rhs - Gu
    for (j = 0; j < n; j++)
      r[off + j] = r[off + j] - Gs[off + j];

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
      for (j = 0; j < n; j++)
        r[off + j] = r[off + j] - Gs[off + j];
    }

    // Apply S^T
    if (d->sagg)
      mat_vec_csr(r + off, r + off, l->St, l->Qst, wrk, bfr);

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

  for (int lvl = d->nlevels - 2; lvl >= 0; lvl--) {
    struct mg_lvl *l = d->levels[lvl];
    off = lvl_off[lvl];
    // J*e
    gs(r + off, gs_double, gs_add, 0, l->J, bfr);

    // Apply S
    if (d->sagg)
      mat_vec_csr(r + off, r + off, l->S, l->Qs, wrk, bfr);

    // u = u + over*S*J*e
    n = lvl_off[lvl + 1] - off;
    for (j = 0; j < n; j++)
      r[off + j] = l->over * r[off + j] + u[off + j];
  }

  // Avoid this
  for (i = 0; i < lvl_off[1]; i++)
    u1[i] = r[i];
}

void mg_free(struct mg *d) {
  if (d != NULL) {
    struct mg_lvl **l = d->levels;
    for (uint i = 0; i < d->nlevels; i++) {
      if (i > 0 && l[i]->M != NULL)
        par_mat_free(l[i]->M), free(l[i]->M);
      if (l[i]->J != NULL)
        gs_free(l[i]->J), l[i]->J = NULL;
      if (l[i]->Q != NULL)
        gs_free(l[i]->Q), l[i]->Q = NULL;
      if (l[i]->Qs != NULL)
        gs_free(l[i]->Qs), l[i]->Qs = NULL;
      if (l[i]->Qst != NULL)
        gs_free(l[i]->Qst), l[i]->Qst = NULL;
      if (l[i]->S != NULL)
        par_mat_free(l[i]->S), l[i]->S = NULL;
      if (l[i]->St != NULL)
        par_mat_free(l[i]->St), l[i]->St = NULL;
      if (l[i] != NULL)
        free(l[i]), l[i] = NULL;
    }

    if (d->levels != NULL)
      free(d->levels), d->levels = NULL;
    if (d->level_off != NULL)
      free(d->level_off), d->level_off = NULL;
    if (d->buf != NULL)
      free(d->buf), d->buf = NULL;
    free(d);
  }
}

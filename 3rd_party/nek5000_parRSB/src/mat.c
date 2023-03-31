#include "mat.h"

#define CSC 0
#define CSR 1

#define FREE(ptr, x)                                                           \
  {                                                                            \
    if (ptr->x != NULL)                                                        \
      free(ptr->x);                                                            \
  }

//------------------------------------------------------------------------------
// `mat` struct for local matrices
//
// Compress array by summing up entries which share the same (r, c) values
// and the array is modified in place. Also, the diagonal entries are modified
// to ensure all the
int compress_nbrs(struct array *eij, struct array *nbr, buffer *bfr) {
  array_init(struct mij, eij, nbr->n);
  if (nbr->n == 0)
    return 1;

  sarray_sort_2(struct nbr, nbr->ptr, nbr->n, r, 1, c, 1, bfr);
  struct nbr *ptr = (struct nbr *)nbr->ptr;

  struct mij m;
  m.idx = 0;

  sint i = 0;
  while (i < nbr->n) {
    m.r = ptr[i].r, m.c = ptr[i].c;

    sint j = i + 1;
    while (j < nbr->n && ptr[j].r == ptr[i].r && ptr[j].c == ptr[i].c)
      j++;

    m.v = i - j; // = - (j - i)
    array_cat(struct mij, eij, &m, 1);
    i = j;
  }

  // Now make sure the row sum is zero
  struct mij *pe = (struct mij *)eij->ptr;
  i = 0;
  while (i < eij->n) {
    sint j = i, k = -1, s = 0;
    while (j < eij->n && pe[j].r == pe[i].r) {
      if (pe[j].r == pe[j].c)
        k = j;
      else
        s += pe[j].v;
      j++;
    }
    assert(k >= 0);
    pe[k].v = -s;
    i = j;
  }
}

int csr_setup(struct mat *mat, struct array *entries, int sep, buffer *buf) {
  uint nnz = entries->n;
  if (nnz == 0) {
    mat->start = mat->n = 0;
    mat->Lp = mat->Li = NULL;
    mat->L = mat->D = NULL;
    return 0;
  }

  // Renumber cols and rows
  sarray_sort(struct mij, entries->ptr, entries->n, c, 1, buf);
  struct mij *ptr = (struct mij *)entries->ptr;

  ulong n = ptr[0].c;
  ptr[0].c = 0;
  for (uint i = 1; i < nnz; i++) {
    if (ptr[i].c == n)
      ptr[i].c = ptr[i - 1].c;
    else
      n = ptr[i].c, ptr[i].c = ptr[i - 1].c + 1;
  }

  sarray_sort(struct mij, entries->ptr, entries->n, r, 1, buf);
  ptr = (struct mij *)entries->ptr;
  n = mat->start = ptr[0].r, ptr[0].r = 0;

  for (uint i = 1; i < nnz; i++)
    if (ptr[i].r == n)
      ptr[i].r = ptr[i - 1].r;
    else
      n = ptr[i].r, ptr[i].r = ptr[i - 1].r + 1;
  uint nr = ptr[nnz - 1].r + 1;

  // Reserve enough memory for work arrays
  buffer_reserve(buf, sizeof(struct mij) * nnz);
  struct mij *unique = (struct mij *)buf->ptr;

  // Setup the matrix, separate diagonal
  mat->n = nr;
  uint *Lp = mat->Lp = (uint *)tcalloc(uint, nr + 1);

  sarray_sort_2(struct mij, entries->ptr, entries->n, r, 1, c, 1, buf);
  ptr = (struct mij *)entries->ptr;

  sep = sep != 0;
  Lp[0] = 0, unique[0] = ptr[0];
  uint i, j;
  for (nr = 1, i = 1, j = 0; i < nnz; i++) {
    if ((unique[j].r != ptr[i].r) || (unique[j].c != ptr[i].c)) {
      if (unique[j].r != ptr[i].r)
        Lp[nr] = j + 1 - sep * nr, nr++;
      unique[++j] = ptr[i];
    } else
      unique[j].v += ptr[i].v;
  }
  Lp[mat->n] = ++j - sep * mat->n;

  mat->Li = (uint *)tcalloc(uint, j - sep * mat->n);
  mat->L = (scalar *)tcalloc(scalar, j - sep * mat->n);

  uint nadj;
  if (sep) {
    mat->D = (scalar *)tcalloc(scalar, mat->n);
    uint ndiag;
    for (i = ndiag = nadj = 0; i < j; i++) {
      if (mat->start + unique[i].r == mat->start + unique[i].c)
        mat->D[ndiag++] = unique[i].v;
      else
        mat->L[nadj] = unique[i].v, mat->Li[nadj++] = unique[i].c;
    }
  } else {
    mat->D = NULL;
    for (i = nadj = 0; i < j; i++)
      mat->L[nadj] = unique[i].v, mat->Li[nadj++] = unique[i].c;
  }

  return 0;
}

int mat_print(struct mat *mat) {
  uint i, j;
  for (i = 0; i < mat->n; i++) {
    for (j = mat->Lp[i]; j < mat->Lp[i + 1]; j++)
      printf("%llu %llu %lf\n", mat->start + i, mat->start + mat->Li[j],
             mat->L[j]);
    if (mat->D != NULL)
      printf("%lld %lld %lf\n", mat->start + i, mat->start + i, mat->D[i]);
    fflush(stdout);
  }

  return 0;
}

void mat_dump(const char *name, struct mat *A, struct crystal *cr,
              buffer *bfr) {
  struct array mijs;
  array_init(struct mij, &mijs, 1024);

  struct mij m = {.r = 0, .c = 0, .idx = 0, .p = 0, .v = 0};

  uint i, j;
  for (i = 0; i < A->n; i++) {
    m.r = A->start + i;
    for (j = A->Lp[i]; j < A->Lp[i + 1]; j++) {
      m.c = A->start + A->Li[j], m.v = A->L[j];
      array_cat(struct mij, &mijs, &m, 1);
    }
    if (A->D != NULL) {
      m.c = m.r, m.v = A->D[i];
      array_cat(struct mij, &mijs, &m, 1);
    }
  }

  sarray_transfer(struct mij, &mijs, p, 1, cr);
  sarray_sort_2(struct mij, mijs.ptr, mijs.n, r, 1, c, 1, bfr);

  struct comm *c = &cr->comm;
  if (c->id == 0) {
    FILE *fp = fopen(name, "w+");
    if (fp != NULL) {
      struct mij *pm = (struct mij *)mijs.ptr;
      for (uint i = 0; i < mijs.n; i++)
        fprintf(fp, "%llu %llu %.15lf\n", pm[i].r, pm[i].c, pm[i].v);
      fclose(fp);
    }
  }

  array_free(&mijs);
}

int mat_free(struct mat *mat) {
  FREE(mat, Lp);
  FREE(mat, Li);
  FREE(mat, L);
  FREE(mat, D);
  return 0;
}

//------------------------------------------------------------------------------
// Find neighbors in the graph
//
void find_nbrs(struct array *arr, const ulong *eid, const slong *vtx,
               const uint nelt, const int nv, struct crystal *cr, buffer *buf) {
  struct array vertices;
  array_init(struct nbr, &vertices, nelt * nv);

  struct comm *c = &cr->comm;
  struct nbr v = {.r = 0, .c = 0, .proc = 0};
  uint i, j;
  for (i = 0; i < nelt; i++) {
    v.r = eid[i];
    assert(v.r > 0);
    for (j = 0; j < nv; j++) {
      v.c = vtx[i * nv + j], v.proc = v.c % c->np;
      array_cat(struct nbr, &vertices, &v, 1);
    }
  }

  sarray_transfer(struct nbr, &vertices, proc, 1, cr);
  sarray_sort(struct nbr, vertices.ptr, vertices.n, c, 1, buf);

  // FIXME: Assumes quads or hexes
  struct nbr *pv = (struct nbr *)vertices.ptr, t = {.r = 0, .c = 0, .proc = 0};
  array_init(struct nbr, arr, vertices.n * 10 + 1);
  uint s = 0, e;
  while (s < vertices.n) {
    e = s + 1;
    while (e < vertices.n && pv[s].c == pv[e].c)
      e++;
    for (i = s; i < e; i++) {
      t = pv[i];
      for (j = s; j < e; j++) {
        t.c = pv[j].r;
        assert(t.r > 0 && t.c > 0);
        array_cat(struct nbr, arr, &t, 1);
      }
    }
    s = e;
  }

  sarray_transfer(struct nbr, arr, proc, 1, cr);
  array_free(&vertices);
}

//------------------------------------------------------------------------------
// `par_mat` matrix for parallel distributed matrices
//
int par_csr_setup(struct par_mat *mat, struct array *entries, int sd,
                  buffer *buf) {
  mat->type = CSR;
  if (entries == NULL || entries->n == 0) {
    mat->cn = mat->rn = 0;
    mat->adj_off = mat->adj_idx = mat->diag_idx = NULL;
    mat->adj_val = mat->diag_val = NULL;
    mat->rows = mat->cols = NULL;
    return 0;
  }

  sarray_sort(struct mij, entries->ptr, entries->n, c, 1, buf);
  struct mij *ptr = (struct mij *)entries->ptr;

  // Reserve enough memory for work arrays
  uint nnz = entries->n;
  buffer_reserve(buf, (sizeof(struct mij) + 2 * sizeof(ulong)) * (nnz + 1));

  ulong *cols = (ulong *)buf->ptr;
  cols[0] = ptr[0].c, ptr[0].idx = 0, mat->cn = 1;
  assert(cols[0] > 0);
  uint i;
  for (i = 1; i < nnz; i++) {
    if (ptr[i - 1].c != ptr[i].c)
      cols[mat->cn] = ptr[i].c, ptr[i].idx = mat->cn++;
    else
      ptr[i].idx = ptr[i - 1].idx;
  }

  mat->cols = (ulong *)tcalloc(ulong, mat->cn);
  memcpy(mat->cols, cols, sizeof(ulong) * mat->cn);

  sarray_sort_2(struct mij, entries->ptr, entries->n, r, 1, c, 1, buf);
  ptr = entries->ptr;

  sd = (sd != 0); // sd needs to be 1 or 0

  uint *adj_off = (uint *)buf->ptr;
  ulong *rows = (ulong *)(adj_off + nnz + 1);
  struct mij *unique = (struct mij *)(rows + (nnz + 1));

  adj_off[0] = 0, unique[0] = ptr[0], rows[0] = ptr[0].r;
  uint j = 0;
  for (i = mat->rn = 1; i < nnz; i++) {
    if ((unique[j].r != ptr[i].r) || (unique[j].c != ptr[i].c)) {
      if (unique[j].r != ptr[i].r) {
        adj_off[mat->rn] = j + 1 - sd * mat->rn;
        rows[mat->rn++] = ptr[i].r;
      }
      unique[++j] = ptr[i];
    } else
      unique[j].v += ptr[i].v;
  }
  adj_off[mat->rn] = ++j - sd * mat->rn;

  mat->rows = (ulong *)tcalloc(ulong, mat->rn);
  memcpy(mat->rows, rows, mat->rn * sizeof(ulong));

  mat->adj_off = (uint *)tcalloc(uint, mat->rn + 1);
  memcpy(mat->adj_off, adj_off, (mat->rn + 1) * sizeof(uint));
  mat->adj_idx = (uint *)tcalloc(uint, j - sd * mat->rn);
  mat->adj_val = (scalar *)tcalloc(scalar, j - sd * mat->rn);

  uint nadj = 0, ndiag = 0;
  i = 0;
  if (sd) {
    mat->diag_val = (scalar *)tcalloc(scalar, mat->rn);
    mat->diag_idx = (uint *)tcalloc(uint, mat->rn);
    for (; i < j; i++) {
      if (unique[i].r == unique[i].c) {
        mat->diag_idx[ndiag] = unique[i].idx;
        mat->diag_val[ndiag++] = unique[i].v;
      } else {
        mat->adj_idx[nadj] = unique[i].idx;
        mat->adj_val[nadj++] = unique[i].v;
      }
    }
  } else {
    mat->diag_val = NULL, mat->diag_idx = NULL;
    for (; i < j; i++)
      mat->adj_idx[nadj] = unique[i].idx, mat->adj_val[nadj++] = unique[i].v;
  }

  return 0;
}

int par_csc_setup(struct par_mat *mat, struct array *entries, int sd,
                  buffer *buf) {
  mat->type = CSC;
  if (entries == NULL || entries->n == 0) {
    mat->cn = mat->rn = 0;
    mat->adj_off = mat->adj_idx = mat->diag_idx = NULL;
    mat->adj_val = mat->diag_val = NULL;
    mat->rows = mat->cols = NULL;
    return 0;
  }

  sarray_sort(struct mij, entries->ptr, entries->n, r, 1, buf);
  struct mij *ptr = (struct mij *)entries->ptr;

  // Reserve enough memory for work arrays
  uint nnz = entries->n;
  buffer_reserve(buf, (sizeof(struct mij) + 2 * sizeof(ulong)) * (nnz + 1));

  ulong *rows = (ulong *)buf->ptr;
  rows[0] = ptr[0].r, ptr[0].idx = 0, mat->rn = 1;
  uint i;
  for (i = 1; i < nnz; i++) {
    if (ptr[i - 1].r != ptr[i].r)
      rows[mat->rn] = ptr[i].r, ptr[i].idx = mat->rn++;
    else
      ptr[i].idx = ptr[i - 1].idx;
  }

  mat->rows = tcalloc(ulong, mat->rn);
  memcpy(mat->rows, rows, sizeof(ulong) * mat->rn);

  sarray_sort_2(struct mij, entries->ptr, entries->n, c, 1, r, 1, buf);
  ptr = entries->ptr;

  sd = sd != 0; // sd needs to be 1 or 0

  uint *adj_off = (uint *)buf->ptr;
  ulong *cols = (ulong *)(adj_off + nnz + 1);
  struct mij *unique = (struct mij *)(cols + (nnz + 1));

  adj_off[0] = 0, unique[0] = ptr[0], cols[0] = ptr[0].c;
  uint j = 0;
  for (i = mat->cn = 1; i < nnz; i++) {
    if ((unique[j].r != ptr[i].r) || (unique[j].c != ptr[i].c)) {
      if (unique[j].c != ptr[i].c) {
        adj_off[mat->cn] = j + 1 - sd * mat->cn;
        cols[mat->cn++] = ptr[i].c;
      }
      unique[++j] = ptr[i];
    } else
      unique[j].v += ptr[i].v;
  }
  adj_off[mat->cn] = ++j - sd * mat->cn;

  mat->cols = (ulong *)tcalloc(ulong, mat->cn);
  memcpy(mat->cols, cols, mat->cn * sizeof(ulong));

  mat->adj_off = (uint *)tcalloc(uint, mat->cn + 1);
  memcpy(mat->adj_off, adj_off, (mat->cn + 1) * sizeof(uint));
  mat->adj_idx = (uint *)tcalloc(uint, j - sd * mat->cn);
  mat->adj_val = (scalar *)tcalloc(scalar, j - sd * mat->cn);

  uint nadj = 0, ndiag = 0;
  i = 0;
  if (sd) {
    mat->diag_idx = (uint *)tcalloc(uint, mat->cn);
    mat->diag_val = (scalar *)tcalloc(scalar, mat->cn);
    for (; i < j; i++) {
      if (unique[i].r == unique[i].c) {
        mat->diag_idx[ndiag] = unique[i].idx;
        mat->diag_val[ndiag++] = unique[i].v;
      } else {
        mat->adj_idx[nadj] = unique[i].idx;
        mat->adj_val[nadj++] = unique[i].v;
      }
    }
  } else {
    mat->diag_idx = NULL, mat->diag_val = NULL;
    for (; i < j; i++)
      mat->adj_idx[nadj] = unique[i].idx, mat->adj_val[nadj++] = unique[i].v;
  }

  return 0;
}

int par_mat_setup(struct par_mat *M, struct array *mijs, const int type,
                  const int sd, buffer *bfr) {
  assert(type == CSC || type == CSR);

  M->type = type;
  if (mijs == NULL || mijs->n == 0) {
    M->cn = M->rn = 0;
    M->adj_off = M->adj_idx = M->diag_idx = NULL;
    M->adj_val = M->diag_val = NULL;
    M->rows = M->cols = NULL;
    return 0;
  }

#define INIT_IDX_FLD(f, fn, fld)                                               \
  do {                                                                         \
    sarray_sort(struct mij, mijs->ptr, mijs->n, f, 1, bfr);                    \
    struct mij *ptr = (struct mij *)mijs->ptr;                                 \
    tmp[0] = ptr[0].f, ptr[0].idx = 0, M->fn = 1;                              \
    for (uint i = 1; i < nnz; i++) {                                           \
      if (ptr[i - 1].f != ptr[i].f)                                            \
        tmp[M->fn] = ptr[i].f, ptr[i].idx = M->fn++;                           \
      else                                                                     \
        ptr[i].idx = ptr[i - 1].idx;                                           \
    }                                                                          \
    M->fld = tcalloc(ulong, M->fn);                                            \
    memcpy(M->fld, tmp, sizeof(ulong) * M->fn);                                \
  } while (0)

#define INIT_ADJ_OFF(f, fn, fld)                                               \
  do {                                                                         \
    tmp[0] = ptr[0].f;                                                         \
    uint i, j;                                                                 \
    for (i = M->fn = 1, j = 0; i < nnz; i++) {                                 \
      if ((unique[j].r != ptr[i].r) || (unique[j].c != ptr[i].c)) {            \
        if (unique[j].f != ptr[i].f) {                                         \
          off[M->fn] = j + 1 - sd * M->fn;                                     \
          tmp[M->fn++] = ptr[i].f;                                             \
        }                                                                      \
        unique[++j] = ptr[i];                                                  \
      } else                                                                   \
        unique[j].v += ptr[i].v;                                               \
    }                                                                          \
    off[M->fn] = ++j - sd * M->fn;                                             \
                                                                               \
    M->fld = tcalloc(ulong, M->fn);                                            \
    memcpy(M->fld, tmp, sizeof(ulong) * M->fn);                                \
                                                                               \
    M->adj_off = tcalloc(uint, M->fn + 1);                                     \
    memcpy(M->adj_off, off, (M->fn + 1) * sizeof(uint));                       \
                                                                               \
    M->adj_idx = tcalloc(uint, j - sd * M->fn);                                \
    M->adj_val = tcalloc(scalar, j - sd * M->fn);                              \
  } while (0)

  uint nnz = mijs->n;
  buffer_reserve(bfr, (sizeof(struct mij) + 2 * sizeof(ulong)) * (nnz + 1));
  ulong *tmp = (ulong *)bfr->ptr;

  if (type == CSC) {
    INIT_IDX_FLD(r, rn, rows);
    sarray_sort_2(struct mij, mijs->ptr, mijs->n, c, 1, r, 1, bfr);
  } else {
    INIT_IDX_FLD(c, cn, cols);
    sarray_sort_2(struct mij, mijs->ptr, mijs->n, r, 1, c, 1, bfr);
  }

  uint *off = (uint *)bfr->ptr;
  tmp = (ulong *)(off + nnz + 1);
  struct mij *unique = (struct mij *)(tmp + nnz + 1),
             *ptr = (struct mij *)mijs->ptr;

  off[0] = 0, unique[0] = ptr[0];
  uint nn;
  if (type == CSC) {
    INIT_ADJ_OFF(c, cn, cols);
    nn = M->cn;
  } else {
    INIT_ADJ_OFF(r, rn, rows);
    nn = M->rn;
  }

  uint nadj = 0, ndiag = 0, j = off[nn] + sd * nn, i = 0;
  if (sd != 0) {
    M->diag_idx = (uint *)tcalloc(uint, nn);
    M->diag_val = (scalar *)tcalloc(scalar, nn);
    for (; i < j; i++) {
      if (unique[i].r == unique[i].c) {
        M->diag_idx[ndiag] = unique[i].idx;
        M->diag_val[ndiag++] = unique[i].v;
      } else {
        M->adj_idx[nadj] = unique[i].idx;
        M->adj_val[nadj++] = unique[i].v;
      }
    }
  } else {
    M->diag_idx = NULL, M->diag_val = NULL;
    for (; i < j; i++)
      M->adj_idx[nadj] = unique[i].idx, M->adj_val[nadj++] = unique[i].v;
  }

#undef INIT_ADJ_OFF
#undef INIT_IDX_FLD

  return 0;
}

void par_csr_to_csc(struct par_mat *N, const struct par_mat *M, int diag,
                    struct crystal *cr, buffer *bfr) {
  assert(IS_CSR(M));

  slong *cols = tcalloc(slong, M->cn + M->rn);
  for (uint i = 0; i < M->cn; i++)
    cols[i] = -M->cols[i];
  for (uint i = 0; i < M->rn; i++)
    cols[M->cn + i] = M->rows[i];

  struct comm *c = &cr->comm;
  struct gs_data *gsh = gs_setup(cols, M->cn + M->rn, c, 0, gs_pairwise, 0);

  sint *owner = (sint *)cols;
  for (uint i = 0; i < M->cn; i++)
    owner[i] = -1;
  for (uint i = 0; i < M->rn; i++)
    owner[M->cn + i] = c->id;

  gs(owner, gs_int, gs_max, 0, gsh, bfr);
  gs_free(gsh);

  uint nnz = (M->rn > 0 ? M->adj_off[M->rn] + M->rn : 0);
  struct array mijs;
  array_init(struct mij, &mijs, nnz);

  struct mij m = {.r = 0, .c = 0, .idx = 0, .p = 0, .v = 0};
  for (uint i = 0; i < M->rn; i++) {
    m.r = M->rows[i];
    for (uint j = M->adj_off[i], je = M->adj_off[i + 1]; j < je; j++) {
      m.c = M->cols[M->adj_idx[j]];
      m.p = owner[M->adj_idx[j]];
      m.v = M->adj_val[j];
      array_cat(struct mij, &mijs, &m, 1);
    }
  }
  if (IS_DIAG(M)) {
    for (uint i = 0; i < M->rn; i++) {
      m.r = m.c = M->rows[i];
      m.v = M->diag_val[i];
      m.p = owner[M->diag_idx[i]];
      array_cat(struct mij, &mijs, &m, 1);
    }
  }
  free(cols);

  sarray_transfer(struct mij, &mijs, p, 0, cr);
  par_mat_setup(N, &mijs, 0, diag, bfr);
  array_free(&mijs);
}

void par_csc_to_csr(struct par_mat *N, const struct par_mat *M, int diag,
                    struct crystal *cr, buffer *bfr) {
  assert(IS_CSC(M) && !IS_DIAG(M));

  slong *rows = tcalloc(slong, M->rn + M->cn);
  for (uint i = 0; i < M->rn; i++)
    rows[i] = -M->rows[i];
  for (uint i = 0; i < M->cn; i++)
    rows[M->rn + i] = M->cols[i];

  struct comm *c = &cr->comm;
  struct gs_data *gsh = gs_setup(rows, M->rn + M->cn, c, 0, gs_pairwise, 0);

  sint *owner = (sint *)rows;
  for (uint i = 0; i < M->rn; i++)
    owner[i] = -1;
  for (uint i = 0; i < M->cn; i++)
    owner[M->rn + i] = c->id;

  gs(owner, gs_int, gs_max, 0, gsh, bfr);
  gs_free(gsh);

  uint nnz = (M->cn > 0 ? M->adj_off[M->cn] + M->cn : 0);
  struct array mijs;
  array_init(struct mij, &mijs, nnz);

  struct mij m = {.r = 0, .c = 0, .idx = 0, .p = 0, .v = 0};
  for (uint i = 0; i < M->cn; i++) {
    m.c = M->cols[i];
    for (uint j = M->adj_off[i], je = M->adj_off[i + 1]; j < je; j++) {
      m.r = M->rows[M->adj_idx[j]];
      m.p = owner[M->adj_idx[j]];
      m.v = M->adj_val[j];
      array_cat(struct mij, &mijs, &m, 1);
    }
  }
  free(rows);

  sarray_transfer(struct mij, &mijs, p, 0, cr);
  par_mat_setup(N, &mijs, 1, diag, bfr);
  array_free(&mijs);
}

void par_mat_print(struct par_mat *A) {
  if (IS_CSR(A)) {
    for (uint i = 0; i < A->rn; i++) {
      for (uint j = A->adj_off[i]; j < A->adj_off[i + 1]; j++)
        printf("%lld %lld %lf\n", A->rows[i], A->cols[A->adj_idx[j]],
               A->adj_val[j]);
      if (IS_DIAG(A))
        printf("%lld %lld %lf\n", A->rows[i], A->rows[i], A->diag_val[i]);
    }
  } else if (IS_CSC(A)) {
    for (uint i = 0; i < A->cn; i++) {
      for (uint j = A->adj_off[i]; j < A->adj_off[i + 1]; j++)
        printf("%lld %lld %lf\n", A->rows[A->adj_idx[j]], A->cols[i],
               A->adj_val[j]);
      if (IS_DIAG(A))
        printf("%lld %lld %lf\n", A->cols[i], A->cols[i], A->diag_val[i]);
    }
  }
  fflush(stdout);
}

void par_mat_dump(const char *name, struct par_mat *A, struct crystal *cr,
                  buffer *bfr) {
  struct array mijs;
  array_init(struct mij, &mijs, 1024);

  struct mij t = {.r = 0, .c = 0, .idx = 0, .p = 0, .v = 0};
  if (IS_CSR(A)) {
    for (uint i = 0; i < A->rn; i++) {
      t.r = A->rows[i];
      for (uint j = A->adj_off[i]; j < A->adj_off[i + 1]; j++) {
        t.c = A->cols[A->adj_idx[j]], t.v = A->adj_val[j];
        array_cat(struct mij, &mijs, &t, 1);
      }
      if (IS_DIAG(A)) {
        t.c = t.r, t.v = A->diag_val[i];
        array_cat(struct mij, &mijs, &t, 1);
      }
    }
  } else if (IS_CSC(A)) {
    for (uint i = 0; i < A->cn; i++) {
      t.c = A->cols[i];
      for (uint j = A->adj_off[i]; j < A->adj_off[i + 1]; j++) {
        t.r = A->rows[A->adj_idx[j]], t.v = A->adj_val[j];
        array_cat(struct mij, &mijs, &t, 1);
      }
      if (IS_DIAG(A)) {
        t.r = t.c, t.v = A->diag_val[i];
        array_cat(struct mij, &mijs, &t, 1);
      }
    }
  }

  sarray_transfer(struct mij, &mijs, p, 1, cr);
  sarray_sort_2(struct mij, mijs.ptr, mijs.n, r, 1, c, 1, bfr);

  struct comm *c = &cr->comm;
  if (c->id == 0) {
    FILE *fp = fopen(name, "w+");
    if (fp != NULL) {
      struct mij *pm = (struct mij *)mijs.ptr;
      for (uint i = 0; i < mijs.n; i++)
        fprintf(fp, "%llu %llu %.15lf\n", pm[i].r, pm[i].c, pm[i].v);
      fclose(fp);
    }
  }

  array_free(&mijs);
}

int par_mat_free(struct par_mat *A) {
  FREE(A, rows);
  FREE(A, cols);
  FREE(A, adj_off);
  FREE(A, adj_idx);
  FREE(A, adj_val);
  FREE(A, diag_idx);
  FREE(A, diag_val);
  A->rn = A->cn = 0;
  A->type = -1;

  return 0;
}

//------------------------------------------------------------------------------
// compress entries of a matrix
//
static int compress_mij(struct array *eij, struct array *entries, buffer *bfr) {
  eij->n = 0;
  if (entries->n == 0)
    return 1;

  sarray_sort_2(struct mij, entries->ptr, entries->n, r, 1, c, 1, bfr);
  struct mij *ptr = (struct mij *)entries->ptr;

  struct mij m;
  m.idx = 0;

  uint i = 0;
  while (i < entries->n) {
    m = ptr[i];
    uint j = i + 1;
    while (j < entries->n && ptr[j].r == ptr[i].r && ptr[j].c == ptr[i].c)
      m.v += ptr[j].v, j++;

    array_cat(struct mij, eij, &m, 1);
    i = j;
  }

  // Now make sure the row sum is zero
  struct mij *pe = (struct mij *)eij->ptr;
  i = 0;
  while (i < eij->n) {
    sint j = i, k = -1;
    scalar s = 0;
    while (j < eij->n && pe[j].r == pe[i].r) {
      if (pe[j].r == pe[j].c)
        k = j;
      else
        s += pe[j].v;
      j++;
    }
    assert(k >= 0);
    pe[k].v = -s;
    i = j;
  }
}

struct par_mat *par_csr_setup_con(const uint nelt, const ulong *eid,
                                  const slong *vtx, int nv, int sep,
                                  struct comm *c, struct crystal *cr,
                                  buffer *bfr) {
  struct array nbrs, eij;
  find_nbrs(&nbrs, eid, vtx, nelt, nv, cr, bfr);
  compress_nbrs(&eij, &nbrs, bfr);

  struct par_mat *M = tcalloc(struct par_mat, 1);
  par_csr_setup(M, &eij, sep, bfr);

  array_free(&eij), array_free(&nbrs);

  return M;
}

struct par_mat *par_csr_setup_ext(struct array *entries, int sep, buffer *bfr) {
  struct array eij;
  array_init(struct mij, &eij, 100);

  compress_mij(&eij, entries, bfr);

  struct par_mat *M = tcalloc(struct par_mat, 1);
  par_csr_setup(M, &eij, sep, bfr);

  array_free(&eij);

  return M;
}

struct par_mat *par_csc_setup_ext(struct array *entries, int sep, buffer *bfr) {
  struct array eij;
  array_init(struct mij, &eij, 100);

  compress_mij(&eij, entries, bfr);

  struct par_mat *M = tcalloc(struct par_mat, 1);
  par_csc_setup(M, &eij, sep, bfr);

  array_free(&eij);

  return M;
}

int IS_CSC(const struct par_mat *A) { return (A->type == CSC); }

int IS_CSR(const struct par_mat *A) { return (A->type == CSR); }

int IS_DIAG(const struct par_mat *A) { return (A->diag_val != NULL); }

struct gs_data *setup_Q(const struct par_mat *M, const struct comm *c,
                        buffer *bfr) {
  uint n;
  ulong *ids, *diag;
  if (IS_CSR(M))
    n = M->rn, ids = M->cols, diag = M->rows;
  else if (IS_CSC(M))
    n = M->cn, ids = M->rows, diag = M->cols;
  else {
    fprintf(stderr, "%s:%d Wrong matrix type !\n", __FILE__, __LINE__);
    exit(1);
  }

  assert(n == 0 || IS_DIAG(M));

  // Setup gs handle for the mat-vec
  uint nnz = n > 0 ? M->adj_off[n] + n : 0;
  buffer_reserve(bfr, nnz * sizeof(slong));
  slong *sids = (slong *)bfr->ptr;

  uint i, j = 0;
  for (i = 0; i < n; i++)
    for (j = M->adj_off[i]; j < M->adj_off[i + 1]; j++)
      sids[j] = -ids[M->adj_idx[j]];
  for (i = 0; i < n; i++)
    sids[j++] = diag[i];

  return gs_setup(sids, nnz, c, 0, gs_crystal_router, 0);
}

void mat_vec_csr(scalar *y, const scalar *x, const struct par_mat *M,
                 struct gs_data *gsh, scalar *buf, buffer *bfr) {
  assert(IS_CSR(M));
  assert(M->rn == 0 || IS_DIAG(M));

  uint n = M->rn, *Lp = M->adj_off, nnz = n > 0 ? Lp[n] : 0;
  uint i, j, je;
  for (i = 0; i < nnz; i++)
    buf[i] = 0.0; // Is this really necessary?
  for (i = 0, j = nnz; i < n; i++, j++)
    y[i] = buf[j] = x[i];

  gs(buf, gs_double, gs_add, 0, gsh, bfr);

  scalar *D = M->diag_val, *L = M->adj_val;
  for (i = 0; i < n; i++) {
    for (y[i] *= D[i], j = Lp[i], je = Lp[i + 1]; j != je; j++)
      y[i] += L[j] * buf[j];
  }
}

//------------------------------------------------------------------------------
// Dump and array in parallel
void par_arr_dump(const char *name, struct array *arr, struct crystal *cr,
                  buffer *bfr) {
  struct mij *ptr = arr->ptr;
  for (uint i = 0; i < arr->n; i++)
    ptr[i].p = 0;
  sarray_transfer(struct mij, arr, p, 0, cr);
  sarray_sort_2(struct mij, arr->ptr, arr->n, r, 1, c, 1, bfr);

  struct comm *c = &cr->comm;
  if (c->id == 0) {
    FILE *fp = fopen(name, "w+");
    if (fp != NULL) {
      struct mij *ptr = (struct mij *)arr->ptr;
      for (uint i = 0; i < arr->n; i++)
        fprintf(fp, "%llu %llu %.15lf\n", ptr[i].r, ptr[i].c, ptr[i].v);
      fclose(fp);
    }
  }
}

#undef CSC
#undef CSR
#undef FREE

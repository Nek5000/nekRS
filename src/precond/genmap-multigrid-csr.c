#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

#define ABS(i) ((i < 0) ? -i : i)

void csr_mat_setup(struct array *entries, struct comm *c, csr_mat *M_) {
  entry *ptr = entries->ptr;

  buffer buf;
  buffer_init(&buf, 1024);
  sarray_sort_2(entry, ptr, entries->n, r, 1, c, 1, &buf);
  buffer_free(&buf);

  uint st = 0, e = 0;
  sint diag;
  while (st < entries->n) {
    diag = -1;
    e = st;
    while (e < entries->n && ptr[st].r == ptr[e].r) {
      if (ptr[e].r == ptr[e].c)
        diag = e;
      else
        ptr[e].v = -1.0;
      e++;
    }
    assert(diag >= 0);
    ptr[diag].v = e - st - 1;
    st = e;
  }

  uint i = 0, j, n = 0;
  while (i < entries->n) {
    j = i + 1;
    while (j < entries->n && ptr[i].r == ptr[j].r)
      j++;
    i = j, n++;
  }

  GenmapMalloc(1, M_);
  csr_mat M = *M_;
  M->rn = n;

  slong out[2][1], bf[2][1], in = M->rn;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bf);
  M->row_start = out[0][0] + 1;

  GenmapMalloc(M->rn + 1, &M->row_off);

  if (n == 0)
    M->col = NULL, M->v = NULL, M->diag = NULL;
  else {
    GenmapMalloc(entries->n, &M->col), GenmapMalloc(entries->n, &M->v);
    GenmapMalloc(M->rn, &M->diag);
  }

  ptr = entries->ptr;
  uint rn = 0;
  for (i = 0; i < entries->n; i++) {
    M->col[i] = ptr[i].c, M->v[i] = ptr[i].v;
    if (ptr[i].r == ptr[i].c)
      M->diag[rn++] = ptr[i].v;
  }
  assert(rn == M->rn);

  M->row_off[0] = 0, i = 0;
  uint nn = 0;
  while (i < entries->n) {
    j = i + 1;
    while (j < entries->n && ptr[i].r == ptr[j].r)
      j++;
    i = M->row_off[++nn] = j;
  }
  assert(n == nn);
  assert(M->row_off[n] == entries->n);

  M->gsh = get_csr_top(M, c);
}

struct gs_data *get_csr_top(csr_mat M, struct comm *c) {
  const uint rn = M->rn;
  const uint n = M->row_off[rn];

  slong *ids;
  if (n > 0)
    GenmapMalloc(n, &ids);

  uint i, j;
  for (i = 0; i < rn; i++)
    for (j = M->row_off[i]; j < M->row_off[i + 1]; j++)
      if (M->row_start + i == M->col[j])
        ids[j] = M->col[j];
      else
        ids[j] = -M->col[j];

  struct gs_data *gsh = gs_setup(ids, n, c, 0, gs_pairwise, 0);

  if (n > 0)
    GenmapFree(ids);

  return gsh;
}

void csr_mat_gather(csr_mat M, struct gs_data *gsh, GenmapScalar *x,
                    GenmapScalar *buf, buffer *bfr) {
  ulong s = M->row_start;
  sint i, j;
  for (i = 0; i < M->rn; i++)
    for (j = M->row_off[i]; j < M->row_off[i + 1]; j++)
      if (M->col[j] == s + i)
        buf[j] = x[i];
      else
        buf[j] = 0.0;

  gs(buf, gs_scalar, gs_add, 0, gsh, bfr);
}

void csr_mat_apply(GenmapScalar *y, csr_mat M, GenmapScalar *x) {
  const uint rn = M->rn;
  if (rn == 0)
    return;

  const uint *offsets = M->row_off;
  const GenmapScalar *v = M->v;

  uint i, j, je;
  for (i = 0; i < rn; i++) {
    for (y[i] = 0.0, j = offsets[i], je = offsets[i + 1]; j < je; j++)
      y[i] += (*v++) * (*x++);
  }
}

void csr_mat_print(csr_mat M, struct comm *c) {
  const sint rn = M->rn;
  const uint *offsets = M->row_off;
  const GenmapScalar *v = M->v;
  const ulong *col = M->col;

  uint i, j, k;

  for (k = 0; k < c->np; k++) {
    comm_barrier(c);
    if (c->id == k) {
      for (i = 0; i < rn; i++) {
        for (j = offsets[i]; j < offsets[i + 1]; j++)
          fprintf(stderr, "(%lld,%lld) -> %.10lf\n", M->row_start + i, col[j],
                  v[j]);
      }
    }
    fflush(stderr);
  }
}

int csr_mat_free(csr_mat M) {
  if (M->col)
    GenmapFree(M->col);
  if (M->v)
    GenmapFree(M->v);
  if (M->diag)
    GenmapFree(M->diag);
  if (M->row_off)
    GenmapFree(M->row_off);
  if (M->gsh)
    gs_free(M->gsh);
  GenmapFree(M);

  return 0;
}

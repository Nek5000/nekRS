#include "mat.h"
#include "multigrid.h"
#include "parrsb-impl.h"

struct laplacian {
  int type, nv;
  uint nel;
  void *data;
};

//------------------------------------------------------------------------------
// Laplacian - as a `struct par_mat` in CSR mat
//
struct csr_laplacian {
  struct par_mat *M;
  struct gs_data *gsh;
  scalar *buf;
};

static void find_nbrs_rsb(struct array *arr, const struct rsb_element *elems,
                          const uint nelt, const int nv, const struct comm *c,
                          struct crystal *cr, buffer *buf) {
  slong out[2][1], bfr[2][1], in = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bfr);
  ulong eid = out[0][0] + 1;

  struct array vertices;
  array_init(struct nbr, &vertices, nelt * nv);

  struct nbr v;
  uint i, j;
  for (i = 0; i < nelt; i++) {
    v.r = eid++;
    for (j = 0; j < nv; j++) {
      v.c = elems[i].vertices[j], v.proc = v.c % c->np;
      array_cat(struct nbr, &vertices, &v, 1);
    }
  }

  sarray_transfer(struct nbr, &vertices, proc, 1, cr);

  sarray_sort(struct nbr, vertices.ptr, vertices.n, c, 1, buf);
  struct nbr *vptr = vertices.ptr;
  uint vn = vertices.n;

  // FIXME: Assumes quads or hexes
  struct nbr t;
  uint s = 0, e;
  array_init(struct nbr, arr, vertices.n * 10);
  while (s < vn) {
    e = s + 1;
    while (e < vn && vptr[s].c == vptr[e].c)
      e++;
    for (i = s; i < e; i++) {
      t = vptr[i];
      for (j = s; j < e; j++) {
        t.c = vptr[j].r;
        array_cat(struct nbr, arr, &t, 1);
      }
    }
    s = e;
  }

  sarray_transfer(struct nbr, arr, proc, 1, cr);
  array_free(&vertices);
}

static int par_csr_init(struct laplacian *l, const struct rsb_element *elems,
                        const uint nelt, const int nv, const struct comm *c,
                        buffer *bfr) {
  struct crystal cr;
  crystal_init(&cr, c);

  struct array nbrs, eij;
  find_nbrs_rsb(&nbrs, elems, nelt, nv, c, &cr, bfr);
  compress_nbrs(&eij, &nbrs, bfr);

  struct csr_laplacian *L = l->data = tcalloc(struct csr_laplacian, 1);
  struct par_mat *M = L->M = par_csr_setup_ext(&eij, 1, bfr);
  L->gsh = setup_Q(L->M, c, bfr);

  uint nnz = M->rn > 0 ? M->adj_off[M->rn] + M->rn : 0;
  L->buf = tcalloc(scalar, nnz);

  crystal_free(&cr);

  array_free(&nbrs);
  array_free(&eij);

  return 0;
}

static int par_csc_init(struct laplacian *l, const struct rsb_element *elems,
                        const uint nelt, const int nv, const struct comm *c,
                        buffer *bfr) {
  struct crystal cr;
  crystal_init(&cr, c);

  struct array nbrs, eij;
  find_nbrs_rsb(&nbrs, elems, nelt, nv, c, &cr, bfr);
  compress_nbrs(&eij, &nbrs, bfr);

  struct csr_laplacian *L = l->data = tcalloc(struct csr_laplacian, 1);
  struct par_mat *M = L->M = par_csc_setup_ext(&eij, 1, bfr);
  L->gsh = setup_Q(L->M, c, bfr);

  uint nnz = M->rn > 0 ? M->adj_off[M->rn] + M->rn : 0;
  L->buf = tcalloc(scalar, nnz);

  crystal_free(&cr);

  array_free(&nbrs);
  array_free(&eij);

  return 0;
}

static int par_csr(scalar *v, const struct laplacian *l, scalar *u,
                   buffer *bfr) {
  struct csr_laplacian *L = (struct csr_laplacian *)l->data;
  if (L != NULL) {
    mat_vec_csr(v, u, L->M, L->gsh, L->buf, bfr);
    return 0;
  }
  return 1;
}

static int par_csc(scalar *v, const struct laplacian *l, scalar *u,
                   buffer *bfr) {
#if 0
  struct csr_laplacian *L = (struct csr_laplacian *)l->data;
  if (L != NULL) {
    mat_vec_csr(v, u, L->M, L->gsh, L->buf, bfr);
    return 0;
  }
  return 1;
#else
  return 1;
#endif
}

static int par_csr_free(struct laplacian *l) {
  if (l->data != NULL) {
    struct csr_laplacian *L = (struct csr_laplacian *)l->data;
    par_mat_free(L->M), gs_free(L->gsh), free(L->buf);
    free(L);
  }
  return 0;
}

//------------------------------------------------------------------------------
// Laplacian - GS
//
struct gs_laplacian {
  scalar *diag, *u;
  struct gs_data *gsh;
};

static int gs_weighted_init(struct laplacian *l, struct rsb_element *elems,
                            uint lelt, int nv, struct comm *c, buffer *buf) {

  uint npts = nv * lelt;
  slong *vertices = tcalloc(slong, npts);
  uint i, j;
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++)
      vertices[i * nv + j] = elems[i].vertices[j];

  struct gs_laplacian *gl = l->data = tcalloc(struct gs_laplacian, 1);
  gl->u = tcalloc(scalar, npts);
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++)
      gl->u[nv * i + j] = 1.0;

  gl->gsh = gs_setup(vertices, npts, c, 0, gs_crystal_router, 0);
  gs(gl->u, gs_double, gs_add, 0, gl->gsh, buf);

  gl->diag = tcalloc(scalar, lelt);
  for (i = 0; i < lelt; i++) {
    gl->diag[i] = 0.0;
    for (j = 0; j < nv; j++)
      gl->diag[i] += gl->u[nv * i + j];
  }

  if (vertices != NULL)
    free(vertices);

  return 0;
}

static int gs_weighted(scalar *v, struct laplacian *l, scalar *u, buffer *buf) {
  uint lelt = l->nel;
  int nv = l->nv;

  struct gs_laplacian *gl = l->data;

  uint i, j;
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++)
      gl->u[nv * i + j] = u[i];

  gs(gl->u, gs_double, gs_add, 0, gl->gsh, buf);

  for (i = 0; i < lelt; i++) {
    v[i] = gl->diag[i] * u[i];
    for (j = 0; j < nv; j++)
      v[i] -= gl->u[nv * i + j];
  }

  return 0;
}

static int gs_weighted_free(struct laplacian *l) {
  struct gs_laplacian *gl = l->data;
  if (gl->u != NULL)
    free(gl->u);
  if (gl->diag != NULL)
    free(gl->diag);
  gs_free(gl->gsh);
  free(l->data);
  return 0;
}

//------------------------------------------------------------------------------
// Laplacian
//
struct laplacian *laplacian_init(struct rsb_element *elems, uint nel, int nv,
                                 int type, struct comm *c, buffer *buf) {
  struct laplacian *l = tcalloc(struct laplacian, 1);
  l->type = type;
  l->nv = nv;
  l->nel = nel;

  if (type & CSR)
    par_csr_init(l, elems, nel, nv, c, buf);
  else if (type & CSC)
    par_csc_init(l, elems, nel, nv, c, buf);
  else if (type & GS)
    gs_weighted_init(l, elems, nel, nv, c, buf);
  else
    return NULL;

  return l;
}

int laplacian(scalar *v, struct laplacian *l, scalar *u, buffer *buf) {
  if (l->type & CSR)
    par_csr(v, l, u, buf);
  else if (l->type & CSC)
    par_csc(v, l, u, buf);
  else if (l->type & GS)
    gs_weighted(v, l, u, buf);
  else
    return 1;

  return 0;
}

void laplacian_free(struct laplacian *l) {
  if (l) {
    if (l->type & CSR)
      par_csr_free(l);
    else if (l->type & GS) {
      gs_weighted_free(l);
    }
    free(l);
  }
}

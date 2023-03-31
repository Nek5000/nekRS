#include "coarse-impl.h"
#include "metrics.h"
#include <float.h>

//------------------------------------------------------------------------------
// Setup coarse grid system. Initial dumb API.
//
// Number rows, local first then interface. Returns global number of local
// elements.
struct rcb_t {
  uint i, s;
  double coord[3];
  slong vtx[8];
};

static void nmbr_local_rcb(struct array *a, uint s, uint e, const unsigned nc,
                           const unsigned ndim, const unsigned level,
                           struct comm *c, buffer *bfr) {
  sint size = e - s;
  if (size <= 1)
    return;

  double max[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX},
         min[3] = {DBL_MAX, DBL_MAX, DBL_MAX};

  struct rcb_t *pa = (struct rcb_t *)a->ptr;
  for (uint i = s; i < e; i++) {
    for (int j = 0; j < ndim; j++) {
      if (pa[i].coord[j] < min[j])
        min[j] = pa[i].coord[j];
      if (pa[i].coord[j] > max[j])
        max[j] = pa[i].coord[j];
    }
  }

  double len = max[0] - min[0];
  int axis = 0;
  for (int j = 1; j < ndim; j++) {
    if (max[j] - min[j] > len)
      axis = j, len = max[j] - min[j];
  }

  struct rcb_t *ps = pa + s;
  switch (axis) {
  case 0:
    sarray_sort(struct rcb_t, ps, size, coord[0], 3, bfr);
    break;
  case 1:
    sarray_sort(struct rcb_t, ps, size, coord[1], 3, bfr);
    break;
  case 2:
    sarray_sort(struct rcb_t, ps, size, coord[2], 3, bfr);
    break;
  default:
    break;
  }

  // Number the elements in the interface
  uint npts = size * nc;
  slong *vtx = tcalloc(slong, npts);
  for (uint i = s, k = 0; i < e; i++) {
    for (int j = 0; j < nc; j++, k++)
      vtx[k] = pa[i].vtx[j];
  }

  struct gs_data *gsh = gs_setup(vtx, npts, c, 0, gs_pairwise, 0);

  sint *dof = tcalloc(sint, npts);
  uint mid = (s + e) / 2;
  for (uint i = mid, k = (mid - s) * nc; i < e; i++) {
    for (int j = 0; j < nc; j++, k++)
      dof[k] = 1;
  }

  gs(dof, gs_int, gs_add, 0, gsh, bfr);

  for (uint i = mid, k = (mid - s) * nc; i < e; i++) {
    for (int j = 0; j < nc; j++, k++)
      dof[k] = 0;
  }

  gs(dof, gs_int, gs_add, 0, gsh, bfr);

  for (uint i = s, k = 0; i < e; i++, k++) {
    for (int j = 0; j < nc; j++) {
      if (dof[k * nc + j] > 0 && pa[i].s == INT_MAX) {
        pa[i].s = level;
        break;
      }
    }
  }

  gs_free(gsh);
  free(dof), free(vtx);

  nmbr_local_rcb(a, s, mid, nc, ndim, level + 1, c, bfr);
  nmbr_local_rcb(a, mid, e, nc, ndim, level + 1, c, bfr);
}

// Number the DOFs internal first, faces second and all the rest (wire basket)
// next. This keeps zeros as is and renumber the positive entries in `ids`
// array.
static void number_dual_graph_dofs(ulong *dofs, struct coarse *crs, uint n,
                                   const slong *ids, uint nelt, unsigned ndim,
                                   const scalar *coord, buffer *bfr) {
  int nnz = (n > 0);
  struct comm c;
  comm_split(&crs->c, nnz, crs->c.id, &c);

  unsigned nc = n / nelt;
  uint i, j;
  if (nnz) {
    sint *dof = tcalloc(sint, n);
    int level = 1;
    while (c.np > 1) {
      struct gs_data *gsh = gs_setup(ids, n, &c, 0, gs_pairwise, 0);

      int bin = (c.id >= (c.np + 1) / 2);
      for (i = 0; i < n; i++)
        dof[i] = bin;

      gs(dof, gs_int, gs_add, 0, gsh, bfr);

      if (bin == 1) {
        for (i = 0; i < n; i++)
          dof[i] = 0;
      }

      gs(dof, gs_int, gs_add, 0, gsh, bfr);

      for (i = 0; i < nelt; i++) {
        for (j = 0; j < nc; j++) {
          if (dof[i * nc + j] > 0 && !dofs[i]) {
            dofs[i] = level;
            break;
          }
        }
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

  for (i = crs->n[0] = crs->n[1] = 0; i < nelt; i++) {
    if (dofs[i] > 0)
      crs->n[1]++;
    else
      crs->n[0]++;
  }

  slong in[2] = {crs->n[0], crs->n[1]}, out[2][2], wrk[2][2];
  comm_scan(out, &crs->c, gs_long, gs_add, in, 2, wrk);
  crs->s[0] = out[0][0] + 1, crs->ng[0] = out[1][0];
  crs->s[1] = out[0][1] + 1, crs->ng[1] = out[1][1];

  struct array local;
  array_init(struct rcb_t, &local, crs->n[0]);

  struct rcb_t t = {.s = INT_MAX};
  ulong s = crs->ng[0] + crs->s[1];
  for (uint i = 0; i < nelt; i++) {
    if (dofs[i] > 0)
      dofs[i] = s++;
    else {
      t.i = i;
      memcpy(t.coord, &coord[i * ndim], ndim * sizeof(scalar));
      memcpy(t.vtx, &ids[i * nc], nc * sizeof(slong));
      array_cat(struct rcb_t, &local, &t, 1);
    }
  }

  if (local.n > 0) {
    nmbr_local_rcb(&local, 0, local.n, nc, ndim, 1, &c, bfr);
    sarray_sort(struct rcb_t, local.ptr, local.n, s, 0, bfr);
    struct rcb_t *pl = (struct rcb_t *)local.ptr;
    ulong s = crs->s[0];
    for (sint i = local.n - 1; i >= 0; i--)
      dofs[pl[i].i] = s++;
  }

  comm_free(&c);
  array_free(&local);
}

struct coarse *coarse_setup(unsigned n, unsigned nc, const long long *vl,
                            const scalar *coord, unsigned null_space,
                            unsigned type, struct comm *c) {
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

  uint size = n * nc;
  slong *tid = tcalloc(slong, size);
  for (uint i = 0; i < size; i++)
    tid[i] = vl[i];

  ulong *nid = tcalloc(ulong, n);
  unsigned ndim = (nc == 8) ? 3 : 2;
  number_dual_graph_dofs(nid, crs, size, tid, crs->un, ndim, coord, &crs->bfr);

  // Find unique ids and user vector to compressed vector mapping.
  // In the case of dual-graph Laplacian, all the ids are unique.
  // But here we arrange them in the sorted order.
  ulong *uid = tcalloc(ulong, n);
  crs->u2c = tcalloc(sint, n);
  crs->cn = unique_ids(crs->u2c, uid, crs->un, nid, &crs->bfr);
  crs->an = crs->cn;

  struct crystal cr;
  crystal_init(&cr, &crs->c);

  struct array nbrs, eij;
  find_nbrs(&nbrs, nid, tid, n, nc, &cr, &crs->bfr);
  // Convert `struct nbr` -> `struct mij` and compress entries which share the
  // same (r, c) values. Set the diagonal element to have zero row sum
  compress_nbrs(&eij, &nbrs, &crs->bfr);
  array_free(&nbrs);

  switch (type) {
  case 0:
    schur_setup(crs, &eij, &cr, &crs->bfr);
    break;
  default:
    break;
  }

  array_free(&eij), crystal_free(&cr);
  free(tid), free(nid), free(uid);

  return crs;
}

void coarse_solve(scalar *x, struct coarse *crs, scalar *b, scalar tol) {
  metric_init();

  scalar *rhs = tcalloc(scalar, 2 * crs->an), *xx = rhs + crs->an;
  for (uint i = 0; i < crs->un; i++) {
    if (crs->u2c[i] >= 0)
      rhs[crs->u2c[i]] += b[i];
  }

  switch (crs->type) {
  case 0:
    schur_solve(xx, crs, rhs, tol, &crs->bfr);
    break;
  default:
    break;
  }

  for (uint i = 0; i < crs->un; i++) {
    if (crs->u2c[i] >= 0)
      x[i] = xx[crs->u2c[i]];
  }
  free(rhs);

  metric_push_level();
  metric_crs_print(&crs->c, 1);
  metric_finalize();
}

void coarse_free(struct coarse *crs) {
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
    comm_free(&crs->c), buffer_free(&crs->bfr);
    free(crs), crs = NULL;
  }
}

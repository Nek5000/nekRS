#include "ilu.h"
#include <math.h>

#define CSC 0
#define CSR 1

//=============================================================================
// ILU levels
//
// Currently there are two methods of finding levels
//   1. Based on final element distribution among processors (dst_lvls)
//   2. Based on RSB levels while partitioning (rsb_lvls)
struct key_t {
  ulong e;
  uint p;
};

struct e2n_t {
  ulong e, n;
};

struct request_t {
  ulong r;
  uint p, o;
};

static int find_unique_nbrs(struct array *e2nm, uint n, int nv,
                            const ulong *ids, const slong *vtx,
                            struct crystal *cr, buffer *bfr) {
  struct array nbrs;
  find_nbrs(&nbrs, ids, vtx, n, nv, cr, bfr);

  array_init(struct e2n_t, e2nm, n * 10);
  if (nbrs.n > 0) {
    sarray_sort_2(struct nbr, nbrs.ptr, nbrs.n, r, 1, c, 1, bfr);
    struct nbr *pn = (struct nbr *)nbrs.ptr;

    struct e2n_t en;
    uint i, j;
    for (i = 1, j = 0; i < nbrs.n; i++) {
      if ((pn[i].r != pn[j].r) || (pn[i].c != pn[j].c)) {
        en.e = pn[j].r, en.n = pn[j].c;
        array_cat(struct e2n_t, e2nm, &en, 1);
        j = i;
      }
    }
    en.e = pn[j].r, en.n = pn[j].c;
    array_cat(struct e2n_t, e2nm, &en, 1);
    sarray_sort_2(struct e2n_t, e2nm->ptr, e2nm->n, e, 1, n, 1, bfr);
  }
  array_free(&nbrs);

  return 0;
}

static int local_dof(const ulong *rows, const ulong I, const uint n) {
  for (uint i = 0; i < n; i++)
    if (rows[i] == I)
      return i;
  return n;
}

// Fill dofs array with unique dofs found in this processr
static int update_keys(struct array *keys, struct array *nbrs, const uint ln,
                       const ulong *lids, struct crystal *cr, buffer *bfr) {
  uint i, j;
  struct array temp, rqst;
  array_init(struct request_t, &temp, nbrs->n);
  array_init(struct request_t, &rqst, nbrs->n);

  struct comm *c = &cr->comm;
  struct e2n_t *pn = (struct e2n_t *)nbrs->ptr;
  struct request_t t;
  for (i = 0; i < nbrs->n; i++) {
    t.r = pn[i].n, t.p = t.r % c->np;
    t.o = (local_dof(lids, t.r, ln) < ln);
    array_cat(struct request_t, &temp, &t, 1);
  }

  struct request_t *pt = (struct request_t *)temp.ptr;
  if (temp.n > 0) {
    sarray_sort(struct request_t, temp.ptr, temp.n, r, 1, bfr);
    for (i = 1, j = 0; i < temp.n; i++) {
      if (pt[i].r != pt[j].r) {
        array_cat(struct request_t, &rqst, &pt[j], 1);
        j = i;
      }
    }
    array_cat(struct request_t, &rqst, &pt[j], 1);
  }

  sarray_transfer(struct request_t, &rqst, p, 1, cr);
  sarray_sort_2(struct request_t, rqst.ptr, rqst.n, r, 1, o, 0, bfr);

  struct request_t *pr = (struct request_t *)rqst.ptr;
  if (rqst.n > 0) {
    for (i = 1, j = 0; i < rqst.n; i++) {
      if (pr[i].r != pr[j].r) {
        // owner for dof j, j + 1, ... i - 1 is pr[i - 1].p
        assert(pr[i - 1].o == 1);
        for (; j < i; j++)
          pr[j].o = pr[i - 1].p;
        // j = i at the end
      }
    }
    assert(pr[i - 1].o == 1);
    for (; j < i; j++)
      pr[j].o = pr[i - 1].p;
  }

  sarray_transfer(struct request_t, &rqst, o, 0, cr);
  sarray_sort_2(struct request_t, rqst.ptr, rqst.n, r, 1, p, 0, bfr);

  // All the requests are forwarded correctly. Send the data back
  // to the requesting processors. Note that the requests are unique.
  struct key_t *pk = (struct key_t *)keys->ptr;
  pr = (struct request_t *)rqst.ptr;
  temp.n = 0;
  for (i = j = 0; i < rqst.n; i++) {
    while (pk[j].e < pr[i].r)
      j++;
    // Sanity check
    assert(pk[j].e == pr[i].r);
    t.o = pr[i].p;
    for (uint k = j; k < keys->n && pk[k].e == pk[j].e; k++) {
      t.r = pk[k].e, t.p = pk[k].p;
      array_cat(struct request_t, &temp, &t, 1);
    }
  }
  array_free(&rqst);

  sarray_transfer(struct request_t, &temp, o, 0, cr);
  sarray_sort_2(struct request_t, temp.ptr, temp.n, r, 1, p, 0, bfr);

  // Update the keys array. Update here is a complete rewrite.
  struct array keyt;
  array_init(struct key_t, &keyt, temp.n);

  struct key_t s;
  pt = (struct request_t *)temp.ptr;
  for (i = 0; i < ln; i++) {
    ulong e = lids[i];
    // Find `e` in the nbrs array
    for (j = 0; j < nbrs->n && pn[j].e < e; j++)
      ;
    assert(j < nbrs->n && pn[j].e == e);
    // Now go through all the neighbors and update the keys
    for (; j < nbrs->n && pn[j].e == e; j++) {
      ulong n = pn[j].n;
      // find the key of `n` in temp
      uint k = 0;
      for (; k < temp.n && pt[k].r < n; k++)
        ;
      assert(k < temp.n && pt[k].r == n);
      for (; k < temp.n && pt[k].r == n; k++) {
        s.e = e, s.p = pt[k].p;
        array_cat(struct key_t, &keyt, &s, 1);
      }
    }
  }
  array_free(&temp);

  keys->n = 0;
  if (keyt.n > 0) {
    sarray_sort_2(struct key_t, keyt.ptr, keyt.n, e, 1, p, 0, bfr);
    pk = (struct key_t *)keyt.ptr;
    for (i = 1, j = 0; i < keyt.n; i++) {
      if ((pk[i].e != pk[j].e) || (pk[i].p != pk[j].p)) {
        array_cat(struct key_t, keys, &pk[j], 1);
        j = i;
      }
    }
    array_cat(struct key_t, keys, &pk[j], 1);
  }

  array_free(&keyt);

  return 0;
}

// This routine will update `lvl_n`, `lvl_off` and `lvl_ids` with the DOF
// belongig to current level. In the process, it will remove the DOFs and their
// connectivity from ids, and vtx arrays. `n` will be adjusted to reflect
// changes.
static int dst_lvls_aux(int *lvl_n, uint *lvl_off, uint *lvl_owner,
                        ulong *lvl_ids, uint *n, ulong *ids, slong *vtx, int nv,
                        struct array *keys, struct comm *c, int verbose) {
  // Find the min key size locally.
  uint i, j, k;
  sint min = INT_MAX;
  struct key_t *pk = (struct key_t *)keys->ptr;
  if (keys->n > 0) {
    for (i = 1, j = 0; i < keys->n; i++) {
      if (pk[i].e != pk[j].e) {
        // Different element, update min key size if required
        min = (min > i - j ? i - j : min);
        j = i;
      }
    }
    min = (min > i - j ? i - j : min);
  }

  sint buf[2];
  comm_allreduce(c, gs_int, gs_min, &min, 1, buf);
  if (min == INT_MAX)
    return 0;

  int lvl = *lvl_n;
  uint off = lvl_off[lvl];
  if (keys->n > 0) {
    for (i = 1, j = 0; i < keys->n; i++) {
      if (pk[i].e != pk[j].e) {
        if (i - j == min)
          lvl_ids[off] = pk[j].e, lvl_owner[off] = pk[i - 1].p, off++;
        j = i;
      }
    }
    if (i - j == min)
      lvl_ids[off] = pk[j].e, lvl_owner[off] = pk[i - 1].p, off++;
  }

  assert(lvl < 50);
  lvl++, lvl_off[lvl] = off;
  if (verbose > 1) {
    printf("id: %d |key| = %d lvl = %d size = %u\n", c->id, min, lvl,
           lvl_off[lvl] - lvl_off[lvl - 1]);
    fflush(stdout);
  }

  // Now we have to update ids and vtx. This can be done in place.
  for (i = lvl_off[lvl - 1], j = 0, k = 0; i < lvl_off[lvl]; i++, j++) {
    for (; j < *n && ids[j] < lvl_ids[i]; j++, k++) {
      ids[k] = ids[j];
      for (int v = 0; v < nv; v++)
        vtx[k * nv + v] = vtx[j * nv + v];
    }
    assert(j < *n && ids[j] == lvl_ids[i]);
  }
  for (; j < *n; j++, k++) {
    ids[k] = ids[j];
    for (int v = 0; v < nv; v++)
      vtx[k * nv + v] = vtx[j * nv + v];
  }

  *n -= lvl_off[lvl] - lvl_off[lvl - 1], *lvl_n = lvl;

  return 0;
}

static int dst_lvls(uint *lvl_off, uint *lvl_owner, ulong *lvl_ids,
                    const uint n_, const int nv, const ulong *ids_,
                    const slong *vtx_, struct crystal *cr, int verbose,
                    buffer *bfr) {
  // Copy ids and vtx since we are going to modify them
  uint n = n_;
  ulong *ids = tcalloc(ulong, n);
  slong *vtx = tcalloc(slong, n * nv);
  for (uint i = 0, j = 0; i < n; i++) {
    ids[i] = ids_[i];
    for (int v = 0; v < nv; v++, j++)
      vtx[j] = vtx_[j];
  }

  struct comm *c = &cr->comm;

  // Initialize keys: set key of each dof to the current MPI rank.
  // keys array should has unique entries and should be sorted first
  // by .e and then by .p.
  struct array keys;
  array_init(struct key_t, &keys, n);
  struct key_t e2p = {.e = 0, .p = c->id};
  for (uint i = 0; i < n; i++) {
    e2p.e = ids[i];
    array_cat(struct key_t, &keys, &e2p, 1);
  }
  sarray_sort_2(struct key_t, keys.ptr, keys.n, e, 1, p, 0, bfr);

  slong ng = n, buf[2];
  comm_allreduce(c, gs_long, gs_add, &ng, 1, buf);

  int nlvls = 0;
  struct array nbrs;
  while (ng > 0) {
    // Find unique neighbors of a DOF. DOF is a neighbor of itself.
    find_unique_nbrs(&nbrs, n, nv, ids, vtx, cr, bfr);

    // Send and receive key to/from neighbors. We forward all the requests
    // for the key of a DOF to the processor that owns the DOF and then that
    // processor takes care of the request. To do that, we first find all the
    // unique requests.
    update_keys(&keys, &nbrs, n, ids, cr, bfr);

    // Find the min key size
    // Add all the dofs with key size equal to min key size to current level
    // Update ids and vtx by removing the dofs with min key size
    dst_lvls_aux(&nlvls, lvl_off, lvl_owner, lvl_ids, &n, ids, vtx, nv, &keys,
                 c, verbose);

    ng = n;
    comm_allreduce(c, gs_long, gs_add, &ng, 1, buf);
    if (verbose > 1) {
      if (c->id == 0)
        printf("lvl = %d ng = %lld\n", nlvls, ng);
      fflush(stdout);
    }
    array_free(&nbrs);
  }

  free(ids), free(vtx);

  return nlvls;
}

static int rsb_lvls(uint *lvl_off, uint *lvl_owner, ulong *lvl_ids,
                    const uint n, const int nv, const ulong *ids,
                    const slong *vtx, struct comm *ci, int verbose,
                    buffer *bfr) {
  slong ng = n, buf[2];
  comm_allreduce(ci, gs_long, gs_add, &ng, 1, buf);

  // What we are going to do is identify the elements in the interface at each
  // level. These elements constitute the level of ILU. Owner of the element is
  // the processor which at least own a single vertex (possibly duplicated) of
  // the element.

  uint size = n * nv;
  sint *in = tcalloc(sint, size);
  sint *lvl = tcalloc(sint, n);
  sint *owner = tcalloc(sint, n);
  if (owner == NULL || lvl == NULL || in == NULL) {
    fprintf(stderr, "Failed to allocate lvl, owner or in !\n");
    exit(1);
  }

  struct comm c, t;
  comm_dup(&c, ci);

  uint i;
  sint nlvls = 1, j;
  while (c.np > 1) {
    struct gs_data *gsh = gs_setup(vtx, size, &c, 0, gs_pairwise, 0);

    int bin = (c.id >= (c.np + 1) / 2);
    for (i = 0; i < size; i++)
      in[i] = bin;

    gs(in, gs_int, gs_max, 0, gsh, bfr);

    if (bin == 1) {
      for (i = 0; i < size; i++)
        in[i] = 0;
    }

    gs(in, gs_int, gs_max, 0, gsh, bfr);

    sint ownr = 0;
    for (i = 0; i < n; i++) {
      for (j = 0; j < nv; j++) {
        if (in[i * nv + j] > 0) {
          if (lvl[i] == 0) {
            lvl[i] = nlvls;
            ownr = ci->id + 1;
          }
          break;
        }
      }
    }

    comm_allreduce(&c, gs_int, gs_max, &ownr, 1, buf);

    for (i = 0; i < n; i++) {
      if (lvl[i] == nlvls)
        owner[i] = ownr - 1;
    }

    nlvls++;

    gs_free(gsh);
    comm_split(&c, bin, c.id, &t), comm_free(&c);
    comm_dup(&c, &t), comm_free(&t);
  }
  comm_free(&c);

  int rem = 0;
  for (uint i = 0; i < n; i++) {
    if (lvl[i] == 0) {
      lvl[i] = nlvls;
      owner[i] = ci->id;
      rem = 1;
    }
  }
  nlvls += rem;
  comm_allreduce(ci, gs_int, gs_max, &nlvls, 1, buf);

  // Reverse the level numbers
  for (uint i = 0; i < n; i++)
    lvl[i] = nlvls - lvl[i];

  struct linfo_t {
    uint lvl, owner;
    ulong id;
  };

  struct array linfos;
  array_init(struct linfo_t, &linfos, n);

  struct linfo_t linfo = {.lvl = 0, .owner = 0, .id = 0};
  for (uint i = 0; i < n; i++) {
    linfo.lvl = lvl[i], linfo.owner = owner[i], linfo.id = ids[i];
    array_cat(struct linfo_t, &linfos, &linfo, 1);
  }
  sarray_sort(struct linfo_t, linfos.ptr, linfos.n, lvl, 0, bfr);

  if (linfos.n > 0) {
    struct linfo_t *pl = (struct linfo_t *)linfos.ptr;
    for (uint l = 0, i = 0; l < nlvls; l++) {
      for (; i < linfos.n && pl[i].lvl == l; i++)
        lvl_ids[i] = pl[i].id, lvl_owner[i] = pl[i].owner;
      lvl_off[l + 1] = i;
    }
  }

  array_free(&linfos);
  free(owner), free(lvl), free(in);

  return nlvls;
}

static int find_lvls(uint *lvl_off, uint *lvl_owner, ulong *lvl_ids,
                     const uint n, const int nv, const ulong *ids,
                     const slong *vtx, int type, struct crystal *cr,
                     int verbose, buffer *bfr) {
  int nlvls = 0;
  switch (type) {
  case 0:
    nlvls = dst_lvls(lvl_off, lvl_owner, lvl_ids, n, nv, ids, vtx, cr, verbose,
                     bfr);
    break;
  case 1:
    nlvls = rsb_lvls(lvl_off, lvl_owner, lvl_ids, n, nv, ids, vtx, &cr->comm,
                     verbose, bfr);
    break;
  default:
    break;
  }
  return nlvls;
}

//=============================================================================
// ILU
//
struct ilu {
  int pivot, verbose;
  // 1st dropping rule: An entry a_ij is dropped abs(a_ij) < tol
  scalar tol;
  // 2nd dropping rule: Entries are dropped so that total nnz per row/col < p
  uint nnz_per_row;

  // Calculated values internal to ILU
  uint nlvls, *lvl_off;
  ulong *perm;
  struct par_mat A, L, U;
  struct crystal cr;
};

//=============================================================================
// ILU(0)
//
static int ilu0_get_rows(struct par_mat *E, int lvl, uint *lvl_off,
                         struct par_mat *A, struct crystal *cr, buffer *bfr) {
  struct owner {
    ulong ri;
    uint rp, p;
  };

  assert(IS_CSR(A) && !IS_DIAG(A));

  struct array owners, requests;
  array_init(struct owner, &owners, A->rn * 30);
  array_init(struct owner, &requests, A->rn * 30);

  struct comm *c = &cr->comm;
  struct owner t;
  for (uint i = lvl_off[lvl - 1]; i < lvl_off[lvl]; i++) {
    ulong I = A->rows[i];
    for (uint j = A->adj_off[i];
         j < A->adj_off[i + 1] && A->cols[A->adj_idx[j]] < I; j++) {
      t.ri = A->cols[A->adj_idx[j]], t.rp = c->np, t.p = t.ri % c->np;
      array_cat(struct owner, &owners, &t, 1);
    }
  }

  for (uint i = lvl_off[0]; i < lvl_off[lvl]; i++) {
    t.ri = A->rows[i], t.rp = c->id, t.p = t.ri % c->np;
    array_cat(struct owner, &owners, &t, 1);
  }

  sarray_sort_2(struct owner, owners.ptr, owners.n, ri, 1, rp, 0, bfr);
  struct owner *ptr = (struct owner *)owners.ptr;
  uint i, j;
  for (i = 0; i < owners.n; i = j) {
    for (j = i + 1; j < owners.n && ptr[j].ri == ptr[i].ri; j++)
      ;
    array_cat(struct owner, &requests, &ptr[i], 1);
  }
  array_free(&owners);

  // Match row ids and set `p` to the original processor
  sarray_transfer(struct owner, &requests, p, 1, cr);

  // Set rp to the owner
  sarray_sort_2(struct owner, requests.ptr, requests.n, ri, 1, rp, 0, bfr);
  ptr = (struct owner *)requests.ptr;
  for (i = 0; i < requests.n; i = j) {
    assert(ptr[i].rp < c->np);
    for (j = i + 1; j < requests.n && ptr[j].ri == ptr[i].ri; j++) {
      assert(ptr[j].rp == c->np);
      ptr[j].rp = ptr[i].rp;
    }
  }

  // Forward requests to the owner processor
  sarray_transfer(struct owner, &requests, rp, 0, cr);

  sarray_sort_2(struct owner, requests.ptr, requests.n, ri, 1, p, 0, bfr);
  ptr = (struct owner *)requests.ptr;

  struct array sends;
  array_init(struct mij, &sends, A->rn * 30);

  for (i = 0; i < requests.n; i = j) {
    ulong ri = ptr[i].ri;
    uint ro = local_dof(A->rows, ri, A->rn);
    assert(ro < A->rn);
    for (j = i; j < requests.n && ptr[j].ri == ri; j++) {
      // No need to send to owner
      if (ptr[j].p != c->id) {
        // copy_row(&sends, ro, ptr[j].p, A);
        struct mij m = {.r = A->rows[ro], .idx = 0, .p = ptr[j].p};
        for (uint k = A->adj_off[ro], ke = A->adj_off[ro + 1]; k < ke; k++) {
          m.c = A->cols[A->adj_idx[k]], m.v = A->adj_val[k];
          array_cat(struct mij, &sends, &m, 1);
        }
      }
    }
  }
  array_free(&requests);

  sarray_transfer(struct mij, &sends, p, 1, cr);
  par_csr_setup(E, &sends, 0, bfr);
  array_free(&sends);

  return 0;
}

static void ilu0_update_row(const uint io, const uint k, struct par_mat *A,
                            struct par_mat *E, int verbose, int lvl) {
  uint *off = A->adj_off, *idx = A->adj_idx;
  uint *koff = A->adj_off, *kidx = A->adj_idx;
  ulong *cols = A->cols, *kcols = A->cols;
  scalar *val = A->adj_val, *kval = A->adj_val;

  const ulong K = cols[idx[k]];
  const ulong I = A->rows[io];

  // Find offsets of K in A
  sint ko = -1;
  uint j;
  for (j = 0; j < A->rn; j++) {
    if (A->rows[j] == K) {
      ko = j;
      break;
    }
  }

  // Search in E if K is not found in A
  if (ko == -1 && E != NULL) {
    koff = E->adj_off, kidx = E->adj_idx;
    kval = E->adj_val, kcols = E->cols;
    for (j = 0; j < E->rn; j++) {
      if (E->rows[j] == K) {
        ko = j;
        break;
      }
    }
  }

  // Oops, K is no where to be found
  if (ko == -1) {
    fprintf(stderr, "%s:%d lvl = %d, k = %u ko = %d\n", __FILE__, __LINE__, lvl,
            k, ko);
    exit(1);
  }

  // Calculate a_ik = a_ik / a_kk
  scalar a_kk = 0;
  for (j = koff[ko]; j < koff[ko + 1]; j++) {
    if (kcols[kidx[j]] == K) {
      a_kk = kval[j];
      break;
    }
  }

  if (fabs(a_kk) < 1e-10) {
    fprintf(stderr, "%s:%d ilu0: Diagonal is zero ! k = %llu\n", __FILE__,
            __LINE__, K);
    exit(1);
  }

  // cols[idx[k]] = K and val[k] = a_ik
  scalar a_ik = val[k] / a_kk;
  if (verbose) {
    printf("a_kk = %lf a_ik = %lf a_ik/a_kk = %lf\n", a_kk, val[j], a_ik);
    fflush(stdout);
  }
  val[k] = a_ik;

  uint kj;
  scalar a_kj;
  for (j = k + 1; j < off[io + 1]; j++) {
    for (kj = koff[ko]; kj < koff[ko + 1] && kcols[kidx[kj]] < cols[idx[j]];
         kj++)
      ;
    if (kj < koff[ko + 1] && kcols[kidx[kj]] == cols[idx[j]])
      a_kj = kval[kj];
    else
      a_kj = 0;

    if (verbose) {
      printf("a_ij = %lf a_ik = %lf a_kj = %lf\n", val[j], a_ik, a_kj);
      fflush(stdout);
    }
    // a_ij = a_ij - a_ik * a_kj
    val[j] -= a_ik * a_kj;
  }
}

static void ilu0_level(int lvl, uint *lvl_off, struct par_mat *A,
                       struct par_mat *E, int verbose) {
  ulong *cols = A->cols, *rows = A->rows;
  uint *off = A->adj_off, *idx = A->adj_idx, i, k;
  for (i = lvl_off[lvl - 1] + (lvl == 1); i < lvl_off[lvl]; i++)
    for (k = off[i]; k < off[i + 1] && cols[idx[k]] < rows[i]; k++)
      ilu0_update_row(i, k, A, E, verbose, lvl);
}

static void ilu0(struct ilu *ilu, buffer *bfr) {
  ilu0_level(1, ilu->lvl_off, &ilu->A, NULL, 0);
  struct par_mat E;
  for (int l = 2; l <= ilu->nlvls; l++) {
    ilu0_get_rows(&E, l, ilu->lvl_off, &ilu->A, &ilu->cr, bfr);
    ilu0_level(l, ilu->lvl_off, &ilu->A, &E, 0);
    par_mat_free(&E);
  }
}

//=============================================================================
// ILUC
//
struct eij_t {
  ulong r, c;
  uint p;
  scalar v;
};

// We are going to separate A matrix to L and U where L is the strictly lower
// triangular part of A and U is the upper triangular part of A (including the
// diagonal). Since A is in CSR format, extracting U (in CSR format) is easy.
// L will be distributed by columns and we need to figure out the owner of a
// given column.
static void iluc_sep_lu(struct ilu *ilu, buffer *bfr) {
  // Recover the communicator
  struct crystal *cr = &ilu->cr;
  struct comm *c = &cr->comm;

  // Setup U
  struct par_mat *A = &ilu->A;
  struct array uijs, lijs;
  array_init(struct mij, &uijs, A->rn * 30);
  array_init(struct mij, &lijs, A->rn * 30);

  struct mij m = {.r = 0, .c = 0, .idx = 0, .p = 0, .v = 0};
  uint i, j, je;
  for (i = 0; i < A->rn; i++) {
    m.r = A->rows[i];
    j = A->adj_off[i], je = A->adj_off[i + 1];
    for (; j < je && A->cols[A->adj_idx[j]] < m.r; j++) {
      m.c = A->cols[A->adj_idx[j]], m.v = A->adj_val[j];
      m.p = m.c % c->np, m.idx = (local_dof(A->rows, m.c, A->rn) < A->rn);
      array_cat(struct mij, &lijs, &m, 1);
    }
    // Add the unit diagonal to L (We actually don't need to send this)
    m.c = m.r, m.v = 1, m.p = m.c % c->np, m.idx = 1;
    array_cat(struct mij, &lijs, &m, 1);

    for (; j < je; j++) {
      m.c = A->cols[A->adj_idx[j]], m.v = A->adj_val[j];
      array_cat(struct mij, &uijs, &m, 1);
    }
  }

  par_mat_setup(&ilu->U, &uijs, CSR, 0, bfr);
  array_free(&uijs);

  // Setup L
  sarray_transfer(struct mij, &lijs, p, 1, cr);
  if (lijs.n > 0) {
    sarray_sort_2(struct mij, lijs.ptr, lijs.n, c, 1, idx, 0, bfr);
    struct mij *pl = (struct mij *)lijs.ptr;
    for (i = 1, j = 0; i < lijs.n; i++) {
      if (pl[i].c != pl[j].c) {
        assert(pl[i - 1].idx == 1);
        for (; j < i; j++)
          pl[j].p = pl[i - 1].p;
        // j == i at the end
      }
    }
    // residual
    assert(pl[i - 1].idx == 1);
    for (; j < i; j++)
      pl[j].p = pl[i - 1].p;
  }

  sarray_transfer(struct mij, &lijs, p, 0, cr);
  par_mat_setup(&ilu->L, &lijs, CSC, 0, bfr);
  array_free(&lijs);
}

static void iluc_fwrd_rqsts(struct array *fwds, struct array *rqsts,
                            const int type, const ulong K,
                            const struct array *A, struct crystal *cr,
                            buffer *bfr) {
  fwds->n = rqsts->n = 0;
  struct request_t t = {.r = 0, .p = 0, .o = 1};

  struct comm *c = &cr->comm;

#define INIT_RQST(f, g, arr)                                                   \
  do {                                                                         \
    if (A->n > 0) {                                                            \
      sarray_sort_2(struct mij, A->ptr, A->n, f, 1, g, 1, bfr);                \
      struct mij *pa = (struct mij *)A->ptr;                                   \
      uint i = 1, j = 0;                                                       \
      for (; i < A->n; i++) {                                                  \
        if (pa[i].f != pa[j].f) {                                              \
          t.r = pa[j].f, t.p = t.r % c->np;                                    \
          array_cat(struct request_t, arr, &t, 1);                             \
          j = i;                                                               \
        }                                                                      \
      }                                                                        \
      if (j < i) {                                                             \
        t.r = pa[j].f, t.p = t.r % c->np;                                      \
        array_cat(struct request_t, arr, &t, 1);                               \
      }                                                                        \
    }                                                                          \
  } while (0)

  if (type == CSC)
    INIT_RQST(r, c, rqsts);
  else
    INIT_RQST(c, r, rqsts);
#undef INIT_RQST

  if (K > 0) {
    t.r = K, t.p = K % c->np, t.o = 0;
    array_cat(struct request_t, rqsts, &t, 1);
  }

  sarray_transfer(struct request_t, rqsts, p, 1, cr);

  // Okay, we got all the requests (if any) and non-zero row/col ids in the same
  // processor. Now we forward the requests to the original owners.
  if (rqsts->n > 0) {
    sarray_sort_2(struct request_t, rqsts->ptr, rqsts->n, r, 1, o, 0, bfr);
    struct request_t *pr = (struct request_t *)rqsts->ptr;
    uint s = 0, e = 1;
    for (; e < rqsts->n; e++) {
      if (pr[e].r != pr[s].r) {
        if (pr[s].o == 0) { // This is a request
          uint p = pr[s].p;
          for (s = s + 1; s < e; s++) {
            pr[s].o = p;
            array_cat(struct request_t, fwds, &pr[s], 1);
          }
        }
        s = e;
      }
    }
    if (s < e && pr[s].o == 0) {
      uint p = pr[s].p;
      for (s = s + 1; s < e; s++) {
        pr[s].o = p;
        array_cat(struct request_t, fwds, &pr[s], 1);
      }
    }
  }

  sarray_transfer(struct request_t, fwds, p, 0, cr);
}

static void iluc_send_data(struct array *data, const int type, struct array *A,
                           struct array *work, struct crystal *cr,
                           buffer *bfr) {
  if (type == CSC) {
    sarray_sort_2(struct mij, A->ptr, A->n, c, 1, r, 1, bfr);
    sarray_sort_2(struct eij_t, work->ptr, work->n, r, 1, c, 1, bfr);
  } else {
    sarray_sort_2(struct mij, A->ptr, A->n, r, 1, c, 1, bfr);
    sarray_sort_2(struct eij_t, work->ptr, work->n, c, 1, r, 1, bfr);
  }

  // We only have one request per processor, so sorting by processor is the
  // same as sorting by row id. But just to be safe we will sort by row id.
  data->n = 0;
  if (work->n > 0) {
    struct eij_t *pw = (struct eij_t *)work->ptr;
    uint i = 1, j = 0;
    for (; i < work->n; i++) {
      if ((pw[i].r != pw[j].r) || (pw[i].c != pw[j].c)) {
        array_cat(struct eij_t, data, &pw[j], 1);
        j = i;
      } else
        pw[j].v += pw[i].v;
    }
    if (j < i)
      array_cat(struct eij_t, data, &pw[j], 1);
  }

  sarray_transfer(struct eij_t, data, p, 0, cr);
}

static void iluc_get_data(struct array *data, ulong K, int type,
                          struct array *A, struct array *B, struct crystal *cr,
                          struct array *rqsts, struct array *fwds,
                          struct array *work, buffer *bfr) {
  iluc_fwrd_rqsts(fwds, rqsts, type, K, A, cr, bfr);

  work->n = 0;
  if (fwds->n > 0) {
    sarray_sort(struct request_t, fwds->ptr, fwds->n, r, 1, bfr);
    struct request_t *pf = (struct request_t *)fwds->ptr;

    uint i, j, k, l, n;
    scalar v;
    struct eij_t m = {.r = 0, .c = 0, .p = 0, .v = 0};

#define FILL_RQST(f, g, nd)                                                    \
  do {                                                                         \
    struct mij *pa = (struct mij *)A->ptr;                                     \
    struct mij *pb = (struct mij *)B->ptr;                                     \
    for (i = 0, j = 0; i < fwds->n; i++) {                                     \
      l = 0;                                                                   \
      m.f = pf[i].r, m.p = pf[i].o;                                            \
      for (; j < A->n && pa[j].f < m.f; j++)                                   \
        ;                                                                      \
      assert(j < A->n && pa[j].f == m.f);                                      \
      for (k = j; k < A->n && pa[k].f == m.f && pa[k].g < m.f; k++) {          \
        v = pa[k].v;                                                           \
        for (; l < B->n && pb[l].f < pa[k].g; l++)                             \
          ;                                                                    \
        assert(l < B->n && pb[l].f == pa[k].g);                                \
        for (n = l; n < B->n && pb[n].f == pa[k].g && (pb[n].g < m.f + nd);    \
             n++)                                                              \
          ;                                                                    \
        for (; n < B->n && pb[n].f == pa[k].g; n++) {                          \
          m.g = pb[n].g, m.v = -v * pb[n].v;                                   \
          array_cat(struct eij_t, work, &m, 1);                                \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  } while (0)

    if (type == CSC)
      FILL_RQST(r, c, 0);
    else
      FILL_RQST(c, r, 1);

#undef FILL_RQST
  }

  iluc_send_data(data, type, A, work, cr, bfr);
}

static void iluc_update(struct array *tij, ulong K, struct array *data, int row,
                        buffer *bfr) {
  // FIXME: This can be done more efficiently
  struct mij m = {.r = 0, .c = 0, .idx = 0, .p = 0, .v = 0};
  uint j;
  if (K) {
    if (row) {
      sarray_sort(struct eij_t, data->ptr, data->n, c, 1, bfr);
      struct eij_t *pd = (struct eij_t *)data->ptr;
      m.r = K;
      for (j = 0; j < data->n; j++) {
        m.c = pd[j].c, m.v = pd[j].v;
        array_cat(struct mij, tij, &m, 1);
      }
    } else {
      sarray_sort(struct eij_t, data->ptr, data->n, r, 1, bfr);
      struct eij_t *pd = (struct eij_t *)data->ptr;
      m.c = K;
      for (; j < data->n; j++) {
        m.r = pd[j].r, m.v = pd[j].v;
        array_cat(struct mij, tij, &m, 1);
      }
    }
  }

  struct array tmp;
  array_init(struct mij, &tmp, tij->n + 1);

  if (tij->n > 0) {
    uint i = 1, j = 0;
    struct mij *pt = NULL;
    if (row) {
      sarray_sort(struct mij, tij->ptr, tij->n, c, 1, bfr);
      pt = (struct mij *)tij->ptr;
      for (; i < tij->n; i++) {
        if (pt[i].c != pt[j].c) {
          array_cat(struct mij, &tmp, &pt[j], 1);
          j = i;
        } else
          pt[j].v += pt[i].v;
      }
    } else {
      sarray_sort(struct mij, tij->ptr, tij->n, r, 1, bfr);
      pt = (struct mij *)tij->ptr;
      for (; i < tij->n; i++) {
        if (pt[i].r != pt[j].r) {
          array_cat(struct mij, &tmp, &pt[j], 1);
          j = i;
        } else
          pt[j].v += pt[i].v;
      }
    }
    if (j < i && pt)
      array_cat(struct mij, &tmp, &pt[j], 1);

    tij->n = 0;
    array_cat(struct mij, tij, tmp.ptr, tmp.n);
  }

  array_free(&tmp);
}

static void iluc_level(struct array *lij, struct array *uij, int lvl,
                       struct ilu *ilu, struct array *data, struct array *work,
                       buffer *bfr) {
  // Work arrays
  struct array rij, cij;
  array_init(struct mij, &rij, 30);
  array_init(struct mij, &cij, 30);

  struct array rqst, fwds;
  array_init(struct request_t, &rqst, 30);
  array_init(struct request_t, &fwds, 30);

  struct par_mat *L = &ilu->L, *U = &ilu->U;
  struct crystal *cr = &ilu->cr;

  // Figure out start and end of the level and agree on a range
  uint *lvl_off = ilu->lvl_off, s = lvl_off[lvl - 1];
  sint buf[2], size = lvl_off[lvl] - s;
  comm_allreduce(&cr->comm, gs_int, gs_max, &size, 1, buf);
  uint e = s + size;

  uint i, j, je, k;
  for (k = s; k < e; k++) {
    ulong K = (k < lvl_off[lvl]) ? U->rows[k] : 0;

    // Fetch required data (combine with the other call below)
    iluc_get_data(data, K, CSC, lij, uij, cr, &rqst, &fwds, work, bfr);

    // Init z[1:K] = 0, z[K:n] = a_{K, K:n}, i.e., z = u_{K,:}
    rij.n = 0;
    if (K) {
      struct mij m = {.r = K, .c = 0, .idx = 0, .p = 0, .v = 0};
      for (j = U->adj_off[k], je = U->adj_off[k + 1]; j < je; j++) {
        m.c = U->cols[U->adj_idx[j]], m.v = U->adj_val[j];
        array_cat(struct mij, &rij, &m, 1);
      }
    }
    // Update z if l_KI != 0 for all I, 1 <= I < K
    iluc_update(&rij, K, data, 1, bfr);

    // Fetch required data (combine with the other call above)
    iluc_get_data(data, K, CSR, uij, lij, cr, &rqst, &fwds, work, bfr);

    // Init w[1:K] = 0, w[K] = 1, w[K+1:n] = a_{K+1:n, K}, i.e., w = l_{:, K}
    cij.n = 0;
    if (K) {
      struct mij m = {.r = 0, .c = K, .idx = 0, .p = 0, .v = 0};
      for (j = L->adj_off[k] + 1, je = L->adj_off[k + 1]; j < je; j++) {
        m.r = L->rows[L->adj_idx[j]], m.v = L->adj_val[j];
        array_cat(struct mij, &cij, &m, 1);
      }
    }
    // Update w if u_IK != 0 for all I, 1 <= I < K
    iluc_update(&cij, K, data, 0, bfr);

    // Set u_{k, :} = z and find u_kk
    scalar u_kk = 1;
    struct mij *pt = (struct mij *)rij.ptr;
    if (K) {
      if (rij.n > 0 && fabs(pt[0].v) > 1e-12)
        u_kk = pt[0].v;
      array_cat(struct mij, uij, rij.ptr, rij.n);
    }

    // Set l_{:, K} = w/u_KK and l_KK = 1
    pt = (struct mij *)cij.ptr;
    for (j = 0; j < cij.n; j++)
      pt[j].v /= u_kk;

    if (K) {
      struct mij m = {.r = K, .c = K, .idx = 0, .p = 0, .v = 1};
      array_cat(struct mij, &cij, &m, 1);
      array_cat(struct mij, lij, cij.ptr, cij.n);
    }
  }

  array_free(&rij), array_free(&cij);
  array_free(&rqst), array_free(&fwds);
}

//=============================================================================
// ILUCP
//
struct pivot_t {
  ulong k;
  uint p, pivot;
};

static void ilucp_get_data(struct array *data, ulong P, int type,
                           struct array *A, ulong K, struct array *B,
                           struct array *pvts, struct crystal *cr,
                           struct array *rqsts, struct array *fwds,
                           struct array *work, buffer *bfr) {
  iluc_fwrd_rqsts(fwds, rqsts, type, P, A, cr, bfr);

  work->n = 0;
  if (fwds->n > 0) {
    sarray_sort(struct request_t, fwds->ptr, fwds->n, r, 1, bfr);
    struct request_t *pf = (struct request_t *)fwds->ptr;

    uint i, j, k, l, n, o;
    scalar v;
    struct eij_t m = {.r = 0, .c = 0, .p = 0, .v = 0};

#define FILL_RQST(f, g, nd)                                                    \
  do {                                                                         \
    struct mij *pa = (struct mij *)A->ptr;                                     \
    struct mij *pb = (struct mij *)B->ptr;                                     \
    for (i = 0, j = 0; i < fwds->n; i++) {                                     \
      l = 0;                                                                   \
      m.f = pf[i].r, m.p = pf[i].o;                                            \
      for (; j < A->n && pa[j].f < m.f; j++)                                   \
        ;                                                                      \
      assert(j < A->n && pa[j].f == m.f);                                      \
      for (k = j; k < A->n && pa[k].f == m.f && pa[k].g < K; k++) {            \
        v = pa[k].v;                                                           \
        for (; l < B->n && pb[l].f < pa[k].g; l++)                             \
          ;                                                                    \
        assert(l < B->n && pb[l].f == pa[k].g);                                \
        for (n = l; n < B->n && pb[n].f == pa[k].g && (pb[n].g < K + nd); n++) \
          ;                                                                    \
        if (pvts != NULL) {                                                    \
          struct pivot_t *pp = (struct pivot_t *)pvts->ptr;                    \
          o = 0;                                                               \
          for (; n < B->n && pb[n].f == pa[k].g; n++) {                        \
            m.g = pb[n].g, m.v = -v * pb[n].v;                                 \
            while (o < pvts->n && pp[o].k < m.g)                               \
              o++;                                                             \
            assert(o < pvts->n && pp[o].k == m.g);                             \
            if (!pp[o].pivot)                                                  \
              array_cat(struct eij_t, work, &m, 1);                            \
          }                                                                    \
        } else {                                                               \
          for (; n < B->n && pb[n].f == pa[k].g; n++) {                        \
            m.g = pb[n].g, m.v = -v * pb[n].v;                                 \
            array_cat(struct eij_t, work, &m, 1);                              \
          }                                                                    \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  } while (0)

    if (type == CSC)
      FILL_RQST(r, c, 0);
    else
      FILL_RQST(c, r, 1);

#undef FILL_RQST
  }

  iluc_send_data(data, type, A, work, cr, bfr);
}

static ulong ilucp_find_pvt(ulong *perm, uint k, int lvl, uint *lvl_off,
                            struct array *row, struct crystal *cr,
                            buffer *bfr) {
  // First sort by the absolute value and then setup a gs handle to iteratively
  // select a pivot
  ulong p = 0;
  if (k < lvl_off[lvl]) {
    scalar v = 0;
    struct mij *pr = (struct mij *)row->ptr;
    for (uint i = 0; i < row->n && pr[i].c < lvl_off[lvl]; i++) {
      if (fabs(pr[i].v) > v) {
        v = fabs(pr[i].v);
        p = pr[i].c;
      }
    }
    perm[k] = p;
  }
  return p;
}

static void ilucp_update_pvts(struct array *pvts, struct array *rij,
                              ulong *perm, uint k, int lvl, uint *lvl_off,
                              struct crystal *cr, buffer *bfr) {
  struct comm *c = &cr->comm;

  struct pivot_t t = {.k = 0, .pivot = 0};
  struct mij *pr = (struct mij *)rij->ptr;
  for (uint i = 0; i < rij->n; i++) {
    t.k = pr[i].c, t.p = t.k % c->np;
    array_cat(struct pivot_t, pvts, &t, 1);
  }

  uint e = (k < lvl_off[lvl] ? k : lvl_off[lvl]);
  t.pivot = 1;
  for (uint i = 0; i < e; i++) {
    t.k = perm[i], t.p = t.k % c->np;
    array_cat(struct pivot_t, pvts, &t, 1);
  }

  if (pvts->n > 0) {
    struct array temp;
    array_init(struct pivot_t, &temp, pvts->n + 1);

    sarray_sort_2(struct pivot_t, pvts->ptr, pvts->n, k, 1, pivot, 1, bfr);
    struct pivot_t *pp = (struct pivot_t *)pvts->ptr;
    uint i = 1, j = 0;
    for (; i < pvts->n; i++) {
      if (pp[i].k != pp[j].k) {
        array_cat(struct pivot_t, &temp, &pp[i - 1], 1);
        j = i;
      }
    }
    if (j < i)
      array_cat(struct pivot_t, &temp, &pp[i - 1], 1);
    pvts->n = 0;
    array_cat(struct pivot_t, pvts, temp.ptr, temp.n);
    array_free(&temp);
  }

  sarray_transfer(struct pivot_t, pvts, p, 1, cr);
  sarray_sort_2(struct pivot_t, pvts->ptr, pvts->n, k, 1, pivot, 0, bfr);

  if (pvts->n > 0) {
    struct pivot_t *pp = (struct pivot_t *)pvts->ptr;
    uint i = 1, j = 0;
    for (; i < pvts->n; i++) {
      if (pp[i].k != pp[j].k) {
        for (; j < i - 1; j++)
          pp[j].pivot = pp[i - 1].pivot;
        j = i;
      }
    }
    if (j < i) {
      for (; j < i - 1; j++)
        pp[j].pivot = pp[i - 1].pivot;
    }
  }

  sarray_transfer(struct pivot_t, pvts, p, 1, cr);
  sarray_sort(struct pivot_t, pvts->ptr, pvts->n, k, 1, bfr);
}

static void ilucp_level(struct array *lij, struct array *uij, int lvl,
                        struct ilu *ilu, struct array *pvts, struct array *data,
                        struct array *work, buffer *bfr) {
  // Work arrays
  struct array rij, cij;
  array_init(struct mij, &rij, 30);
  array_init(struct mij, &cij, 30);

  struct array rqst, fwds;
  array_init(struct request_t, &rqst, 30);
  array_init(struct request_t, &fwds, 30);

  struct par_mat *L = &ilu->L, *U = &ilu->U;
  struct crystal *cr = &ilu->cr;

  // Figure out start and end of the level and agree on a range
  uint *lvl_off = ilu->lvl_off, s = lvl_off[lvl - 1];
  sint buf[2], size = lvl_off[lvl] - s;
  comm_allreduce(&cr->comm, gs_int, gs_max, &size, 1, buf);
  uint e = s + size;

  uint i, j, je, k, l;
  for (k = s; k < e; k++) {
    ulong K = (k < lvl_off[lvl]) ? U->rows[k] : 0;

    // Fetch required data. We will skip the data in the  columns which were
    // choosen as pivots.
    ilucp_get_data(data, K, CSC, lij, K, uij, pvts, cr, &rqst, &fwds, work,
                   bfr);

    // Init z[1:K] = 0, z[K:n] = a_{K, K:n}, i.e., z = u_{K,:} and skip the
    // columns which have been choosen as pivots.
    rij.n = 0;
    if (K) {
      struct mij m = {.r = K, .c = 0, .idx = 0, .p = 0, .v = 0};
      struct pivot_t *pp = (struct pivot_t *)pvts->ptr;
      for (j = U->adj_off[k], je = U->adj_off[k + 1], l = 0; j < je; j++) {
        m.c = U->cols[U->adj_idx[j]], m.v = U->adj_val[j];
        while (l < pvts->n && pp[l].k < m.c)
          l++;
        assert(pp[l].k == m.c);
        if (!pp[l].pivot)
          array_cat(struct mij, &rij, &m, 1);
      }
    }

    // Update z if l_KI != 0 for all I, 1 <= I < K
    iluc_update(&rij, K, data, 1, bfr);

    // Select the pivot now -- all the active processors have to agree on their
    // own pivot. If two processors share the same pivot, smallest one wins and
    // others have to concede and find another one. So we will send a pivot
    // candidate list and make each processor pick one. Right now the candidate
    // list = updated row.
    ulong P = ilucp_find_pvt(ilu->perm, k, lvl, lvl_off, &rij, cr, bfr);

    // Sync the pivots: Basically everyone gets updated about which cols of U
    // have become pivots. Can't be done through a gs call, will have to send
    // all the cols in U and the current row a_k along with the info if its
    // a pivot.
    ilucp_update_pvts(pvts, &rij, ilu->perm, k, lvl, lvl_off, cr, bfr);

    // Fetch required data for col updated. Can't combine with above call when
    // we pivot? Will need to reimplement this part
    ilucp_get_data(data, P, CSR, uij, K, lij, NULL, cr, &rqst, &fwds, work,
                   bfr);

    // Init w[1:K] = 0, w[K] = 1, w[K+1:n] = a_{K+1:n, K}, i.e., w = l_{:, K}
    cij.n = 0;
    if (K) {
      struct mij m = {.r = 0, .c = K, .idx = 0, .p = 0, .v = 0};
      for (j = L->adj_off[k] + 1, je = L->adj_off[k + 1]; j < je; j++) {
        m.r = L->rows[L->adj_idx[j]], m.v = L->adj_val[j];
        array_cat(struct mij, &cij, &m, 1);
      }
    }
    // Update w if u_IK != 0 for all I, 1 <= I < K
    iluc_update(&cij, K, data, 0, bfr);

    // Set u_{k, :} = z and find u_kk
    // FIXME: This should u_{perm[k],perm[k]}}, not u_kk
    scalar u_kk = 1;
    struct mij *pt = (struct mij *)rij.ptr;
    if (K) {
      if (rij.n > 0 && fabs(pt[0].v) > 1e-12)
        u_kk = pt[0].v;
      array_cat(struct mij, uij, rij.ptr, rij.n);
    }

    // Set l_{:, K} = w/u_KK and l_KK = 1
    pt = (struct mij *)cij.ptr;
    for (j = 0; j < cij.n; j++)
      pt[j].v /= u_kk;

    if (K) {
      struct mij m = {.r = K, .c = K, .idx = 0, .p = 0, .v = 1};
      array_cat(struct mij, &cij, &m, 1);
      array_cat(struct mij, lij, cij.ptr, cij.n);
    }
  }

  array_free(&rij), array_free(&cij), array_free(&rqst), array_free(&fwds);
}

static void iluc(struct ilu *ilu, buffer *bfr) {
  struct crystal *cr = &ilu->cr;
  struct comm *c = &cr->comm;

  // Setup L and U
  iluc_sep_lu(ilu, bfr);

  struct par_mat *A = &ilu->A, *L = &ilu->L, *U = &ilu->U;

  struct array uij, lij, data, work;
  array_init(struct mij, &uij, A->rn * 30 + 1);
  array_init(struct mij, &lij, A->rn * 30 + 1);
  array_init(struct eij_t, &data, A->rn * 30 + 1);
  array_init(struct eij_t, &work, A->rn * 30 + 1);

  struct array pvts;
  array_init(struct pivot_t, &pvts, L->cn + 1);

  if (ilu->pivot) {
    ilu->perm = tcalloc(ulong, A->rn);
    // Initialize with the columns of U, i.e, columns of L
    struct pivot_t t = {.k = 0, .p = 0, .pivot = 0};
    for (uint i = 0; i < U->cn; i++) {
      t.k = U->cols[i], t.p = t.k % c->np;
      array_cat(struct pivot_t, &pvts, &t, 1);
    }

    for (int l = 1; l <= ilu->nlvls; l++)
      ilucp_level(&lij, &uij, l, ilu, &pvts, &data, &work, bfr);
  } else {
    for (int l = 1; l <= ilu->nlvls; l++)
      iluc_level(&lij, &uij, l, ilu, &data, &work, bfr);
  }

  par_mat_free(L), par_mat_free(U);
  par_mat_setup(U, &uij, CSR, 0, bfr);
  par_mat_setup(L, &lij, CSC, 0, bfr);

  const char *val = getenv("PARRSB_DUMP_ILU");
  if (val != NULL && atoi(val) != 0) {
    par_mat_dump("LL.txt", L, cr, bfr);
    par_mat_dump("UU.txt", U, cr, bfr);
  }

  array_free(&pvts);
  array_free(&lij), array_free(&uij);
  array_free(&work), array_free(&data);
}

//=============================================================================
// ILU API related functions
//
// `vtx` array is in the order of sorted element ids
static int ilu_setup_aux(struct ilu *ilu, int nlvls, uint *lvl_off,
                         uint *lvl_owner, ulong *lvl_ids, const uint n,
                         const int nv, const slong *vtx, const int verbose,
                         buffer *bfr) {
  struct elm {
    slong vtx[8];
    uint p, lvl;
    ulong e;
  };

  struct crystal *cr = &ilu->cr;
  struct comm *c = &cr->comm;

  // Send the elements in each level to the owner
  struct array elms;
  array_init(struct elm, &elms, n);

  struct elm elm;
  for (int l = 0; l < nlvls; l++) {
    for (uint i = lvl_off[l]; i < lvl_off[l + 1]; i++) {
      elm.lvl = l + 1, elm.e = lvl_ids[i], elm.p = lvl_owner[i];
      array_cat(struct elm, &elms, &elm, 1);
    }
  }
  sarray_sort(struct elm, elms.ptr, elms.n, e, 1, bfr);

  struct elm *pe = (struct elm *)elms.ptr;
  if (elms.n > 0) {
    // Sanity check
    assert(elms.n == n);
    for (uint i = 0; i < n; i++) {
      for (int v = 0; v < nv; v++)
        pe[i].vtx[v] = vtx[i * nv + v];
    }
  }

  sarray_transfer(struct elm, &elms, p, 1, cr);
  sarray_sort_2(struct elm, elms.ptr, elms.n, lvl, 0, e, 1, bfr);

  // Setup the ILU structure: allocate ILU data structures.
  ilu->nlvls = nlvls;
  ilu->lvl_off = (uint *)tcalloc(uint, ilu->nlvls + 1);

  uint s = 0, e = 0;
  ilu->lvl_off[0] = s;
  pe = (struct elm *)elms.ptr;
  for (int l = 1; l <= ilu->nlvls; l++) {
    while (e < elms.n && pe[e].lvl == l)
      e++;
    ilu->lvl_off[l] = ilu->lvl_off[l - 1] + e - s;
    s = e;
  }

  // Number rows now: All the elements in Level 0 are numbered before Level
  // 1 and so on.
  ulong *ids = trealloc(ulong, ids, elms.n);
  ulong ng = 0;
  for (int l = 0; l < ilu->nlvls; l++) {
    e = ilu->lvl_off[l + 1], s = ilu->lvl_off[l];
    slong out[2][1], buf[2][1], in = e - s;
    comm_scan(out, c, gs_long, gs_add, &in, 1, buf);
    ulong start = ng + out[0][0] + 1;
    for (; s < e; s++)
      ids[s] = start++;
    ng += out[1][0];
  }

  slong *vrt = tcalloc(slong, elms.n * nv);
  for (uint i = 0; i < elms.n; i++) {
    for (int j = 0; j < nv; j++)
      vrt[i * nv + j] = pe[i].vtx[j];
  }

  if (verbose > 1) {
    for (uint i = 0; i < elms.n; i++) {
      printf("fid = %llu, ", ids[i]);
      for (int v = 0; v < nv; v++)
        printf("%lld, ", vrt[i * nv + v]);
      printf("\n");
      fflush(stdout);
    }
  }

  // Find and compress neighbors in order to form the Laplacian
  struct array nbrs, eij;
  find_nbrs(&nbrs, ids, vrt, elms.n, nv, cr, bfr);
  compress_nbrs(&eij, &nbrs, bfr);
  free(ids), free(vrt);
  array_free(&elms), array_free(&nbrs);

  // Setup the parallel CSR matrix
  par_csr_setup(&ilu->A, &eij, 0, bfr);
  array_free(&eij);

  return 0;
}

struct ilu *ilu_setup(const uint n, const int nv, const long long *llvtx,
                      const ilu_options *options, MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  struct ilu *ilu = tcalloc(struct ilu, 1);
  ilu->pivot = options->pivot, ilu->verbose = options->verbose;
  ilu->tol = options->tol, ilu->nnz_per_row = options->nnz_per_row;
  ilu->lvl_off = NULL, ilu->perm = NULL;
  crystal_init(&ilu->cr, &c);

  slong *vtx = tcalloc(slong, n * nv);
  for (uint i = 0; i < n * nv; i++)
    vtx[i] = llvtx[i];

  // Establish a numbering based on input
  slong out[2][1], buf[2][1], in = n;
  comm_scan(out, &c, gs_long, gs_add, &in, 1, buf);
  ulong s = out[0][0], ng = out[1][0];

  ulong *ids = tcalloc(ulong, n);
  for (uint i = 0; i < n; i++)
    ids[i] = s + i + 1;

  buffer bfr;
  buffer_init(&bfr, 1024);

  uint *lvl_off = tcalloc(uint, 100 + n), *lvl_owner = lvl_off + 100;
  ulong *lvl_ids = tcalloc(ulong, n);
  int nlvls = find_lvls(lvl_off, lvl_owner, lvl_ids, n, nv, ids, vtx, 1,
                        &ilu->cr, ilu->verbose, &bfr);
  ilu_setup_aux(ilu, nlvls, lvl_off, lvl_owner, lvl_ids, n, nv, vtx,
                ilu->verbose, &bfr);

  char *val = getenv("PARRSB_DUMP_ILU");
  if (val != NULL && atoi(val) != 0)
    par_mat_dump("A.txt", &ilu->A, &ilu->cr, &bfr);

  // Setup the ILU factors
  switch (options->type) {
  case 0:
    ilu0(ilu, &bfr);
    break;
  case 1:
    iluc(ilu, &bfr);
    break;
  default:
    break;
  }

  val = getenv("PARRSB_DUMP_ILU");
  if (val != NULL && atoi(val) != 0)
    par_mat_dump("B.txt", &ilu->A, &ilu->cr, &bfr);

  free(ids), free(vtx), free(lvl_off), free(lvl_ids);
  buffer_free(&bfr), comm_free(&c);

  return ilu;
}

void ilu_free(struct ilu *ilu) {
  if (ilu) {
    crystal_free(&ilu->cr);
    if (ilu->nlvls > 0) {
      par_mat_free(&ilu->A);
      // FIXME: Cleanup L and U
      // par_mat_free(&ilu->L);
      // par_mat_free(&ilu->U);
    }
    if (ilu->lvl_off)
      free(ilu->lvl_off), ilu->lvl_off = NULL;
    if (ilu->perm)
      free(ilu->perm), ilu->perm = NULL;
    free(ilu);
  }
}

#undef CSC
#undef CSR

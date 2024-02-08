#ifdef SUFFIX
#define SUFFIXED_NAME(x) TOKEN_PASTE(x, SUFFIX)
#else
#define SUFFIXED_NAME(x) x
#endif

#if defined(tensor_dot)
#undef tensor_dot
#endif
#define tensor_dot SUFFIXED_NAME(tensor_dot)
T tensor_dot(const T *a, const T *b, uint n) {
  T sum = 0;
  for (; n; --n)
    sum += *a++ * *b++;
  return sum;
}

/* factors: L is in CSR format
            D is a diagonal matrix stored as a vector
   actual factorization is:

                  -1      T
     A   = (I-L) D   (I-L)

      -1        -T        -1
     A   = (I-L)   D (I-L)

   (triangular factor is unit diagonal; the diagonal is not stored)
*/
#define sparse_cholesky SUFFIXED_NAME(sparse_cholesky)
struct sparse_cholesky {
  uint n, *Lrp, *Lj;
  T *L, *D;
};

/*
  symbolic factorization: finds the sparsity structure of L

  uses the concept of elimination tree:
    the parent of node j is node i when L(i,j) is the first
      non-zero in column j below the diagonal (i>j)
    L's structure is discovered row-by-row; the first time
      an entry in column j is set, it must be the parent

  the nonzeros in L are the nonzeros in A + paths up the elimination tree

  linear in the number of nonzeros of L
*/
#define factor_symbolic SUFFIXED_NAME(factor_symbolic)
static void factor_symbolic(uint n, const uint *Arp, const uint *Aj,
                            struct sparse_cholesky *out, buffer *buf) {
  uint *visit = tmalloc(uint, 2 * n), *parent = visit + n;
  uint *Lrp, *Lj;
  uint i, nz = 0;

  out->n = n;

  for (i = 0; i < n; ++i) {
    uint p = Arp[i], pe = Arp[i + 1];
    visit[i] = i, parent[i] = n;
    for (; p != pe; ++p) {
      uint j = Aj[p];
      if (j >= i)
        break;
      for (; visit[j] != i; j = parent[j]) {
        ++nz, visit[j] = i;
        if (parent[j] == n) {
          parent[j] = i;
          break;
        }
      }
    }
  }

  Lrp = out->Lrp = tmalloc(uint, n + 1 + nz);
  Lj = out->Lj = Lrp + n + 1;

  Lrp[0] = 0;
  for (i = 0; i < n; ++i) {
    uint p = Arp[i], pe = Arp[i + 1], count = 0, *Ljr = &Lj[Lrp[i]];
    visit[i] = i;
    for (; p != pe; ++p) {
      uint j = Aj[p];
      if (j >= i)
        break;
      for (; visit[j] != i; j = parent[j])
        Ljr[count++] = j, visit[j] = i;
    }
    sortv(Ljr, Ljr, count, sizeof(uint), buf);
    Lrp[i + 1] = Lrp[i] + count;
  }
  free(visit);
}

/*
  numeric factorization:

  L is built row-by-row, using:    ( ' indicates transpose )


  [ A  r ]  = [ (I-L)   ] [ D^(-1)  ] [ (I-L)' -s ]
  [ r' a ]    [  -s'  1 ] [     1/d ] [         1 ]

            = [ A   (I-L) D^(-1) (-s)  ]
              [ r'  s' D^(-1) s + 1/d  ]

  so, if r' is the next row of A, up to but excluding the diagonal,
  then the next row of L, s', obeys

     r = - (I-L) D^(-1) s

  let y = (I-L)^(-1) (-r)
  then s = D y, and d = 1/(s' y)

*/
#define factor_numeric SUFFIXED_NAME(factor_numeric)
static void factor_numeric(uint n, const uint *Arp, const uint *Aj, const T *A,
                           struct sparse_cholesky *out, uint *visit, T *y) {
  const uint *Lrp = out->Lrp, *Lj = out->Lj;
  T *D, *L;
  uint i;

  out->D = D = tmalloc(T, n + Lrp[n]);
  out->L = L = D + n;

  for (i = 0; i < n; ++i) {
    uint p, pe;
    T a = 0;
    visit[i] = n;
    for (p = Lrp[i], pe = Lrp[i + 1]; p != pe; ++p) {
      uint j = Lj[p];
      y[j] = 0, visit[j] = i;
    }
    for (p = Arp[i], pe = Arp[i + 1]; p != pe; ++p) {
      uint j = Aj[p];
      if (j >= i) {
        if (j == i)
          a = A[p];
        break;
      }
      y[j] = -A[p];
    }
    for (p = Lrp[i], pe = Lrp[i + 1]; p != pe; ++p) {
      uint q, qe, j = Lj[p];
      T lij, yj = y[j];
      for (q = Lrp[j], qe = Lrp[j + 1]; q != qe; ++q) {
        uint k = Lj[q];
        if (visit[k] == i)
          yj += L[q] * y[k];
      }
      y[j] = yj;
      L[p] = lij = D[j] * yj;
      a -= yj * lij;
    }
    D[i] = 1 / a;
  }
}

/* x = A^(-1) b;  works when x and b alias */
#define sparse_cholesky_solve SUFFIXED_NAME(sparse_cholesky_solve)
static void sparse_cholesky_solve(T *x, const struct sparse_cholesky *fac,
                                  T *b) {
  const uint n = fac->n, *Lrp = fac->Lrp, *Lj = fac->Lj;
  uint i, p, pe;
  const T *L = fac->L, *D = fac->D;
  for (i = 0; i < n; ++i) {
    T xi = b[i];
    for (p = Lrp[i], pe = Lrp[i + 1]; p != pe; ++p)
      xi += L[p] * x[Lj[p]];
    x[i] = xi;
  }
  for (i = 0; i < n; ++i)
    x[i] *= D[i];
  for (i = n; i;) {
    T xi = x[--i];
    for (p = Lrp[i], pe = Lrp[i + 1]; p != pe; ++p)
      x[Lj[p]] += L[p] * xi;
  }
}

#define sparse_cholesky_factor SUFFIXED_NAME(sparse_cholesky_factor)
static void sparse_cholesky_factor(uint n, const uint *Arp, const uint *Aj,
                                   const T *A, struct sparse_cholesky *out,
                                   buffer *buf) {
  const uint n_uints_as_T = (n * sizeof(uint) + sizeof(T) - 1) / sizeof(T);
  buffer_reserve(buf, (n_uints_as_T + n) * sizeof(T));
  factor_symbolic(n, Arp, Aj, out, buf);
  factor_numeric(n, Arp, Aj, A, out, (uint *)buf->ptr,
                 n_uints_as_T + (T *)buf->ptr);
}

#define sparse_cholesky_free SUFFIXED_NAME(sparse_cholesky_free)
static void sparse_cholesky_free(struct sparse_cholesky *fac) {
  free(fac->Lrp);
  fac->Lj = fac->Lrp = 0;
  free(fac->D);
  fac->L = fac->D = 0;
}

#define csr_mat SUFFIXED_NAME(csr_mat)
struct csr_mat {
  uint n, *Arp, *Aj;
  T *A;
};

#define xxt SUFFIXED_NAME(xxt)
struct xxt {
  /* communication */

  struct comm comm;
  uint pcoord;      /* coordinate in communication tree */
  unsigned plevels; /* # of stages of communication */
  sint *pother;     /* let p = pother[i], then during stage i of fan-in,
                           if p>=0, receive from p
                           if p< 0, send to (-p-1)
                       fan-out is just the reverse ...
                       on proc 0, pother is never negative
                       on others, pother is negative for the last stage only */
  comm_req *req;

  /* separators */

  unsigned nsep;  /* number of separators */
  uint *sep_size; /* # of dofs on each separator,
                     ordered from the bottom to the top of the tree:
                     separator 0 is the bottom-most one (dofs not shared)
                     separator nsep-1 is the root of the tree */

  unsigned null_space;
  T *share_weight;

  /* vector sizes */

  uint un; /* user's vector size */

  /* xxt_solve works with "condensed" vectors;
     same dofs as in user's vectors, but no duplicates and no Dirichlet nodes,
     and also ordered topologically (children before parents) according to the
     separator tree */

  uint cn;        /* size of condensed vectors */
  sint *perm_u2c; /* permutation from user vector to condensed vector,
                     p=perm_u2c[i]; xu[i] = p=-1 ? 0 : xc[p];          */
  uint ln, sn;    /* xc[0 ... ln-1] are not shared   (ln=sep_size[0])
                     xc[ln ... ln+sn-1] are shared
                     ln+sn = cn                    */

  uint xn; /* # of columns of x = sum_i(sep_size[i]) - sep_size[0] */

  /* data */
  struct sparse_cholesky fac_A_ll;
  struct csr_mat A_sl;
  uint *Xp;
  T *X; /* column i of X starts at X[Xp[i]] */

  /* execution buffers */
  T *vl, *vc, *vx, *combuf;
};

/*
  for the binary communication tree, the procs are divided in half
  at each level, with the second half always the larger

  e.g., for np = 13:

       +------13-------+
       |               |
   +---6---+       +---7---+
   |       |       |       |
 +-3-+   +-3-+   +-3-+   +-4-+
 1   2   1   2   1   2   2   2
    1^1     1^1     1^1 1^1 1^1

  plevels is the number of levels in the tree
    = np==1 ? 1 : ( lg(np-1)+2 )

  labelling the nodes with proc id's gives this communication tree:

       +-------0-------+
       |               |
   +---0---+       +---6---+
   |       |       |       |
 +-0-+   +-3-+   +-6-+   +-9-+
 0   1   3   4   6   7   9   b
    1^2     4^5     7^8 9^a b^c

  consider proc 7 (pid = 7);
  pcoord gives the position of the leaf labelled 7:
    Root Right Left Right Left -> RRLRL -> 11010
    so pcoord = 11010 binary
  note the parent coordinate can be found by bit shifting right
    (i.e. dividing by 2)
*/

/* sets: pcoord, nsep, plevels, pother, req */
#define locate_proc SUFFIXED_NAME(locate_proc)
static void locate_proc(struct xxt *data) {
  const uint id = data->comm.id;
  uint n = data->comm.np, c = 1, odd = 0, base = 0;
  unsigned level = 0;
  while (n > 1) {
    ++level;
    odd = (odd << 1) | (n & 1);
    c <<= 1, n >>= 1;
    if (id >= base + n)
      c |= 1, base += n, n += (odd & 1);
  }
  data->pcoord = c;
  data->nsep = level + 1;
  data->plevels = data->nsep - 1;
  data->pother = tmalloc(sint, data->plevels);
  data->req = tmalloc(comm_req, data->plevels);
  for (level = 0; level < data->plevels; ++level) {
    if ((c & 1) == 1) {
      uint targ = id - (n - (odd & 1));
      data->pother[level] = -(sint)(targ + 1);
      data->plevels = level + 1;
      break;
    } else {
      data->pother[level] = id + n;
      c >>= 1, n = (n << 1) + (odd & 1), odd >>= 1;
    }
  }
}

/* determine the size of each separator;
   sums the separator sizes following the fan-in, fan-out comm. pattern
   uses the share-counts to avoid counting dofs more than once */
/* sets: xn, sep_size, ln, sn */
#define discover_sep_sizes SUFFIXED_NAME(discover_sep_sizes)
static void discover_sep_sizes(struct xxt *data, struct array *dofa,
                               buffer *buf) {
  const unsigned ns = data->nsep, nl = data->plevels;
  const uint n = dofa->n;
  float *v, *recv;
  unsigned i, lvl;
  uint j;
  const struct dof *dof = (const struct dof *)dofa->ptr;

  buffer_reserve(buf, 2 * ns * sizeof(float));
  v = (float *)buf->ptr, recv = v + ns;

  for (i = 0; i < ns; ++i)
    v[i] = 0;
  for (j = 0; j < n; ++j)
    v[dof[j].level] += 1 / (float)dof[j].count;

  /* fan-in */
  for (lvl = 0; lvl < nl; ++lvl) {
    sint other = data->pother[lvl];
    unsigned s = ns - (lvl + 1);
    if (other < 0) {
      comm_send(&data->comm, v + lvl + 1, s * sizeof(float), -other - 1, s);
    } else {
      comm_recv(&data->comm, recv + lvl + 1, s * sizeof(float), other, s);
      for (i = lvl + 1; i < ns; ++i)
        v[i] += recv[i];
    }
  }
  /* fan-out */
  for (; lvl;) {
    sint other = data->pother[--lvl];
    unsigned s = ns - (lvl + 1);
    if (other < 0)
      comm_recv(&data->comm, v + lvl + 1, s * sizeof(float), -other - 1, s);
    else
      comm_send(&data->comm, v + lvl + 1, s * sizeof(float), other, s);
  }

  data->xn = 0;
  data->sep_size = tmalloc(uint, ns);
  for (i = 0; i < ns; ++i) {
    uint s = v[i] + .1f;
    data->sep_size[i] = s;
    data->xn += s;
  }
  data->ln = data->sep_size[0];
  data->sn = data->cn - data->ln;
  data->xn -= data->ln;
}

/* assuming [A,Aend) is sorted,
   removes 0's and any duplicate entries,
   returns new end */
#define unique_nonzero SUFFIXED_NAME(unique_nonzero)
static ulong *unique_nonzero(ulong *A, ulong *Aend) {
  if (Aend == A)
    return A;
  else {
    ulong *end = Aend - 1, last = *end, *p = A, *q = A, v = 0;
    *end = 1;
    while (*q == 0)
      ++q; /*  *q==0 => q!=end since *end==0 */
    *end = 0;
    while (q != end) {
      v = *q++, *p++ = v;
      while (*q == v)
        ++q; /*  *q==v => q!=end since *end==0 */
    }
    if (last != v)
      *p++ = last;
    return p;
  }
}

#define merge_sep_ids SUFFIXED_NAME(merge_sep_ids)
static void merge_sep_ids(struct xxt *data, ulong *sep_id, ulong *other,
                          ulong *work, unsigned s0, buffer *buf) {
  const unsigned ns = data->nsep;
  unsigned s;
  ulong *p = sep_id, *q = other;
  for (s = s0; s < ns; ++s) {
    ulong *end;
    uint size = data->sep_size[s];
    memcpy(work, p, size * sizeof(ulong));
    memcpy(work + size, q, size * sizeof(ulong));
    sortv_long(work, work, 2 * size, sizeof(ulong), buf);
    end = unique_nonzero(work, work + 2 * size);
    memcpy(p, work, (end - work) * sizeof(ulong));
    p += size, q += size;
  }
}

#define init_sep_ids SUFFIXED_NAME(init_sep_ids)
static void init_sep_ids(struct xxt *data, struct array *dofa, ulong *xid) {
  const unsigned ns = data->nsep;
  const uint n = data->cn, *sep_size = data->sep_size;
  unsigned s = 1;
  uint i, size;
  const struct dof *dof = (const struct dof *)dofa->ptr;
  if (ns == 1)
    return;
  size = sep_size[s];
  for (i = data->ln; i < n; ++i) {
    unsigned si = dof[i].level;
    while (s != si) {
      memset(xid, 0, size * sizeof(ulong));
      xid += size;
      if (++s != ns)
        size = data->sep_size[s];
    }
    *xid++ = dof[i].id, --size;
  }
  while (s != ns) {
    memset(xid, 0, size * sizeof(ulong));
    xid += size;
    if (++s != ns)
      size = data->sep_size[s];
  }
}

#define find_perm_x2c SUFFIXED_NAME(find_perm_x2c)
static void find_perm_x2c(uint ln, uint cn, const struct array *dofc, uint xn,
                          const ulong *xid, sint *perm) {
  const struct dof *dof = (const struct dof *)dofc->ptr, *dof_end = dof + cn;
  const ulong *xid_end = xid + xn;
  uint i = ln;
  dof += ln;
  while (dof != dof_end) {
    ulong v = dof->id;
    while (*xid != v)
      ++xid, *perm++ = -1;
    *perm++ = i++, ++dof, ++xid;
  }
  while (xid != xid_end)
    ++xid, *perm++ = -1;
}

/* sets: perm_x2c */
#define discover_sep_ids SUFFIXED_NAME(discover_sep_ids)
static sint *discover_sep_ids(struct xxt *data, struct array *dofa,
                              buffer *buf) {
  const unsigned ns = data->nsep, nl = data->plevels;
  const uint xn = data->xn, *sep_size = data->sep_size;
  ulong *xid, *recv, *work, *p;
  unsigned lvl;
  uint size, ss;
  sint *perm_x2c;

  size = 0;
  for (lvl = 1; lvl < ns; ++lvl)
    if (sep_size[lvl] > size)
      size = sep_size[lvl];
  xid = tmalloc(ulong, 2 * xn + 2 * size), recv = xid + xn, work = recv + xn;

  init_sep_ids(data, dofa, xid);

  if (nl) {
    /* fan-in */
    p = xid, size = xn;
    for (lvl = 0; lvl < nl; ++lvl) {
      sint other = data->pother[lvl];
      if (other < 0) {
        comm_send(&data->comm, p, size * sizeof(ulong), -other - 1, size);
      } else {
        comm_recv(&data->comm, recv, size * sizeof(ulong), other, size);
        merge_sep_ids(data, p, recv, work, lvl + 1, buf);
      }
      ss = data->sep_size[lvl + 1];
      if (ss >= size || lvl == nl - 1)
        break;
      p += ss, size -= ss;
    }
    /* fan-out */
    for (;;) {
      sint other = data->pother[lvl];
      if (other < 0)
        comm_recv(&data->comm, p, size * sizeof(ulong), -other - 1, size);
      else
        comm_send(&data->comm, p, size * sizeof(ulong), other, size);
      if (lvl == 0)
        break;
      ss = data->sep_size[lvl];
      p -= ss, size += ss, --lvl;
    }
  }

  perm_x2c = tmalloc(sint, xn);
  find_perm_x2c(data->ln, data->cn, dofa, xn, xid, perm_x2c);
  free(xid);

  return perm_x2c;
}

#define apply_QQt SUFFIXED_NAME(apply_QQt)
static void apply_QQt(struct xxt *data, T *v, uint n, uint tag) {
  const unsigned nl = data->plevels;
  T *p = v, *recv = data->combuf;
  unsigned lvl, nsend = 0;
  uint size = n, ss;

  if (n == 0 || nl == 0)
    return;

  tag = tag * 2 + 0;
  /* fan-in */
  for (lvl = 0; lvl < nl; ++lvl) {
    sint other = data->pother[lvl];
    if (other < 0) {
      comm_send(&data->comm, p, size * sizeof(T), -other - 1, tag);
    } else {
      uint i;
      comm_recv(&data->comm, recv, size * sizeof(T), other, tag);
      for (i = 0; i < size; ++i)
        p[i] += recv[i];
    }
    ss = data->sep_size[lvl + 1];
    if (ss >= size || lvl == nl - 1)
      break;
    p += ss, size -= ss;
  }
  /* fan-out */
  for (;;) {
    sint other = data->pother[lvl];
    if (other < 0) {
      comm_recv(&data->comm, p, size * sizeof(T), -other - 1, tag);
    } else {
      comm_isend(&data->req[nsend++], &data->comm, p, size * sizeof(T), other,
                 tag);
    }
    if (lvl == 0)
      break;
    ss = data->sep_size[lvl];
    p -= ss, size += ss, --lvl;
  }
  if (nsend)
    comm_wait(data->req, nsend);
}

#define sum SUFFIXED_NAME(sum)
static double sum(struct xxt *data, double v, uint n, uint tag) {
  const unsigned nl = data->plevels;
  double r;
  unsigned lvl, nsend = 0;
  uint size = n, ss;

  tag = tag * 2 + 1;
  if (n == 0 || nl == 0)
    return v;
  /* fan-in */
  for (lvl = 0; lvl < nl; ++lvl) {
    sint other = data->pother[lvl];
    if (other < 0) {
      comm_send(&data->comm, &v, sizeof(double), -other - 1, tag);
    } else {
      comm_recv(&data->comm, &r, sizeof(double), other, tag);
      v += r;
    }
    ss = data->sep_size[lvl + 1];
    if (ss >= size || lvl == nl - 1)
      break;
    size -= ss;
  }
  /* fan-out */
  for (;;) {
    sint other = data->pother[lvl];
    if (other < 0) {
      comm_recv(&data->comm, &v, sizeof(double), -other - 1, tag);
    } else {
      comm_isend(&data->req[nsend++], &data->comm, &v, sizeof(double), other,
                 tag);
    }
    if (lvl == 0)
      break;
    ss = data->sep_size[lvl];
    size += ss, --lvl;
  }
  if (nsend)
    comm_wait(data->req, nsend);
  return v;
}

/* sorts an array of ids, removes 0's and duplicates;
   just returns the permutation */
#define unique_ids SUFFIXED_NAME(unique_ids)
static uint unique_ids(uint n, const ulong *id, sint *perm, buffer *buf) {
  uint *p, i, un = 0;
  ulong last = 0;
  p = sortp_long(buf, 0, id, n, sizeof(ulong));
  for (i = 0; i < n; ++i) {
    uint j = p[i];
    ulong v = id[j];
    if (v == 0)
      perm[j] = -1;
    else {
      if (v != last)
        last = v, ++un;
      perm[j] = un - 1;
    }
  }
  buf->n = 0;
  return un;
}

/* given user's list of dofs (as id's)
   uses gather-scatter to find share-count and separator # for each
   outputs as a list, sorted topologically (children before parents)
                      according to the sep. tree (and without duplicates),
           as well as the permutation to get there from the user's list */
/* sets: un, cn, perm_u2c */
#define discover_dofs SUFFIXED_NAME(discover_dofs)
static void discover_dofs(struct xxt *data, uint n, const ulong *id,
                          struct array *dofa, buffer *buf,
                          const struct comm *comm) {
  const uint pcoord = data->pcoord, ns = data->nsep;
  sint *perm;
  uint i, cn, *p, *pi;
  ulong *bid;
  struct gs_data *gsh;
  sint *v;
  struct dof *dof;

  data->un = n;
  data->perm_u2c = perm = tmalloc(sint, n);
  data->cn = cn = unique_ids(n, id, perm, buf);
  array_init(struct dof, dofa, cn), dofa->n = cn, dof = (struct dof *)dofa->ptr;
  buffer_reserve(buf, cn * sizeof(ulong)), bid = (ulong *)buf->ptr;
  for (i = 0; i < n; ++i)
    if (perm[i] >= 0)
      bid[perm[i]] = dof[perm[i]].id = id[i];

  gsh = gs_setup((const slong *)bid, cn, comm, 0, gs_crystal_router, 0);
  v = tmalloc(sint, cn);

  for (i = 0; i < cn; ++i)
    v[i] = pcoord;
  gs(v, gs_sint, gs_bpr, 0, gsh, buf);
  for (i = 0; i < cn; ++i)
    dof[i].level = ns - 1 - lg((uint)v[i]);

  for (i = 0; i < cn; ++i)
    v[i] = 1;
  gs(v, gs_sint, gs_add, 0, gsh, buf);
  for (i = 0; i < cn; ++i)
    dof[i].count = v[i];

  free(v);
  gs_free(gsh);

  if (!cn)
    return;
  buffer_reserve(buf, 2 * cn * sizeof(uint));
  p = sortp(buf, 0, &dof[0].level, cn, sizeof(struct dof));
  pi = p + cn;
  for (i = 0; i < cn; ++i)
    pi[p[i]] = i;
  for (i = 0; i < n; ++i)
    if (perm[i] >= 0)
      perm[i] = pi[perm[i]];
  sarray_permute_buf(struct dof, dof, cn, buf);
}

/* vl += A_ls * vs */
#define apply_p_Als SUFFIXED_NAME(apply_p_Als)
static void apply_p_Als(T *vl, struct xxt *data, const T *vs, uint ns) {
  const uint *Arp = data->A_sl.Arp, *Aj = data->A_sl.Aj;
  const T *A = data->A_sl.A;
  uint i, p, pe;
  for (i = 0; i < ns; ++i)
    for (p = Arp[i], pe = Arp[i + 1]; p != pe; ++p)
      vl[Aj[p]] += A[p] * vs[i];
}

/* vs -= A_sl * vl */
#define apply_m_Asl SUFFIXED_NAME(apply_m_Asl)
static void apply_m_Asl(T *vs, uint ns, struct xxt *data, const T *vl) {
  const uint *Arp = data->A_sl.Arp, *Aj = data->A_sl.Aj;
  const T *A = data->A_sl.A;
  uint i, p, pe;
  for (i = 0; i < ns; ++i)
    for (p = Arp[i], pe = Arp[i + 1]; p != pe; ++p)
      vs[i] -= A[p] * vl[Aj[p]];
}

/* returns a column of S : vs = -S(0:ei-1,ei) */
#define apply_S_col SUFFIXED_NAME(apply_S_col)
static void apply_S_col(T *vs, struct xxt *data, struct csr_mat *A_ss, uint ei,
                        T *vl) {
  const uint ln = data->ln;
  const uint *Asl_rp = data->A_sl.Arp, *Ass_rp = A_ss->Arp,
             *Asl_j = data->A_sl.Aj, *Ass_j = A_ss->Aj;
  const T *Ass = A_ss->A, *Asl = data->A_sl.A;
  uint i, p, pe;
  for (i = 0; i < ei; ++i)
    vs[i] = 0;
  for (p = Ass_rp[ei], pe = Ass_rp[ei + 1]; p != pe; ++p) {
    uint j = Ass_j[p];
    if (j >= ei)
      break;
    vs[j] = -Ass[p];
  }
  for (i = 0; i < ln; ++i)
    vl[i] = 0;
  for (p = Asl_rp[ei], pe = Asl_rp[ei + 1]; p != pe; ++p)
    vl[Asl_j[p]] = -Asl[p];
  sparse_cholesky_solve(vl, &data->fac_A_ll, vl);
  apply_m_Asl(vs, ei, data, vl);
}

#define apply_S SUFFIXED_NAME(apply_S)
static void apply_S(T *Svs, uint ns, struct xxt *data, struct csr_mat *A_ss,
                    const T *vs, T *vl) {
  const uint ln = data->ln;
  const uint *Ass_rp = A_ss->Arp, *Ass_j = A_ss->Aj;
  const T *Ass = A_ss->A;
  uint i, p, pe;
  for (i = 0; i < ns; ++i) {
    T sum = 0;
    for (p = Ass_rp[i], pe = Ass_rp[i + 1]; p != pe; ++p) {
      uint j = Ass_j[p];
      if (j >= ns)
        break;
      sum += Ass[p] * vs[j];
    }
    Svs[i] = sum;
  }
  for (i = 0; i < ln; ++i)
    vl[i] = 0;
  apply_p_Als(vl, data, vs, ns);
  sparse_cholesky_solve(vl, &data->fac_A_ll, vl);
  apply_m_Asl(Svs, ns, data, vl);
}

/* vx = X' * vs */
#define apply_Xt SUFFIXED_NAME(apply_Xt)
static void apply_Xt(T *vx, uint nx, const struct xxt *data, const T *vs) {
  const T *X = data->X;
  const uint *Xp = data->Xp;
  uint i;
  for (i = 0; i < nx; ++i)
    vx[i] = tensor_dot(vs, X + Xp[i], Xp[i + 1] - Xp[i]);
}

/* vs = X * vx */
#define apply_X SUFFIXED_NAME(apply_X)
static void apply_X(T *vs, uint ns, const struct xxt *data, const T *vx,
                    uint nx) {
  const T *X = data->X;
  const uint *Xp = data->Xp;
  uint i, j;
  for (i = 0; i < ns; ++i)
    vs[i] = 0;
  for (i = 0; i < nx; ++i) {
    const T v = vx[i];
    const T *x = X + Xp[i];
    uint n = Xp[i + 1] - Xp[i];
    for (j = 0; j < n; ++j)
      vs[j] += x[j] * v;
  }
}

#define allocate_X SUFFIXED_NAME(allocate_X)
static void allocate_X(struct xxt *data, sint *perm_x2c) {
  uint xn = data->xn;
  uint i, h = 0;
  if (data->null_space && xn)
    --xn;
  data->Xp = tmalloc(uint, xn + 1);
  data->Xp[0] = 0;
  for (i = 0; i < xn; ++i) {
    if (perm_x2c[i] != -1)
      ++h;
    data->Xp[i + 1] = data->Xp[i] + h;
  }
  data->X = tmalloc(T, data->Xp[xn]);
}

#define orthogonalize SUFFIXED_NAME(orthogonalize)
static void orthogonalize(struct xxt *data, struct csr_mat *A_ss,
                          sint *perm_x2c, buffer *buf) {
  uint ln = data->ln, sn = data->sn, xn = data->xn;
  T *vl, *vs, *vx, *Svs;
  uint i, j;

  allocate_X(data, perm_x2c);

  buffer_reserve(buf, (ln + 2 * sn + xn) * sizeof(T));
  vl = (T *)buf->ptr, vs = vl + ln, Svs = vs + sn, vx = Svs + sn;

  if (data->null_space && xn)
    --xn;
  for (i = 0; i < xn; ++i) {
    uint ns = data->Xp[i + 1] - data->Xp[i];
    sint ui = perm_x2c[i];
    T ytsy, *x;

    if (ui == -1) {
      for (j = 0; j < i; ++j)
        vx[j] = 0;
    } else {
      ui -= ln;
      apply_S_col(vs, data, A_ss, ui, vl);
      apply_Xt(vx, i, data, vs);
    }
    apply_QQt(data, vx, i, xn - i);
    apply_X(vs, ns, data, vx, i);
    if (ui != -1)
      vs[ui] = 1;
    apply_S(Svs, ns, data, A_ss, vs, vl);
    ytsy = tensor_dot(vs, Svs, ns);
    ytsy = sum(data, ytsy, i + 1, xn - (i + 1));
    if (ytsy < DBL_EPSILON / 128)
      ytsy = 0;
    else
      ytsy = 1 / sqrt(ytsy);
    x = &data->X[data->Xp[i]];
    for (j = 0; j < ns; ++j)
      x[j] = ytsy * vs[j];
  }
}

/* produces CSR matrix from Yale-like format, summing duplicates */
#define condense_matrix SUFFIXED_NAME(condense_matrix)
static void condense_matrix(struct array *mat, uint nr, struct csr_mat *out,
                            buffer *buf) {
  uint k, nz = mat->n;
  struct yale_mat *p, *q;
  sarray_sort_2(struct yale_mat, mat->ptr, mat->n, i, 0, j, 0, buf);

  p = (struct yale_mat *)mat->ptr;
  for (k = 0; k + 1 < nz; ++k, ++p)
    if (p[0].i == p[1].i && p[0].j == p[1].j)
      break;
  if (++k < nz) {
    uint i = p->i, j = p->j;
    q = p + 1;
    for (; k < nz; ++k, ++q) {
      if (i == q->i && j == q->j)
        p->v += q->v, --mat->n;
      else
        ++p, p->i = i = q->i, p->j = j = q->j, p->v = q->v;
    }
  }

  nz = mat->n;
  out->n = nr;
  out->Arp = tmalloc(uint, nr + 1 + mat->n);
  out->Aj = out->Arp + nr + 1;
  T *A = out->A = tmalloc(T, nz);
  for (k = 0; k < nr; ++k)
    out->Arp[k] = 0;
  for (p = (struct yale_mat *)mat->ptr, k = 0; k < nz; ++k, ++p)
    out->Arp[p->i]++, out->Aj[k] = p->j, A[k] = p->v;
  nz = 0;
  for (k = 0; k <= nr; ++k) {
    uint t = out->Arp[k];
    out->Arp[k] = nz, nz += t;
  }
}

#define separate_matrix SUFFIXED_NAME(separate_matrix)
static void separate_matrix(uint nz, const uint *Ai, const uint *Aj,
                            const double *A, const sint *perm, uint ln, uint sn,
                            struct csr_mat *out_ll, struct csr_mat *out_sl,
                            struct csr_mat *out_ss, buffer *buf) {
  uint k, n;
  struct array mat_ll, mat_sl, mat_ss;
  struct yale_mat *mll, *msl, *mss;
  array_init(struct yale_mat, &mat_ll, 2 * nz),
      mll = (struct yale_mat *)mat_ll.ptr;
  array_init(struct yale_mat, &mat_sl, 2 * nz),
      msl = (struct yale_mat *)mat_sl.ptr;
  array_init(struct yale_mat, &mat_ss, 2 * nz),
      mss = (struct yale_mat *)mat_ss.ptr;
  for (k = 0; k < nz; ++k) {
    sint i = perm[Ai[k]], j = perm[Aj[k]];
    if (i < 0 || j < 0 || A[k] == 0)
      continue;
    if ((uint)i < ln) {
      if ((uint)j < ln)
        n = mat_ll.n++, mll[n].i = i, mll[n].j = j, mll[n].v = A[k];
    } else {
      if ((uint)j < ln)
        n = mat_sl.n++, msl[n].i = i - ln, msl[n].j = j, msl[n].v = A[k];
      else
        n = mat_ss.n++, mss[n].i = i - ln, mss[n].j = j - ln, mss[n].v = A[k];
    }
  }
  condense_matrix(&mat_ll, ln, out_ll, buf);
  condense_matrix(&mat_sl, sn, out_sl, buf);
  condense_matrix(&mat_ss, sn, out_ss, buf);
  array_free(&mat_ll);
  array_free(&mat_sl);
  array_free(&mat_ss);
}

#define crs_xxt_setup SUFFIXED_NAME(crs_xxt_setup)
struct xxt *crs_xxt_setup(uint n, const ulong *id, uint nz, const uint *Ai,
                          const uint *Aj, const double *A, uint null_space,
                          const struct comm *comm) {
  struct xxt *data = tmalloc(struct xxt, 1);
  sint *perm_x2c;
  struct array dofa;
  struct csr_mat A_ll, A_ss;
  buffer buf;

  comm_dup(&data->comm, comm);

  locate_proc(data);

  data->null_space = null_space;

  buffer_init(&buf, 1024);

  discover_dofs(data, n, id, &dofa, &buf, &data->comm);
  discover_sep_sizes(data, &dofa, &buf);

  perm_x2c = discover_sep_ids(data, &dofa, &buf);
  if (data->null_space) {
    uint i;
    double count = 0;
    struct dof *dof = (struct dof *)dofa.ptr;
    for (i = 0; i < data->cn; ++i)
      count += 1 / (double)dof[i].count;
    count = 1 / sum(data, count, data->xn, 0);
    data->share_weight = tcalloc(T, data->cn);
    for (i = 0; i < data->cn; ++i)
      data->share_weight[i] = count / dof[i].count;
  }
  array_free(&dofa);

  if (!data->null_space || data->xn != 0) {
    separate_matrix(nz, Ai, Aj, A, data->perm_u2c, data->ln, data->sn, &A_ll,
                    &data->A_sl, &A_ss, &buf);
  } else {
    separate_matrix(nz, Ai, Aj, A, data->perm_u2c, data->ln - 1, 1, &A_ll,
                    &data->A_sl, &A_ss, &buf);
  }

  sparse_cholesky_factor(A_ll.n, A_ll.Arp, A_ll.Aj, A_ll.A, &data->fac_A_ll,
                         &buf);
#if defined(ENABLE_BLAS)
  A_ll_inv = (double *)memalign(64, sizeof(double) * A_ll.n * A_ll.n);
  for (uint i = 0; i < A_ll.n; i++) {
    for (uint j = A_ll.Arp[i]; j < A_ll.Arp[i + 1]; j++)
      A_ll_inv[i * A_ll.n + A_ll.Aj[j]] = A_ll.A[j];
  }

  int *ipiv = tmalloc(int, A_ll.n);
  int info =
      LAPACKE_dgetrf(LAPACK_ROW_MAJOR, A_ll.n, A_ll.n, A_ll_inv, A_ll.n, ipiv);
  assert(info == 0);
  info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, A_ll.n, A_ll_inv, A_ll.n, ipiv);
  assert(info == 0);
  free(ipiv);

  A_ll_inv_f32 = (float *)memalign(64, sizeof(float) * A_ll.n * A_ll.n);
  for (uint i = 0; i < A_ll.n; i++) {
    for (uint j = 0; j < A_ll.n; j++)
      A_ll_inv_f32[i * A_ll.n + j] = (float)A_ll_inv[i * A_ll.n + j];
  }

  y_inv = (double *)memalign(64, sizeof(double) * A_ll.n);
  y_inv_f32 = (float *)memalign(64, sizeof(float) * A_ll.n);
#endif

  free(A_ll.Arp);
  free(A_ll.A);

  data->vl = tmalloc(T, data->ln + data->cn + 2 * data->xn);
  data->vc = data->vl + data->ln;
  data->vx = data->vc + data->cn;
  data->combuf = data->vx + data->xn;

  orthogonalize(data, &A_ss, perm_x2c, &buf);
  free(A_ss.Arp);
  free(A_ss.A);
  free(perm_x2c);
  buffer_free(&buf);

  return data;
}

#define crs_xxt_solve SUFFIXED_NAME(crs_xxt_solve)
void crs_xxt_solve(T *x, struct xxt *data, const T *b) {
  uint cn = data->cn, un = data->un, ln = data->ln, sn = data->sn;
  uint xn = data->xn;

  MPI_Barrier(data->comm.c);
  double t = MPI_Wtime();
  T *vl = data->vl, *vc = data->vc, *vx = data->vx;
  uint i;
  for (i = 0; i < cn; ++i)
    vc[i] = 0;
  for (i = 0; i < un; ++i) {
    sint p = data->perm_u2c[i];
    if (p >= 0)
      vc[p] += b[i];
  }
  local_time += MPI_Wtime() - t;

  if (xn > 0 && (!data->null_space || xn > 1)) {
    if (data->null_space)
      --xn;

    t = MPI_Wtime();
#if defined(ENABLE_BLAS)
    if (sizeof(T) == sizeof(double)) {
      cblas_dgemv(CblasRowMajor, CblasNoTrans, ln, ln, 1.0, A_ll_inv, ln,
                  (double *)vc, 1, 0.0, y_inv, 1);
      for (uint i = 0; i < ln; ++i)
        vc[i] = y_inv[i];
    } else {
      cblas_sgemv(CblasRowMajor, CblasNoTrans, ln, ln, 1.0, A_ll_inv_f32, ln,
                  (float *)vc, 1, 0.0, y_inv_f32, 1);
      for (uint i = 0; i < ln; ++i)
        vc[i] = y_inv_f32[i];
    }
#else
    sparse_cholesky_solve(vc, &data->fac_A_ll, vc);
#endif
    cholesky_time += MPI_Wtime() - t;

    t = MPI_Wtime();
    apply_m_Asl(vc + ln, sn, data, vc);
    local_time += MPI_Wtime() - t;

    t = MPI_Wtime();
    apply_Xt(vx, xn, data, vc + ln);
    xxt_time += MPI_Wtime() - t;

    t = MPI_Wtime();
    apply_QQt(data, vx, xn, 0);
    qqt_time += MPI_Wtime() - t;

    t = MPI_Wtime();
    apply_X(vc + ln, sn, data, vx, xn);
    xxt_time += MPI_Wtime() - t;

    t = MPI_Wtime();
    for (i = 0; i < ln; ++i)
      vl[i] = 0;
    apply_p_Als(vl, data, vc + ln, sn);
    local_time += MPI_Wtime() - t;

    t = MPI_Wtime();
#if defined(ENABLE_BLAS)
    if (sizeof(T) == sizeof(double)) {
      cblas_dgemv(CblasRowMajor, CblasNoTrans, ln, ln, 1.0, A_ll_inv, ln,
                  (double *)vl, 1, 0.0, y_inv, 1);
      for (uint i = 0; i < ln; ++i)
        vl[i] = y_inv[i];
    } else {
      cblas_sgemv(CblasRowMajor, CblasNoTrans, ln, ln, 1.0, A_ll_inv_f32, ln,
                  (float *)vl, 1, 0.0, y_inv_f32, 1);
      for (uint i = 0; i < ln; ++i)
        vl[i] = y_inv_f32[i];
    }
#else
    sparse_cholesky_solve(vl, &data->fac_A_ll, vl);
#endif
    cholesky_time += MPI_Wtime() - t;

    t = MPI_Wtime();
    for (i = 0; i < ln; ++i)
      vc[i] -= vl[i];
    local_time += MPI_Wtime() - t;
  } else {
#if defined(ENABLE_BLAS)
    if (sizeof(T) == sizeof(double)) {
      cblas_dgemv(CblasRowMajor, CblasNoTrans, ln, ln, 1.0, A_ll_inv, ln,
                  (double *)vc, 1, 0.0, y_inv, 1);
      for (uint i = 0; i < ln; ++i)
        vc[i] = y_inv[i];
    } else {
      cblas_sgemv(CblasRowMajor, CblasNoTrans, ln, ln, 1.0, A_ll_inv_f32, ln,
                  (float *)vc, 1, 0.0, y_inv_f32, 1);
      for (uint i = 0; i < ln; ++i)
        vc[i] = y_inv_f32[i];
    }
#else
    sparse_cholesky_solve(vc, &data->fac_A_ll, vc);
#endif
    if (data->null_space) {
      if (xn == 0)
        vc[ln - 1] = 0;
      else if (sn == 1)
        vc[ln] = 0;
    }
  }
  if (data->null_space) {
    T s = 0;
    for (i = 0; i < cn; ++i)
      s += data->share_weight[i] * vc[i];
    s = sum(data, s, data->xn, 0);
    for (i = 0; i < cn; ++i)
      vc[i] -= s;
  }

  t = MPI_Wtime();
  for (i = 0; i < un; ++i) {
    sint p = data->perm_u2c[i];
    x[i] = p >= 0 ? vc[p] : 0;
  }
  local_time += MPI_Wtime() - t;
}

#define crs_xxt_stats SUFFIXED_NAME(crs_xxt_stats)
void crs_xxt_stats(struct xxt *data) {
  int a, b;
  slong in[1], out[2][1], buf[2][1];
  uint xcol;
  if (data->comm.id == 0) {
    unsigned s;
    printf("xxt: separator sizes on %d =", (int)data->comm.id);
    for (s = 0; s < data->nsep; ++s)
      printf(" %d", (int)data->sep_size[s]);
    printf("\n");
    printf("xxt: shared dofs on %d = %d\n", (int)data->comm.id, (int)data->sn);
  }

  a = data->ln;
  comm_allreduce(&data->comm, gs_int, gs_max, &a, 1, &b);
  if (data->comm.id == 0)
    printf("xxt: max non-shared dofs = %d\n", a);
  in[0] = data->ln;
  comm_scan(out, &data->comm, gs_long, gs_add, in, 1, buf);
  if (data->comm.id == 0)
    printf("xxt: total non-shared dofs = %lld\n", out[1][0]);

  a = data->sn;
  comm_allreduce(&data->comm, gs_int, gs_max, &a, 1, &b);
  if (data->comm.id == 0)
    printf("xxt: max shared dofs = %d\n", a);
  in[0] = data->sn;
  comm_scan(out, &data->comm, gs_long, gs_add, in, 1, buf);
  if (data->comm.id == 0)
    printf("xxt: total shared dofs = %lld\n", out[1][0]);

  xcol = data->xn;
  if (xcol && data->null_space)
    --xcol;
  a = xcol;
  comm_allreduce(&data->comm, gs_int, gs_max, &a, 1, &b);
  if (data->comm.id == 0)
    printf("xxt: max X cols = %d\n", a);
  in[0] = data->xn;
  comm_scan(out, &data->comm, gs_long, gs_add, in, 1, buf);
  if (data->comm.id == 0)
    printf("xxt: total X cols = %lld\n", out[1][0] - data->null_space);
  slong gxcol = out[1][0];

  a = data->Xp[xcol];
  comm_allreduce(&data->comm, gs_int, gs_max, &a, 1, &b);
  if (data->comm.id == 0)
    printf("xxt: max X nnz = %d\n", a);
  in[0] = data->Xp[xcol];
  comm_scan(out, &data->comm, gs_long, gs_add, in, 1, buf);
  if (data->comm.id == 0) {
    printf("xxt: total X nnz = %lld\n", out[1][0]);
    printf("xxt: row X nnz = %lld\n", out[1][0] / gxcol);
  }
  fflush(stdout);
}

#define crs_xxt_free SUFFIXED_NAME(crs_xxt_free)
void crs_xxt_free(struct xxt *data) {
  comm_free(&data->comm);
  free(data->pother);
  free(data->req);
  free(data->sep_size);
  free(data->perm_u2c);
  if (data->null_space)
    free(data->share_weight);
  sparse_cholesky_free(&data->fac_A_ll);
  free(data->A_sl.Arp);
  free(data->A_sl.A);
  free(data->Xp);
  free(data->X);
  free(data->vl);
  free(data);

#if defined(ENABLE_BLAS)
  free(A_ll_inv);
  free(A_ll_inv_f32);
  free(y_inv);
  free(y_inv_f32);
#endif
}

#undef crs_xxt_setup
#undef crs_xxt_solve
#undef crs_xxt_stats
#undef crs_xxt_free
#undef xxt
#undef csr_mat
#undef SUFFIXED_NAME

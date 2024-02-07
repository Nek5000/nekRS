#ifdef SUFFIX
#define SUFFIXED_NAME(x) TOKEN_PASTE(x, SUFFIX)
#else
#define SUFFIXED_NAME(x) x
#endif

#define sparse_cholmod_factor SUFFIXED_NAME(sparse_cholmod_factor)
static struct cholmod_csr *sparse_cholmod_factor(uint n, const uint *Arp,
                                                 const uint *Aj, const TT *A,
                                                 buffer *bfr) {
  struct cholmod_csr *B = tcalloc(struct cholmod_csr, 1);

  cholmod_start(&B->cm);
  B->cm.itype = CHOLMOD_INT;
  B->cm.dtype = PRECISION;
  B->cm.error_handler = &error_handler;
  B->cm.postorder = 0;
  B->nr = n;

  const uint nnz = Arp[n];
  cholmod_triplet *T =
      cholmod_allocate_triplet(n, n, nnz, -1, CHOLMOD_REAL, &B->cm);

  int32_t *Ti = (int32_t *)T->i, *Tj = (int32_t *)T->j;
  TT *Tx = (TT *)T->x;

  uint nnz_ = 0;
  for (uint i = 0; i < n; i++) {
    uint j, je = Arp[i + 1];
    for (j = Arp[i]; j < je && Aj[j] < i; j++)
      ;
    for (; j < je; j++, nnz_++) {
      Ti[nnz_] = i;
      Tj[nnz_] = Aj[j];
      Tx[nnz_] = A[j];
    }
  }
  T->nnz = nnz_;

  B->A = cholmod_triplet_to_sparse(T, T->nnz, &B->cm);
  cholmod_free_triplet(&T, &B->cm);

  B->cm.nmethods = 1;
  const int methods[3] = {CHOLMOD_NATURAL, CHOLMOD_AMD, CHOLMOD_NESDIS};
  uint min_nnz = UINT_MAX, min_method = 0;
  for (uint i = 0; i < 3; i++) {
    B->cm.method[0].ordering = methods[i];
    B->L = cholmod_analyze(B->A, &B->cm);
    cholmod_factorize(B->A, B->L, &B->cm);
    cholmod_change_factor(CHOLMOD_REAL, 0, 0, 0, 0, B->L, &B->cm);
    if (B->L->nzmax < min_nnz) {
      min_nnz = B->L->nzmax;
      min_method = i;
    }
    cholmod_free_factor(&B->L, &B->cm);
  }

  B->cm.method[0].ordering = methods[min_method];
  B->L = cholmod_analyze(B->A, &B->cm);
  cholmod_factorize(B->A, B->L, &B->cm);
  cholmod_change_factor(CHOLMOD_REAL, 0, 0, 0, 0, B->L, &B->cm);

  B->r = cholmod_zeros(n, 1, CHOLMOD_REAL, &B->cm);

  return B;
}
#undef sparse_cholmod_factor

#define sparse_cholmod_solve SUFFIXED_NAME(sparse_cholmod_solve)
static void sparse_cholmod_solve(TT *x, struct cholmod_csr *B, const TT *r) {

  TT *rx = (TT *)B->r->x;
  for (uint i = 0; i < B->nr; i++)
    rx[i] = r[i];

  cholmod_dense *xd = cholmod_solve(CHOLMOD_A, B->L, B->r, &B->cm);

  TT *xx = (TT *)xd->x;
  for (uint i = 0; i < B->nr; i++)
    x[i] = xx[i];

  cholmod_free_dense(&xd, &B->cm);
}
#undef sparse_cholmod_solve

#define sparse_cholmod_free SUFFIXED_NAME(sparse_cholmod_free)
static void sparse_cholmod_free(struct cholmod_csr *factor) {
  if (!factor)
    return;
  cholmod_free_sparse(&factor->A, &factor->cm);
  cholmod_free_factor(&factor->L, &factor->cm);
  cholmod_free_dense(&factor->r, &factor->cm);
  cholmod_finish(&factor->cm);
  free(factor);
}
#undef sparse_cholmod_free

#define csr_setup SUFFIXED_NAME(csr_setup)
static struct cholmod_csr *csr_setup(struct csr *A, const int null,
                                     const struct comm *const comm) {
  struct cholmod_csr *B = tcalloc(struct cholmod_csr, 1);

  cholmod_start(&B->cm);
  B->cm.itype = CHOLMOD_INT;
  B->cm.dtype = PRECISION;
  B->cm.error_handler = &error_handler;
  B->cm.postorder = 0;
  B->nr = A->nr;

  uint nnz = A->offs[A->nr];
  cholmod_triplet *T =
      cholmod_allocate_triplet(A->nr, A->nr, nnz, -1, CHOLMOD_REAL, &B->cm);
  int32_t *Ti = (int32_t *)T->i, *Tj = (int32_t *)T->j;

  uint z = 0;
  TT *Tx = (TT *)T->x;
  for (uint i = 0; i < A->nr; i++) {
    uint j;
    for (j = A->offs[i]; A->cols[j] < i; j++)
      ;
    for (uint je = A->offs[i + 1]; j < je; j++)
      Ti[z] = i, Tj[z] = A->cols[j], Tx[z] = A->vals[j], z++;
  }
  T->nnz = z;

  // Convert triplet to CSC matrix.
  B->A = cholmod_triplet_to_sparse(T, T->nnz, &B->cm);
  cholmod_free_triplet(&T, &B->cm);

  B->cm.nmethods = 1;
  const int methods[3] = {CHOLMOD_NATURAL, CHOLMOD_AMD, CHOLMOD_NESDIS};

  uint min_nnz = UINT_MAX, min_method = 0;
  for (uint i = 0; i < 3; i++) {
    B->cm.method[0].ordering = methods[i];
    B->L = cholmod_analyze(B->A, &B->cm);
    cholmod_factorize(B->A, B->L, &B->cm);
    cholmod_change_factor(CHOLMOD_REAL, 0, 0, 0, 0, B->L, &B->cm);
    if (B->L->nzmax < min_nnz) {
      min_nnz = B->L->nzmax;
      min_method = i;
    }
    cholmod_free_factor(&B->L, &B->cm);
  }

  B->cm.method[0].ordering = methods[min_method];
  B->L = cholmod_analyze(B->A, &B->cm);
  cholmod_factorize(B->A, B->L, &B->cm);
  cholmod_change_factor(CHOLMOD_REAL, 0, 0, 0, 0, B->L, &B->cm);

  // Find min, max and average number of nonzeros per row.
  sint min = B->L->nzmax, max = B->L->nzmax;
  slong avg = B->L->nzmax, wrk;
  comm_allreduce(comm, gs_int, gs_min, &min, 1, &wrk);
  comm_allreduce(comm, gs_int, gs_max, &max, 1, &wrk);
  comm_allreduce(comm, gs_long, gs_add, &avg, 1, &wrk);
  if (comm->id == 0) {
    printf("\tBOX CHOLMOD min_nnz = %d max_nnz = %d avg_nnz = %ld\n", min, max,
           avg / comm->np);
    fflush(stdout);
  }

  // Allocate space for RHS.
  B->r = cholmod_zeros(A->nr, 1, CHOLMOD_REAL, &B->cm);

  return B;
}

#define solve SUFFIXED_NAME(solve)
static void solve(TT *x, struct box *const box, const TT *r) {
  struct cholmod_csr *B = (struct cholmod_csr *)box->ss;

  TT *rx = (TT *)B->r->x;
  for (uint i = 0; i < B->nr; i++)
    rx[i] = 0;

  for (uint i = 0; i < box->sn; i++) {
    if (box->u2c[i] >= 0)
      rx[box->u2c[i]] += r[i];
  }

  cholmod_dense *xd = cholmod_solve(CHOLMOD_A, B->L, B->r, &B->cm);

#if 0
  double one[2] = {1, 0}, m1[2] = {-1, 0};
  cholmod_dense *rd = cholmod_copy_dense(B->r, &B->cm);
  cholmod_sdmult(B->A, 0, m1, one, xd, rd, &B->cm);
  printf("rank = %d nr = %d norm(b-Ax) = %e\n", box->global.id, B->nr,
         cholmod_norm_dense(rd, 0, &B->cm));
  fflush(stdout);
  cholmod_free_dense(&rd, &B->cm);
#endif

  TT *xx = (TT *)xd->x;
  for (unsigned i = 0; i < box->sn; i++) {
    if (box->u2c[i] >= 0)
      x[i] = xx[box->u2c[i]];
    else
      x[i] = 0;
  }

  cholmod_free_dense(&xd, &B->cm);
}

#undef solve
#undef csr_setup

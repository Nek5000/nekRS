//=============================================================================
// Test Schur complement solver
//
#include "coarse.h"
#include "parRSB.h"

#include <math.h>
#include <time.h>

static double check_err(double *b, double *x, uint nelt, uint nv,
                        const slong *vtx, MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  slong out[2][1], buf[2][1], in = nelt;
  comm_scan(out, &c, gs_long, gs_add, &in, 1, buf);
  ulong start = out[0][0] + 1;

  ulong *eid = tcalloc(ulong, nelt);
  for (uint i = 0; i < nelt; i++)
    eid[i] = start + i;

  struct crystal cr;
  crystal_init(&cr, &c);

  buffer bfr;
  buffer_init(&bfr, 1024);

  struct array nbrs, eij;
  find_nbrs(&nbrs, eid, vtx, nelt, nv, &cr, &bfr);
  compress_nbrs(&eij, &nbrs, &bfr);

  struct par_mat M;
  par_csr_setup(&M, &eij, 1, &bfr);
  assert(M.rn > 0);

  free(eid), array_free(&nbrs), array_free(&eij);

  struct gs_data *gsh = setup_Q(&M, &c, &bfr);
  double *bl = tcalloc(double, nelt);
  double *wrk = tcalloc(double, M.rn + M.adj_off[M.rn]);
  mat_vec_csr(bl, x, &M, gsh, wrk, &bfr);

  crystal_free(&cr), comm_free(&c);
  gs_free(gsh), par_mat_free(&M);

  double norm = 0.0;
  for (uint i = 0; i < nelt; i++)
    norm += (bl[i] - b[i]) * (bl[i] - b[i]);
  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, comm);

  free(wrk), free(bl);
  buffer_free(&bfr);

  return sqrt(norm);
}

static void setup_rhs(double *b, const unsigned int nelt, MPI_Comm comm) {
  srand(time(NULL));
  double sum = 0;
  for (int i = 0; i < nelt; i++) {
    b[i] = (rand() % 50 + 1.0) / 10;
    sum += b[i];
  }
  MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, comm);

  long long ng = nelt;
  MPI_Allreduce(MPI_IN_PLACE, &ng, 1, MPI_LONG_LONG, MPI_SUM, comm);
  sum /= ng;

  double norm = 0;
  for (int i = 0; i < nelt; i++) {
    b[i] -= sum;
    norm += b[i] * b[i];
  }

  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, comm);
  norm = sqrt(norm);

  for (int i = 0; i < nelt; i++)
    b[i] /= norm;
}

static void setup_and_solve(unsigned nelt, unsigned nv, const long long *vl,
                            const scalar *centroids,
                            const parrsb_cmd_line_opts *in, MPI_Comm comm) {
  // Setup the coarse solve with schur complement solver
  struct comm c;
  comm_init(&c, comm);

  comm_barrier(&c);
  double t = comm_time();
  struct coarse *crs =
      coarse_setup(nelt, nv, vl, centroids, 1, in->crs_type, &c);
  double tsetup = comm_time() - t;

  scalar *b = tcalloc(scalar, 2 * nelt);
  setup_rhs(b, nelt, comm);

  comm_barrier(&c);
  t = comm_time();
  scalar *x = b + nelt;
  coarse_solve(x, crs, b, in->crs_tol);
  double tsolve = MPI_Wtime() - t;

  double enorm = check_err(b, x, nelt, nv, vl, comm);
  if (c.id == 0) {
    printf("MPI Ranks = %d\ncoarse_setup: %lf\ncoarse_solve = %lf\nerr = %lf\n",
           c.np, tsetup, tsolve, enorm);
    fflush(stdout);
  }
  int err = (enorm > 10 * in->crs_tol);
  parrsb_check_error(err, comm);

  // Free resources
  coarse_free(crs), free(b);
  comm_free(&c);
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  parrsb_cmd_line_opts *in = parrsb_parse_cmd_opts(argc, argv);
  parrsb_check_error(in == NULL, comm);

  // Read the geometry from the .re2 file, find connectiviy, partition and then
  // distribute the mesh.
  unsigned nelt, nv;
  long long *vl = NULL;
  double *coord = NULL;
  int err = parrsb_setup_mesh(&nelt, &nv, &vl, &coord, in, comm);
  parrsb_check_error(err, comm);

  int ndim = (nv == 8) ? 3 : 2;
  double *centroids = tcalloc(double, nelt *ndim);
  for (uint i = 0; i < nelt; i++) {
    for (int j = 0; j < nv; j++) {
      for (int d = 0; d < ndim; d++)
        centroids[i * ndim + d] += coord[i * ndim * nv + j * ndim + d];
    }
    for (int d = 0; d < ndim; d++)
      centroids[i * ndim + d] /= nv;
  }

  setup_and_solve(nelt, nv, vl, centroids, in, comm);

  free(vl), free(coord), free(centroids);
  parrsb_cmd_opts_free(in);
  MPI_Finalize();

  return 0;
}

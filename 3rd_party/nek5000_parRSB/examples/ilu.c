//=============================================================================
// Test ILU factorization
//
#include "ilu.h"
#include "parRSB.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  struct parrsb_cmd_opts *in = parrsb_parse_cmd_opts(argc, argv);
  int err = (in == NULL);
  parrsb_check_error(err, comm);

  // Read the geometry from the .re2 file, find connectiviy, partition and then
  // distribute the mesh.
  unsigned int nelt, nv;
  long long *vl = NULL;
  double *coord = NULL;
  parrsb_setup_mesh(&nelt, &nv, &vl, &coord, in, comm);

  // Setup ILU
  ilu_options iluopt = {.type = in->ilu_type,
                        .tol = in->ilu_tol,
                        .pivot = in->ilu_pivot,
                        .verbose = in->verbose,
                        .nnz_per_row = 0};
  struct ilu *ilu = ilu_setup(nelt, nv, vl, &iluopt, comm);
  ilu_free(ilu);

  // Free resources
  free(vl), free(coord);
  free(in);
  MPI_Finalize();

  return 0;
}

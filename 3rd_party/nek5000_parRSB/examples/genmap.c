//=============================================================================
// Generate partitions (.ma2) from Nek5000's mesh (.re2) file.
//
#include "parRSB.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm world = MPI_COMM_WORLD;

  struct parrsb_cmd_opts *in = parrsb_parse_cmd_opts(argc, argv);
  int err = (in == NULL);
  parrsb_check_error(err, world);

  int rank, size;
  MPI_Comm_rank(world, &rank);
  MPI_Comm_size(world, &size);

  if (in->nactive > size)
    in->nactive = size;

  MPI_Comm comm;
  MPI_Comm_split(world, rank < in->nactive, rank, &comm);
  if (rank < in->nactive) {
    // Read the geometry from the .re2 file
    unsigned int nelt, nbcs, nv;
    double *coord = NULL;
    long long *bcs = NULL;
    err = parrsb_read_mesh(&nelt, &nv, NULL, &coord, &nbcs, &bcs, in->mesh,
                           comm, 1);
    parrsb_check_error(err, comm);

    // Find connectivity
    long long *vl = (long long *)calloc(nelt * nv, sizeof(long long));
    err = (vl == NULL);
    parrsb_check_error(err, comm);

    unsigned ndim = (nv == 8 ? 3 : 2);
    err = parrsb_conn_mesh(vl, coord, nelt, ndim, bcs, nbcs, in->tol, comm);
    parrsb_check_error(err, comm);

    // Print pre-partition statistics
    int nss[6];
    if (in->verbose > 0) {
      if (rank == 0) {
        printf("Partition statistics before RSB:\n");
        fflush(stdout);
      }
      parrsb_print_part_stat(vl, nelt, nv, comm);
    }
    parrsb_get_part_stat(NULL, NULL, nss, NULL, vl, nelt, nv, comm);

    // Partition the mesh
    int *part = (int *)calloc(nelt, sizeof(int));
    err = (part == NULL);
    parrsb_check_error(err, comm);

    parrsb_options options = parrsb_default_options;
    err = parrsb_part_mesh(part, NULL, vl, coord, nelt, nv, options, comm);
    parrsb_check_error(err, comm);

    // Redistribute data based on identified partitions
    err = parrsb_dist_mesh(&nelt, &vl, &coord, part, nv, comm);
    parrsb_check_error(err, comm);

    if (in->verbose > 0) {
      if (rank == 0) {
        printf("Partition statistics after RSB:\n");
        fflush(stdout);
      }
      parrsb_print_part_stat(vl, nelt, nv, comm);
    }
    parrsb_get_part_stat(NULL, NULL, &nss[3], NULL, vl, nelt, nv, comm);

    // Write partition to .ma2 file
    if (in->dump) {
      err = parrsb_dump_map(in->mesh, nelt, nv, vl, comm);
      parrsb_check_error(err, comm);
    }

    if (in->test && in->nactive > 1) {
      err = (nss[2] < nss[5]);
      parrsb_check_error(err, comm);
    }

    free(part), free(vl), free(coord), free(bcs);
  }

  // Free resources
  parrsb_cmd_opts_free(in);
  MPI_Comm_free(&comm);
  MPI_Finalize();

  return 0;
}

/*
 * Generate partitions (.ma2) from Nek5000's mesh (.re2) file.
 */
#include <parRSB.h>

static int test_parrsb() {}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int id;
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  parrsb_input *in = parrsb_parse_input(argc, argv);

  int active = 0;
  if (id < in->nactive)
    active = 1;

  MPI_Comm comm;
  MPI_Comm_split(MPI_COMM_WORLD, active, id, &comm);

  /* Read the geometry from the .re2 file */
  unsigned int nelt, nbcs;
  double *coord = NULL;
  long long *bcs = NULL;
  int nv;
  int err = 0;
  if (active == 1)
    err = parrsb_read_mesh(&nelt, &nv, NULL, &coord, &nbcs, &bcs, in->mesh,
                           comm, 1);
  parrsb_check_error(err, comm);

  /* Find connectivity */
  long long *vl = (long long *)calloc(nelt * nv, sizeof(long long));
  int ndim = nv == 8 ? 3 : 2;
  if (active == 1)
    err |= parrsb_find_conn(vl, coord, nelt, ndim, bcs, nbcs, in->tol, comm, 0);
  parrsb_check_error(err, comm);

  int nss[6];
  if (active == 1 && in->nactive > 1) {
    if (id == 0)
      printf("Partition statistics before RSB:\n");
    parrsb_print_part_stat(vl, nelt, nv, comm);
    parrsb_get_part_stat(NULL, NULL, &nss[0], NULL, vl, nelt, nv, comm);
  }

  /* Partition the mesh */
  parrsb_options options = parrsb_default_options;
  int *part = (int *)calloc(nelt, sizeof(int));
  if (active == 1)
    err |= parrsb_part_mesh(part, NULL, vl, coord, nelt, nv, options, comm);
  parrsb_check_error(err, comm);

  /* Redistribute data */
  if (active == 1)
    err |= parrsb_distribute_elements(&nelt, &vl, &coord, part, nv, comm);
  parrsb_check_error(err, comm);

  if (active == 1 && in->nactive > 1) {
    if (id == 0)
      printf("Partition statistics after RSB:\n");
    parrsb_print_part_stat(vl, nelt, nv, comm);
    parrsb_get_part_stat(NULL, NULL, &nss[3], NULL, vl, nelt, nv, comm);
  }

  if (active == 1 && in->test && in->nactive > 1)
    err |= nss[2] < nss[5];
  parrsb_check_error(err, comm);

  /* Write partition to .ma2 file */
  if (active == 1 && in->dump == 1)
    err |= parrsb_dump_map(in->mesh, nelt, nv, vl, part, comm);
  parrsb_check_error(err, comm);

  /* Free resources */
  if (part != NULL)
    free(part);
  if (vl != NULL)
    free(vl);
  if (coord != NULL)
    free(coord);
  if (bcs != NULL)
    free(bcs);
  if (in != NULL)
    free(in);

  MPI_Comm_free(&comm);
  MPI_Finalize();

  return 0;
}

/*
 * Generate connectivity (.co2) from Nek5000 mesh (.re2) file.
 */
#include <parRSB.h>

static int test_parcon(unsigned int neltp, long long *vlp, char *name,
                       MPI_Comm comm) {
  unsigned int nelt;
  int nv;
  long long *vls = NULL;
  int err = parrsb_read_mesh(&nelt, &nv, &vls, NULL, NULL, NULL, name,
                             MPI_COMM_WORLD, 2);
  assert(neltp == nelt);

  uint size = nelt * nv;
  slong *minp = tcalloc(slong, size);
  slong *maxp = tcalloc(slong, size);

  buffer buf;
  buffer_init(&buf, 1024);

  struct comm c;
  comm_init(&c, comm);

  struct gs_data *gsh = gs_setup(vls, size, &c, 0, gs_pairwise, 0);

  uint i;
  for (i = 0; i < size; i++)
    minp[i] = maxp[i] = vlp[i];

  gs(minp, gs_long, gs_min, 0, gsh, &buf);
  gs(maxp, gs_long, gs_max, 0, gsh, &buf);

  for (i = 0; i < size; i++)
    if (minp[i] != maxp[i]) {
      err = 1;
      break;
    }

  gs_free(gsh);
  gsh = gs_setup(vlp, size, &c, 0, gs_pairwise, 0);

  for (i = 0; i < size; i++)
    minp[i] = maxp[i] = vls[i];

  gs(minp, gs_long, gs_min, 0, gsh, &buf);
  gs(maxp, gs_long, gs_max, 0, gsh, &buf);

  for (i = 0; i < size; i++)
    if (minp[i] != maxp[i]) {
      err = 1;
      break;
    }

  int np;
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  if (np == 1) {
    for (i = 0; i < size; i++)
      if (vls[i] != vlp[i]) {
        err = 1;
        break;
      }
  }

  gs_free(gsh);
  comm_free(&c);
  buffer_free(&buf);

  if (minp != NULL)
    free(minp);
  if (maxp != NULL)
    free(maxp);
  if (vls != NULL)
    free(vls);

  return err;
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  parrsb_input *in = parrsb_parse_input(argc, argv);

  /* Read the geometry from the .re2 file */
  unsigned int nelt, nbcs;
  double *coord = NULL;
  long long *bcs = NULL;
  int nv;
  int err = parrsb_read_mesh(&nelt, &nv, NULL, &coord, &nbcs, &bcs, in->mesh,
                             MPI_COMM_WORLD, 1);
  parrsb_check_error(err, MPI_COMM_WORLD);

  /* Find connectivity */
  long long *vl = (long long *)calloc(nelt * nv, sizeof(long long));
  int ndim = nv == 8 ? 3 : 2;
  err |= parrsb_find_conn(vl, coord, nelt, ndim, bcs, nbcs, in->tol,
                          MPI_COMM_WORLD, 0);
  parrsb_check_error(err, MPI_COMM_WORLD);

  /* Test if the relevant env. variable is set */
  if (in->test == 1)
    err |= test_parcon(nelt, vl, in->mesh, MPI_COMM_WORLD);
  parrsb_check_error(err, MPI_COMM_WORLD);

  /* Write connectivity to .co2 file */
  if (in->dump == 1)
    err |= parrsb_dump_con(vl, nelt, nv, in->mesh, MPI_COMM_WORLD);
  parrsb_check_error(err, MPI_COMM_WORLD);

  /* Free resources */
  if (vl != NULL)
    free(vl);
  if (coord != NULL)
    free(coord);
  if (bcs != NULL)
    free(bcs);
  if (in != NULL)
    free(in);

  MPI_Finalize();

  return 0;
}

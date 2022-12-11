//=============================================================================
// Generate connectivity (.co2) from Nek5000 mesh (.re2) file.
//
#include "parRSB.h"

static int test_parcon(unsigned int neltp, long long *vlp, char *name,
                       MPI_Comm comm) {
  unsigned int nelt;
  int nv;
  long long *vls = NULL;
  int err = parrsb_read_mesh(&nelt, &nv, &vls, NULL, NULL, NULL, name, comm, 2);
  parrsb_check_error(err, comm);

  uint size = nelt * nv;
  slong *minp = tcalloc(slong, size);
  slong *maxp = tcalloc(slong, size);
  err = (minp == NULL || maxp == NULL);
  parrsb_check_error(err, comm);

  struct comm c;
  comm_init(&c, comm);
  struct gs_data *gsh = gs_setup(vls, size, &c, 0, gs_pairwise, 0);

  uint i;
  for (i = 0; i < size; i++)
    minp[i] = maxp[i] = vlp[i];

  buffer bfr;
  buffer_init(&bfr, 1024);
  gs(minp, gs_long, gs_min, 0, gsh, &bfr);
  gs(maxp, gs_long, gs_max, 0, gsh, &bfr);

  for (i = 0; err == 0 && i < size; i++) {
    if (minp[i] != maxp[i])
      err = 1;
  }

  gs_free(gsh);
  gsh = gs_setup(vlp, size, &c, 0, gs_pairwise, 0);

  for (i = 0; err == 0 && i < size; i++)
    minp[i] = maxp[i] = vls[i];

  gs(minp, gs_long, gs_min, 0, gsh, &bfr);
  gs(maxp, gs_long, gs_max, 0, gsh, &bfr);

  for (i = 0; err == 0 && i < size; i++) {
    if (minp[i] != maxp[i])
      err = 1;
  }

  if (c.np == 1) {
    for (i = 0; err == 0 && i < size; i++)
      if (vls[i] != vlp[i])
        err = 1;
  }

  gs_free(gsh), comm_free(&c), buffer_free(&bfr);
  free(minp), free(maxp), free(vls);

  return err;
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm world = MPI_COMM_WORLD;

  struct parrsb_cmd_opts *in = parrsb_parse_cmd_opts(argc, argv);
  int err = (in == NULL);
  parrsb_check_error(err, world);

  // Read the geometry from the .re2 file
  unsigned int nelt, nbcs, nv;
  double *coord = NULL;
  long long *bcs = NULL;
  err = parrsb_read_mesh(&nelt, &nv, NULL, &coord, &nbcs, &bcs, in->mesh, world,
                         1);
  parrsb_check_error(err, world);

  // Find connectivity
  long long *vl = (long long *)calloc(nelt * nv, sizeof(long long));
  err = (vl == NULL);
  parrsb_check_error(err, world);

  unsigned ndim = (nv == 8 ? 3 : 2);
  err = parrsb_conn_mesh(vl, coord, nelt, ndim, bcs, nbcs, in->tol, world);
  parrsb_check_error(err, world);

  // Write connectivity to .co2 file if dump is on
  if (in->dump) {
    err = parrsb_dump_con(in->mesh, nelt, nv, vl, world);
    parrsb_check_error(err, world);
  }

  // Turns on testing if test is on
  if (in->test) {
    err = test_parcon(nelt, vl, in->mesh, world);
    parrsb_check_error(err, world);
  }

  // Free resources
  free(vl), free(coord), free(bcs);
  parrsb_cmd_opts_free(in);
  MPI_Finalize();

  return 0;
}

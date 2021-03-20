#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include <gencon-impl.h>
#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  struct comm comm;
  comm_init(&comm, MPI_COMM_WORLD);
  int rank = comm.id, size = comm.np;

  if (argc != 2) {
    if (rank == 0)
      printf("Usage: ./%s <co2 file>\n", argv[0]);
    MPI_Finalize();
    exit(1);
  }

  Mesh mesh;
  read_co2_mesh(&mesh, argv[1], &comm);

  genmap_handle gh;
  genmap_init(&gh, MPI_COMM_WORLD);
  GenmapSetNLocalElements(gh, mesh->nelt);
  genmap_set_nvertices(gh, mesh->nVertex);

  /* Setup mesh */
  struct rsb_element *e = genmap_get_elements(gh);
  Element me = MeshGetElements(mesh);
  GenmapInt i, j;
  for (i = 0; i < mesh->nelt; i++)
    for (j = 0; j < mesh->nVertex; j++)
      e[i].vertices[j] = me[i].vertex[j].globalId;

  /* Setup CSR on fine level */
  genmap_comm c = genmap_global_comm(gh);
  struct array *entries = GenmapFindNeighbors(gh, c);
  csr_mat M;
  csr_mat_setup(entries, &c->gsc, &M);
  array_free(entries);
  free(entries);

  /* Setup MG levels */
  mgData d;
  mgSetup(c, M, &d);
  uint nlevels = d->nlevels;

  GenmapScalar *x = d->x, *y = d->y, *buf = d->buf;

  for (i = 0; i < d->level_off[nlevels]; i++)
    y[i] = 10., x[i] = 1.0;

  buffer bfr;
  buffer_init(&bfr, 1024);
  for (i = 0; i < nlevels; i++) {
    csr_mat Mi = d->levels[i]->M;
    csr_mat_gather(Mi, Mi->gsh, x + d->level_off[i], buf, &bfr);
    csr_mat_apply(y + d->level_off[i], Mi, buf);
  }
  buffer_free(&bfr);

  for (j = 0; j < d->level_off[nlevels]; j++)
    assert(fabs(y[j]) < GENMAP_TOL);

  mgFree(d);

  genmap_finalize(gh);
  mesh_free(mesh);

  comm_free(&comm);
  MPI_Finalize();

  return 0;
}

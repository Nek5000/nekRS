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
  struct array *entries = genmap_find_neighbors(gh, c);
  csr_mat M;
  csr_mat_setup(entries, &c->gsc, &M);
  array_free(entries);
  free(entries);

  slong out[2][1], bf[2][1], in = M->rn;
  comm_scan(out, &comm, gs_long, gs_add, &in, 1, bf);
  slong rg = out[1][0];

  /* Setup MG levels */
  mgData d;
  mgSetup(c, M, &d);

  uint nlevels = d->nlevels;
  mgLevel *l = d->levels;
  uint *lvl_off = d->level_off;
  GenmapScalar *x = d->x;

  for (i = lvl_off[0]; i < lvl_off[1]; i++)
    x[i] = 1.0;
  for (i = lvl_off[1]; i < lvl_off[nlevels]; i++)
    x[i] = 0.0;

  uint k;
  buffer buf;
  buffer_init(&buf, 1024);
  for (i = 0; i < nlevels - 1; i++) {
    gs(x + lvl_off[i], gs_double, gs_add, 1, l[i]->J, &buf);
  }

  if (lvl_off[nlevels - 1] < lvl_off[nlevels]) // rank with last 1-dof
    assert(fabs(x[lvl_off[nlevels] - 1] - rg) < GENMAP_TOL);

  for (i = nlevels - 2; i >= 0; i--) {
    gs(x + lvl_off[i], gs_double, gs_add, 0, l[i]->J, &buf);
  }

  for (i = lvl_off[0]; i < lvl_off[nlevels]; i++)
    assert(fabs(x[i] - rg) < GENMAP_TOL);

  buffer_free(&buf);

  mgFree(d);
  genmap_finalize(gh);
  mesh_free(mesh);

  comm_free(&comm);
  MPI_Finalize();

  return 0;
}

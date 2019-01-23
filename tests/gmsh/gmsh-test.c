#include <stdio.h>
#ifdef GENMAP_MPI
#include <mpi.h>
#endif

#include <genmap-impl.h>


int main(int argc, char **argv) {
#if defined(GENMAP_MPI)
  MPI_Init(&argc, &argv);
#else
  int MPI_COMM_WORLD = 0;
#endif

  GenmapHandle h;

  GenmapInit(&h, MPI_COMM_WORLD, "gmsh");
  GenmapRead(h, argv[1]);

  GenmapRSB(h);

  GenmapPartitionQuality(h);
  GenmapFinalize(h);

#if defined(GENMAP_MPI)
  MPI_Finalize();
#endif
}

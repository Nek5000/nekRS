/*
Parition mesh using Nek5000's vertex connectivity (con) file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "parRSB.h"
#include "conReader.h"

MPI_Comm comm = MPI_COMM_WORLD;

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int myid, np;
  MPI_Comm_rank(comm, &myid);
  MPI_Comm_size(comm, &np);

  int options[3];

  struct con c;
  int ierr = conRead(argv[1], &c, comm);
  if(ierr) goto quit;

  int nel_max = c.nelg / np + 1;
  long long *el = (long long*) malloc(nel_max * sizeof(long long));
  long long *vl = (long long*) malloc(nel_max * c.nv * sizeof(long long));

  int nelo = nel_max;
  options[0] = 1; // use custom options
  options[1] = 5; // debug level
  options[2] = 1; // print statistics
  ierr = parRSB_partMesh(el, vl, &nelo, c.el, c.vl, c.nel, c.nv, options,
                         comm);
  if(ierr) goto quit;


  conFree(&c);
  free(el);
  free(vl);

quit:
  MPI_Finalize();
  return ierr;
}

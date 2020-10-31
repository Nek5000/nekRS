/*
Parition mesh using Nek5000's vertex connectivity (con) file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "conReader.h"

#include "gslib.h"
#include "parRSB.h"
#include "quality.h"

#define MAXNV 8 /* maximum number of vertices per element */
typedef struct {
  int proc;
  long long id;
  long long vtx[MAXNV];
} elm_data;

int main(int argc, char *argv[]) {
  struct comm comm;
  struct crystal cr;
  struct array eList;
  elm_data *data;

  int ierr;
  int e, n, nel, nv;
  int options[3];

  MPI_Init(&argc, &argv);
  comm_init(&comm, MPI_COMM_WORLD);

  int color = MPI_UNDEFINED;

  if (argc != 3) {
    if(comm.id) printf("usage: ./example <#nread> <co2 file>\n");
    return EXIT_FAILURE;
  } 

  int nRead = atoi(argv[1]);
  char* conFile = argv[2];

  if (comm.id < nRead) color = 1;
  MPI_Comm comm_read;
  MPI_Comm_split(comm.c, color, 0, &comm_read);

  ierr = 0;
  struct con con = {};
  if (color != MPI_UNDEFINED) ierr = conRead(conFile, &con, comm_read);
  if(ierr) goto quit;

  comm_bcast(&comm, &con.nv, sizeof(int), 0);
  nv = con.nv;
  nel = con.nel;
  int *part = (int*) malloc(nel * sizeof(int));

  options[0] = 1; /* use custom options */
  options[1] = 3; /* debug level        */
  options[2] = 0; /* not used           */

  ierr = parRSB_partMesh(part, con.vl, nel, nv, options, comm.c);
  if(ierr) goto quit;

  /* redistribute data */
  array_init(elm_data, &eList, nel), eList.n = nel;
  for(data = eList.ptr, e = 0; e < nel; ++e) {
    data[e].proc = part[e];
    data[e].id   = con.el[e];
    for(n = 0; n < nv; ++n) {
      data[e].vtx[n] = con.vl[e * nv + n];
    }
  }
  free(part);
  conFree(&con);

  crystal_init(&cr, &comm);
  sarray_transfer(elm_data, &eList, proc, 0, &cr);
  crystal_free(&cr);
  nel = eList.n;

  long long *el = (long long*) malloc(nel * sizeof(long long));
  long long *vl = (long long*) malloc(nv * nel * sizeof(long long));
  for(data = eList.ptr, e = 0; e < nel; ++e) {
    el[e] = data[e].id;
    for(n = 0; n < nv; ++n) {
      vl[e * nv + n] = data[e].vtx[n];
    }
  }

  printPartStat(vl, nel, nv, comm.c);

  free(el);
  free(vl);
  array_free(&eList);
  comm_free(&comm);

quit:
  MPI_Finalize();
  return ierr;
}

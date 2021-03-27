#ifndef _CONREADER_H_
#define _CONREADER_H_

#include <math.h>

struct con {
  int nv;
  int nelg;
  int nel;
  long long *el, *vl;
};

void conFree(struct con *c) {
  free(c->el);
  free(c->vl);
}

int conRead(const char *fname, struct con *c, MPI_Comm comm) {
  int i, j, ierr;
  int myid, np;
  MPI_Comm_rank(comm, &myid);
  MPI_Comm_size(comm, &np);

  MPI_File fh;
  ierr = MPI_File_open(comm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  if (ierr != 0) {
    if (myid == 0)
      printf("ERROR: cannot open con file %s!\n", fname);
    return 1;
  }

  char hdr[132];
  MPI_File_read_all(fh, hdr, 132, MPI_BYTE, MPI_STATUS_IGNORE);
  char ver[6];
  int nelgt, nelgv, nv;
  sscanf(hdr, "%s %d %d %d", ver, &nelgt, &nelgv, &nv);

  float byte_test;
  MPI_File_read_all(fh, &byte_test, 4, MPI_BYTE, MPI_STATUS_IGNORE);
  if (fabs(byte_test - 6.543210) > 1e-7) {
    if (myid == 0)
      printf("ERROR byte_test failed! %f\n", byte_test);
    return 1;
  }

  int nelr = nelgt / np;
  for (i = 0; i < nelgt % np; ++i)
    if (np - i == myid)
      nelr++;

  int nelr_;
  MPI_Scan(&nelr, &nelr_, 1, MPI_INT, MPI_SUM, comm);
  long long off = sizeof(hdr) + sizeof(int);
  off += (long long)(nelr_ - nelr) * (nv + 1) * sizeof(int);
  MPI_File_set_view(fh, off, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

  int *buf = (int *)malloc(nelr * (nv + 1) * sizeof(int));
  ierr =
      MPI_File_read_all(fh, buf, (nv + 1) * nelr, MPI_INT, MPI_STATUS_IGNORE);
  if (ierr != 0) {
    if (myid == 0)
      printf("ERROR: failure while reading!\n");
    return 1;
  }

  MPI_File_close(&fh);

  c->el = (long long *)malloc(nelr * sizeof(long long));
  c->vl = (long long *)malloc(nv * nelr * sizeof(long long));
  for (i = 0; i < nelr; ++i) {
    c->el[i] = buf[i * (nv + 1)];
    for (j = 0; j < nv; ++j)
      c->vl[i * nv + j] = buf[i * (nv + 1) + 1 + j];
  }

  c->nv = nv;
  c->nelg = nelgv;
  c->nel = nelr;

  free(buf);
  return 0;
}

#endif

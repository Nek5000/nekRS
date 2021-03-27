#include <stdio.h>

#include <genmap-impl.h>

#define write_T(dest, val, T, nunits)                                          \
  do {                                                                         \
    T v = val;                                                                 \
    memcpy(dest, &(val), sizeof(T) * nunits);                                  \
    dest += sizeof(T) * nunits;                                                \
  } while (0)

int GenmapFiedlerDump(const char *fname, genmap_handle h, slong g_start,
                      struct comm *c) {
  uint rank = c->id;

  MPI_File file;
  int err = MPI_File_open(c->c, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &file);
  if (err != 0 && rank == 0)
    fprintf(stderr, "%s:%d Error opening file %s for writing.\n", __FILE__,
            __LINE__, fname);
  if (err != 0)
    return err;

  slong nelt = genmap_get_nel(h);
  slong out[2][1], buf[2][1];
  comm_scan(out, c, gs_long, gs_add, &nelt, 1, buf);
  slong start = out[0][0];
  slong nelgt = out[1][0];

  int ndim = (h->nv == 8) ? 3 : 2;
  uint write_size = ((ndim + 1) * sizeof(double) + 1 * sizeof(slong)) * nelt;
  if (rank == 0)
    write_size += sizeof(long) + sizeof(int); // for nelgt and ndim

  char *pbuf, *pbuf0;
  pbuf = pbuf0 = (char *)calloc(write_size, sizeof(char));
  if (rank == 0) {
    write_T(pbuf0, nelgt, long, 1);
    write_T(pbuf0, ndim, int, 1);
  }

  struct rsb_element *elm = genmap_get_elements(h);
  uint i;
  for (i = 0; i < nelt; i++) {
    slong id = g_start + i;
    write_T(pbuf0, id, long, 1);
    write_T(pbuf0, elm[i].coord[0], double, ndim);
    write_T(pbuf0, elm[i].fiedler, double, 1);
  }

  MPI_Status st;
  err = MPI_File_write_ordered(file, pbuf, write_size, MPI_BYTE, &st);
  if (err != 0 && rank == 0)
    fprintf(stderr, "%s:%d Error opening file %s for writing.\n", __FILE__,
            __LINE__, fname);

  err += MPI_File_close(&file);
  MPI_Barrier(c->c);

  free(pbuf);

  return err;
}

int GenmapVectorDump(const char *fname, GenmapScalar *y, uint size,
                     struct comm *c) {
  MPI_File file;
  int err = MPI_File_open(c->c, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &file);
  uint rank = c->id;
  if (err != 0 && rank == 0)
    fprintf(stderr, "%s:%d Error opening file %s for writing.\n", __FILE__,
            __LINE__, fname);
  if (err != 0)
    return err;

  slong nelt = size;
  slong out[2][1], buf[2][1];
  comm_scan(out, c, gs_long, gs_add, &nelt, 1, buf);
  slong start = out[0][0];
  slong nelgt = out[1][0];

  uint write_size = nelt * sizeof(double);
  if (rank == 0)
    write_size += sizeof(long); // nelgt

  char *pbuf, *pbuf0;
  pbuf = pbuf0 = (char *)calloc(write_size, sizeof(char));

  if (rank == 0)
    write_T(pbuf0, nelgt, long, 1);

  uint i;
  for (i = 0; i < nelt; i++)
    write_T(pbuf0, y[i], double, 1);

  MPI_Status st;
  err = MPI_File_write_ordered(file, pbuf, write_size, MPI_BYTE, &st);
  free(pbuf);

  if (err != 0 && rank == 0)
    fprintf(stderr, "%s:%d Error opening file %s for writing.\n", __FILE__,
            __LINE__, fname);

  err += MPI_File_close(&file);
  MPI_Barrier(c->c);

  return err;
}

int GenmapCentroidDump(const char *fname, genmap_handle h, sint g_id,
                       struct comm *c) {
  int rank = c->id;
  int size = c->np;

  MPI_File file;
  int err = MPI_File_open(c->c, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &file);
  if (err != 0 && rank == 0)
    fprintf(stderr, "%s:%d Error opening file %s for writing.\n", __FILE__,
            __LINE__, fname);
  if (err != 0)
    return err;

  slong nelt = genmap_get_nel(h);
  slong out[2][1], buf[2][1];
  comm_scan(out, c, gs_long, gs_add, &nelt, 1, buf);
  slong start = out[0][0];
  slong nelgt = out[1][0];

  int ndim = (h->nv == 8) ? 3 : 2;
  uint write_size = (ndim * sizeof(double) + sizeof(int)) * nelt;
  if (rank == 0)
    write_size += sizeof(slong) + sizeof(int); // for nelgt and ndim

  char *pbuf, *pbuf0;
  pbuf = pbuf0 = (char *)calloc(write_size, sizeof(char));
  if (rank == 0) {
    write_T(pbuf0, nelgt, slong, 1);
    write_T(pbuf0, ndim, int, 1);
  }

  struct rsb_element *elm = genmap_get_elements(h);
  uint i;
  for (i = 0; i < nelt; i++) {
    write_T(pbuf0, elm[i].coord[0], double, ndim);
    write_T(pbuf0, g_id, int, 1);
  }

  MPI_Status st;
  err = MPI_File_write_ordered(file, pbuf, write_size, MPI_BYTE, &st);
  if (err != 0 && rank == 0)
    fprintf(stderr, "%s:%d Error opening file %s for writing.\n", __FILE__,
            __LINE__, fname);

  err += MPI_File_close(&file);
  MPI_Barrier(c->c);

  free(pbuf);

  return err;
}

int GenmapElementIdDump(const char *fname, genmap_handle h, struct comm *c) {
  int rank = c->id;
  int size = c->np;

  MPI_File file;
  int err = MPI_File_open(c->c, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &file);
  if (err != 0 && rank == 0)
    fprintf(stderr, "%s:%d Error opening file %s for writing.\n", __FILE__,
            __LINE__, fname);
  if (err != 0)
    return err;

  slong nelt = genmap_get_nel(h);
  slong out[2][1], buf[2][1];
  comm_scan(out, c, gs_long, gs_add, &nelt, 1, buf);
  slong start = out[0][0];
  slong nelgt = out[1][0];

  int ndim = (h->nv == 8) ? 3 : 2;
  uint write_size = sizeof(GenmapLong) * nelt;
  if (rank == 0)
    write_size += sizeof(slong) + sizeof(int); // for nelgt and ndim

  char *pbuf, *pbuf0;
  pbuf = pbuf0 = (char *)calloc(write_size, sizeof(char));
  if (rank == 0) {
    write_T(pbuf0, nelgt, slong, 1);
    write_T(pbuf0, ndim, int, 1);
  }

  struct rsb_element *elm = genmap_get_elements(h);
  uint i;
  for (i = 0; i < nelt; i++)
    write_T(pbuf0, elm[i].globalId, GenmapLong, 1);

  MPI_Status st;
  err = MPI_File_write_ordered(file, pbuf, write_size, MPI_BYTE, &st);
  if (err != 0 && rank == 0)
    fprintf(stderr, "%s:%d Error opening file %s for writing.\n", __FILE__,
            __LINE__, fname);

  err += MPI_File_close(&file);
  MPI_Barrier(c->c);

  free(pbuf);

  return err;
}
#undef write_T

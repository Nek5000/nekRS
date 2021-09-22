#include <math.h>
#include <stdio.h>

#include <gencon-impl.h>
#include <genmap.h>
#include <parRSB.h>

#define GC_RE2_HEADER_LEN 80
#define GC_CO2_HEADER_LEN 132

static int readRe2Header(unsigned int *nelt_, int *nv_, ulong *nelgt_,
                         ulong *nelgv_, MPI_File file, struct comm *c) {
  int rank = c->id;
  int size = c->np;
  MPI_Comm comm = c->c;

  int err = 0;
  char *buf = (char *)calloc(GC_RE2_HEADER_LEN + 1, sizeof(char));
  MPI_Status st;
  err = MPI_File_read_all(file, buf, GC_RE2_HEADER_LEN, MPI_BYTE, &st);

  int nelgt, nelgv, ndim;
  char version[6];
  sscanf(buf, "%5s %d %d %d", version, &nelgt, &ndim, &nelgv);

  /* TODO: Assert version */
  int nelt = nelgt / size;
  int nrem = nelgt - nelt * size;
  nelt += (rank > (size - 1 - nrem) ? 1 : 0);

  float byte_test;
  MPI_File_read_all(file, &byte_test, 4, MPI_BYTE, &st);
  if (fabs(byte_test - 6.543210) > 1e-7) {
    if (rank == 0)
      printf("ERROR byte_test failed! %f\n", byte_test);
    err = 1;
  }

  *nelt_ = nelt;
  *nv_ = (ndim == 2) ? 4 : 8;
  *nelgt_ = nelgt;
  *nelgv_ = nelgv;

  free(buf);

  return err;
}

static int readRe2Coordinates(double **coord_, unsigned int nelt, int nv,
                              MPI_File file, struct comm *c) {
  uint rank = c->id;
  uint size = c->np;
  MPI_Comm comm = c->c;

  slong out[2][1], bfr[2][1];
  slong in = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bfr);
  slong start = out[0][0];

  int ndim = (nv == 4) ? 2 : 3;
  int elemDataSize = nv * ndim * sizeof(double) + sizeof(double);
  int header_size = GC_RE2_HEADER_LEN + sizeof(float);

  /* calculate read size for element data on each MPI rank */
  int read_size = nelt * elemDataSize;
  if (rank == 0)
    read_size += header_size;

  int err = 0;
  char *buf = (char *)calloc(read_size, sizeof(char));
  char *buf0 = buf;

  MPI_Status st;
  err = MPI_File_read_ordered(file, buf, read_size, MPI_BYTE, &st);

  if (rank == 0)
    buf0 += header_size;

  /* Allocate coord array */
  size_t coord_size = nelt;
  coord_size = coord_size * nv * ndim;
  double *coord = *coord_ = tcalloc(double, coord_size);

  /* Read elements for each rank */
  double x[GC_MAX_VERTICES], y[GC_MAX_VERTICES], z[GC_MAX_VERTICES];
  int i, j, k;
  if (ndim == 3) {
    for (i = 0; i < nelt; i++) {
      // skip group id
      buf0 += sizeof(double);
      READ_T(x, buf0, double, nv);
      buf0 += sizeof(double) * nv;
      READ_T(y, buf0, double, nv);
      buf0 += sizeof(double) * nv;
      READ_T(z, buf0, double, nv);
      buf0 += sizeof(double) * nv;

      for (k = 0; k < nv; k++) {
        coord[i * nv * ndim + k * ndim + 0] = x[k];
        coord[i * nv * ndim + k * ndim + 1] = y[k];
        coord[i * nv * ndim + k * ndim + 2] = z[k];
      }
    }
  } else if (ndim == 2) {
    for (i = 0; i < nelt; i++) {
      // skip group id
      buf0 += sizeof(double);
      READ_T(x, buf0, double, nv);
      buf0 += sizeof(double) * nv;
      READ_T(y, buf0, double, nv);
      buf0 += sizeof(double) * nv;

      for (k = 0; k < nv; k++) {
        coord[i * nv * ndim + k * ndim + 0] = x[k];
        coord[i * nv * ndim + k * ndim + 1] = y[k];
      }
    }
  }

  free(buf);

  return 0;
}

static int readRe2Boundaries(unsigned int *nbcs_, long long **bcs_,
                             unsigned int nelt, int nv, ulong nelgt,
                             MPI_File file, struct comm *c) {
  uint rank = c->id;
  uint size = c->np;
  MPI_Comm comm = c->c;

  int ndim = (nv == 4) ? 2 : 3;
  int elemDataSize = nv * ndim * sizeof(double) + sizeof(double);
  int header_size = GC_RE2_HEADER_LEN + sizeof(float);

  /* Calculate offset for the curve side data */
  MPI_Offset curveOffset = header_size + nelgt * elemDataSize;

  MPI_Status st;
  char bufL[8];
  if (rank == 0)
    MPI_File_read_at(file, curveOffset, bufL, sizeof(long), MPI_BYTE, &st);
  MPI_Bcast(bufL, sizeof(long), MPI_BYTE, 0, comm);

  double ncurvesD;
  READ_T(&ncurvesD, bufL, long, 1);
  long ncurves = ncurvesD;

  /* calculate offset for boundary conditions data */
  MPI_Offset boundaryOffset =
      curveOffset + sizeof(long) + sizeof(long) * 8 * ncurves;
  if (rank == 0)
    MPI_File_read_at(file, boundaryOffset, bufL, sizeof(long), MPI_BYTE, &st);
  MPI_Bcast(bufL, sizeof(long), MPI_BYTE, 0, comm);

  double nbcsD;
  READ_T(&nbcsD, bufL, long, 1);
  long nbcs = nbcsD;

  int nbcsLocal = nbcs / size;
  int nrem = nbcs - nbcsLocal * size;
  nbcsLocal += (rank > (size - 1 - nrem) ? 1 : 0);

  slong out[2][1], bfr[2][1], in = nbcsLocal;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bfr);
  slong start = out[0][0];

  int offset = boundaryOffset + sizeof(long) + start * 8 * sizeof(long);
  int read_size = nbcsLocal * sizeof(long) * 8;
  char *buf = calloc(read_size, sizeof(char)), *buf0 = buf;
  MPI_File_read_at_all(file, offset, buf, read_size, MPI_BYTE, &st);

  struct array barray;
  array_init(struct Boundary_private, &barray, 10);

  double tmp[5];
  char cbc[4];
  struct Boundary_private boundary;
  sint i;
  for (i = 0; i < nbcsLocal; i++) {
    READ_T(tmp, buf0, long, 1);
    buf0 += sizeof(long);
    boundary.elementId = (long)tmp[0];

    READ_T(tmp, buf0, long, 1);
    buf0 += sizeof(long);
    boundary.faceId = (long)tmp[0];

    READ_T(tmp, buf0, long, 5);
    buf0 += 5 * sizeof(long);
    READ_T(cbc, buf0, char, 3);
    buf0 += sizeof(long);
    cbc[3] = '\0';

    if (strcmp(cbc, GC_PERIODIC) == 0) {
      boundary.bc[0] = (long)tmp[0];
      boundary.bc[1] = (long)tmp[1];
      array_cat(struct Boundary_private, &barray, &boundary, 1);
    }
  }

  nbcs = *nbcs_ = barray.n;
  struct Boundary_private *ptr = barray.ptr;
  long long *bcs = *bcs_ = tcalloc(long long, 4 * nbcs);

  for (i = 0; i < nbcs; i++) {
    bcs[4 * i + 0] = ptr[i].elementId;
    bcs[4 * i + 1] = ptr[i].faceId;
    bcs[4 * i + 2] = ptr[i].bc[0];
    bcs[4 * i + 3] = ptr[i].bc[1];
  }

  array_free(&barray);
  free(buf);

  return 0;
}

static int read_geometry(unsigned int *nelt, int *nv, double **coord,
                         unsigned int *nbcs, long long **bcs, char *fname,
                         struct comm *c) {
  int errs = 0;
  uint rank = c->id;
  uint size = c->np;
  MPI_Comm comm = c->c;

  MPI_File file;
  int err = MPI_File_open(comm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
  if (err != 0) {
    if (rank == 0)
      printf("%s:%d Error opening file: %s\n", __FILE__, __LINE__, fname);
    return 1;
  }

  ulong nelgt, nelgv;
  readRe2Header(nelt, nv, &nelgt, &nelgv, file, c);
  readRe2Coordinates(coord, *nelt, *nv, file, c);
  readRe2Boundaries(nbcs, bcs, *nelt, *nv, nelgt, file, c);

  err = MPI_File_close(&file);
  if (err)
    errs++;

  return errs;
}

static int read_connectivity(unsigned int *nelt_, int *nv_, long long **vl_,
                             char *fname, struct comm *c) {
  comm_ext comm = c->c;
  int rank = c->id;
  int size = c->np;

  MPI_File file;
  int err = MPI_File_open(comm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
  if (err != 0) {
    if (rank == 0)
      printf("Error opening %s for reading.\n", fname);
  }

  char *buf = (char *)calloc(GC_CO2_HEADER_LEN + 1, sizeof(char));
  MPI_Status st;
  err = MPI_File_read_all(file, buf, GC_CO2_HEADER_LEN, MPI_BYTE, &st);

  int nelgt, nelgv, nv;
  char version[6];
  sscanf(buf, "%5s%12d%12d%12d", version, &nelgt, &nelgv, &nv);

  /* TODO: Assert version */
  int nelt = nelgt / size;
  int nrem = nelgt - nelt * size;
  nelt += (rank > (size - 1 - nrem) ? 1 : 0);

  if (*nv_ != 0) {
    assert(*nv_ == nv);
    assert(*nelt_ == nelt);
  } else {
    *nv_ = nv;
    *nelt_ = nelt;
  }

  float byte_test;
  MPI_File_read_all(file, &byte_test, 4, MPI_BYTE, &st);
  if (fabs(byte_test - 6.543210) > 1e-7) {
    if (rank == 0)
      printf("ERROR byte_test failed! %f\n", byte_test);
    return 1;
  }

  slong out[2][1], bfr[2][1], in[1];
  in[0] = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bfr);
  slong start = out[0][0];

  int read_size = nelt * (nv + 1) * sizeof(int);
  int header_size = GC_CO2_HEADER_LEN + sizeof(float);
  if (rank == 0)
    read_size += header_size;

  buf = (char *)realloc(buf, read_size * sizeof(char));
  err = MPI_File_read_ordered(file, buf, read_size, MPI_BYTE, &st);
  err = MPI_File_close(&file);

  char *buf0 = buf;
  if (rank == 0)
    buf0 += header_size;

  long long *vl = *vl_ = tcalloc(long long, nv *nelt);
  int i, j, tmp1, tmp2;
  for (i = 0; i < nelt; i++) {
    READ_T(&tmp1, buf0, int, 1);
    buf0 += sizeof(int);
    for (j = 0; j < nv; j++) {
      READ_T(&tmp2, buf0, int, 1);
      buf0 += sizeof(int);
      vl[i * nv + j] = tmp2;
    }
  }

  free(buf);

  return 0;
}

int parrsb_read_mesh(unsigned int *nel, int *nv, long long **vl, double **coord,
                     unsigned int *nbcs, long long **bcs, char *name,
                     MPI_Comm comm, int read) {
  struct comm c;
  comm_init(&c, comm);

  /* Set nv to 0 so we know if .re2 is read before .co2 */
  *nv = 0;

  /* Read geometry from .re2 file */
  if (read & 1) {
    char geom_name[BUFSIZ];
    strncpy(geom_name, name, BUFSIZ);
    strncat(geom_name, ".re2", 5);
    read_geometry(nel, nv, coord, nbcs, bcs, geom_name, &c);
  }

  /* Read connectivity from .co2 file */
  if (read & 2) {
    char conn_name[BUFSIZ];
    strncpy(conn_name, name, BUFSIZ);
    strncat(conn_name, ".co2", 5);
    read_connectivity(nel, nv, vl, conn_name, &c);
  }

  comm_free(&c);

  return 0;
}

#define WRITE_INT(dest, val)                                                   \
  do {                                                                         \
    memcpy(dest, &(val), sizeof(int));                                         \
  } while (0)

int parrsb_dump_con(char *name, unsigned int nelt, int nv, long long *vl,
                    MPI_Comm comm) {
  const char version[6] = "#v001";
  const float test = 6.54321;

  struct comm c;
  comm_init(&c, comm);
  uint id = c.id;

  char co2_name[BUFSIZ];
  strncpy(co2_name, name, BUFSIZ);
  strncat(co2_name, ".co2", 5);

  MPI_File file;
  int err = MPI_File_open(comm, co2_name, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &file);
  if (err) {
    if (id == 0) {
      fprintf(stderr, "%s:%d Error opening file: %s for writing.\n", __FILE__,
              __LINE__, co2_name);
      fflush(stdout);
    }
    return err;
  }

  slong out[2][1], bfr[2][1];
  slong in = nelt;
  comm_scan(out, &c, gs_long, gs_add, &in, 1, bfr);
  slong start = out[0][0];
  slong nelgt = out[1][0];
  slong nelgv = nelgt;

  int write_size = nelt * (nv + 1) * sizeof(int);
  int header_size = GC_CO2_HEADER_LEN + sizeof(float);
  if (id == 0)
    write_size += header_size;

  char *buf = (char *)calloc(write_size, sizeof(char));
  char *buf0 = buf;

  if (id == 0) {
    sprintf(buf0, "%5s%12d%12d%12d", version, (int)nelgt, (int)nelgv, nv);
    memset(buf0 + strlen(buf0), ' ', GC_CO2_HEADER_LEN - strlen(buf0));
    buf0[GC_CO2_HEADER_LEN] = '\0';
    buf0 += GC_CO2_HEADER_LEN;
    memcpy(buf0, &test, sizeof(float));
    buf0 += sizeof(float);
  }

  int i, j, temp;
  for (i = 0; i < nelt; i++) {
    temp = start + i + 1;
    WRITE_INT(buf0, temp);
    buf0 += sizeof(int);
    for (j = 0; j < nv; j++) {
      temp = vl[i * nv + j];
      WRITE_INT(buf0, temp);
      buf0 += sizeof(int);
    }
  }

  MPI_Status st;
  err |= MPI_File_write_ordered(file, buf, write_size, MPI_BYTE, &st);
  err |= MPI_File_close(&file);

  genmap_barrier(&c);
  comm_free(&c);

  free(buf);

  return err;
}

#undef WRITE_INT

#define HEADER_LEN 132

int parrsb_dump_map(char *name, unsigned int nelt, int nv, long long *vtx,
                    int *pmap, MPI_Comm comm) {
  char version[6] = "#v001";
  float test = 6.54321;

  int nelgt = nelt;
  MPI_Allreduce(&nelt, &nelgt, 1, MPI_INT, MPI_SUM, comm);

  const int npts = nelgt * nv;
  const int depth = (int)log2(1.0 * nelgt);
  const int d2 = (int)(pow(2, depth) + 0.5);
  int nactive = nelgt;
  int nrnk = nelgt;
  int noutflow = 0;

  char header[HEADER_LEN];
  header[HEADER_LEN] = '\0';
  sprintf(header, "%5s%12d%12d%12d%12d%12d%12d%12d", version, nelgt, nactive,
          depth, d2, npts, nrnk, noutflow);
  memset(header + strlen(header), ' ', HEADER_LEN - strlen(header));
  header[HEADER_LEN] = '\0';

  MPI_Info infoIn;
  MPI_Info_create(&infoIn);
  MPI_Info_set(infoIn, "access_style", "write_once,random");

  int errs = 0;
  MPI_File file;
  int err = MPI_File_open(comm, name, MPI_MODE_WRONLY | MPI_MODE_CREATE, infoIn,
                          &file);
  if (err) {
    errs++;
    MPI_Abort(comm, 911);
  }

  int rank;
  MPI_Comm_rank(comm, &rank);

  int writeSize = 0;
  if (rank == 0)
    writeSize = HEADER_LEN * sizeof(char) + sizeof(float);
  writeSize += (nv + 1) * nelt * sizeof(int);

  char *buf = (char *)malloc(writeSize), *buf0 = buf;

  if (rank == 0) {
    memcpy(buf0, header, HEADER_LEN * sizeof(char)),
        buf0 += HEADER_LEN * sizeof(char);
    memcpy(buf0, &test, sizeof(float)), buf0 += sizeof(float);
  }

  int ivtx[8];
  int i, j;
  for (i = 0; i < nelt; i++) {
    memcpy(buf0, &pmap[i], sizeof(int));
    buf0 += sizeof(int);

    for (j = 0; j < nv; j++)
      ivtx[j] = vtx[i * nv + j];
    memcpy(buf0, ivtx, sizeof(int) * nv);
    buf0 += nv * sizeof(int);
  }

  MPI_Status status;
  err = MPI_File_write_ordered(file, buf, writeSize, MPI_BYTE, &status);
  if (err)
    errs++;

  err = MPI_File_close(&file);
  if (err)
    errs++;

  free(buf);

  return errs;
}

#undef HEADER_LEN
#undef GC_RE2_HEADER_LEN
#undef GC_CO2_HEADER_LEN

#define WRITE_T(dest, val, T, nunits)                                          \
  do {                                                                         \
    memcpy(dest, val, sizeof(T) * nunits);                                     \
    dest += sizeof(T) * nunits;                                                \
  } while (0)

int parrsb_dump_part(char *name, unsigned int nel, int nv, double *coord,
                     int gid, MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  int rank = c.id, size = c.np;

  MPI_File file;
  int err = MPI_File_open(comm, name, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &file);
  parrsb_check_error(err, comm);

  slong nelt = nel;
  slong out[2][1], buf[2][1];
  comm_scan(out, &c, gs_long, gs_add, &nelt, 1, buf);
  slong start = out[0][0];
  slong nelgt = out[1][0];

  int ndim = nv == 8 ? 3 : 2;
  uint write_size = (ndim * sizeof(double) + sizeof(int)) * nelt;
  if (rank == 0)
    write_size += sizeof(slong) + sizeof(int); // for nelgt and ndim

  char *pbuf, *pbuf0;
  pbuf = pbuf0 = (char *)tcalloc(char, write_size);
  if (rank == 0) {
    WRITE_T(pbuf0, &nelgt, slong, 1);
    WRITE_T(pbuf0, &ndim, int, 1);
  }

  uint i, j, k;
  double tcoord[3];
  for (i = 0; i < nelt; i++) {
    tcoord[0] = tcoord[1] = tcoord[2] = 0.0;
    for (j = 0; j < nv; j++)
      for (k = 0; k < ndim; k++)
        tcoord[k] += coord[i * nv * ndim + j * ndim + k];
    tcoord[0] /= nv;
    tcoord[1] /= nv;
    tcoord[2] /= nv;
    WRITE_T(pbuf0, tcoord, double, ndim);
    WRITE_T(pbuf0, &gid, int, 1);
  }

  MPI_Status st;
  err = MPI_File_write_ordered(file, pbuf, write_size, MPI_BYTE, &st);
  parrsb_check_error(err, comm);

  err += MPI_File_close(&file);
  parrsb_check_error(err, comm);

  free(pbuf);

  return err;
}

#undef WRITE_T

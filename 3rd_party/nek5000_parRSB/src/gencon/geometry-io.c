#include <math.h>
#include <string.h>

#include <gencon-impl.h>
#include <genmap-impl.h>

int PRE_TO_SYM_VERTEX[GC_MAX_VERTICES] = {0, 1, 3, 2, 4, 5, 7, 6};
int PRE_TO_SYM_FACE[GC_MAX_FACES] = {2, 1, 3, 0, 4, 5};

int transferBoundaryFaces(Mesh mesh, struct comm *c) {
  uint size = c->np;

  struct array *boundary = &mesh->boundary;
  BoundaryFace ptr = boundary->ptr;
  int nFaces = boundary->n;

  slong nelgt = mesh->nelgt;
  sint nelt = nelgt / size;
  sint nrem = nelgt - nelt * size;
  slong N = (size - nrem) * nelt;

  sint i;
  slong eid;
  for (i = 0; i < nFaces; i++) {
    eid = ptr[i].elementId;
    if (eid < N)
      ptr[i].proc = eid / nelt;
    else
      ptr[i].proc = (eid - N) / (nelt + 1) + size - nrem;
  }

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct Boundary_private, boundary, proc, 1, &cr);
  crystal_free(&cr);

  return 0;
}

int readRe2Header(Mesh *mesh_, MPI_File file, struct comm *c) {
  int rank = c->id;
  int size = c->np;
  MPI_Comm comm = c->c;

  char *buf = (char *)calloc(GC_RE2_HEADER_LEN + 1, sizeof(char));
  MPI_Status st;
  int err = MPI_File_read_all(file, buf, GC_RE2_HEADER_LEN, MPI_BYTE, &st);

  int nelgt, nelgv, nDim;
  char version[6];
  sscanf(buf, "%5s %d %d %d", version, &nelgt, &nDim, &nelgv);

  // TODO: Assert version
  int nVertex = (nDim == 2) ? 4 : 8;
  int nelt = nelgt / size;
  int nrem = nelgt - nelt * size;
  nelt += (rank > (size - 1 - nrem) ? 1 : 0);

  float byte_test;
  MPI_File_read_all(file, &byte_test, 4, MPI_BYTE, &st);
  if (fabs(byte_test - 6.543210) > 1e-7) {
    if (rank == 0)
      printf("ERROR byte_test failed! %f\n", byte_test);
    return 1;
  }

  mesh_init(mesh_, nelt, nDim);
  Mesh mesh = *mesh_;
  mesh->nelgt = nelgt;
  mesh->nelgv = nelgv;
  mesh->nelt = nelt;
  mesh->nDim = nDim;

#if defined(GENMAP_DEBUG)
  printf("re2 %d: %s %d %d %d %d\n", rank, version, nelgt, nelgv, nelt, nDim);
#endif

  free(buf);

  return 0;
}

int readRe2Coordinates(Mesh mesh, MPI_File file, struct comm *c) {
  uint rank = c->id;
  uint size = c->np;
  MPI_Comm comm = c->c;

  int nelt = mesh->nelt;
  int nelgt = mesh->nelgt;
  int nDim = mesh->nDim;
  int nVertex = (nDim == 2) ? 4 : 8;

  slong out[2][1], bfr[2][1];
  slong in = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bfr);
  slong start = out[0][0];

  int elemDataSize = nVertex * nDim * sizeof(double) + sizeof(double);
  int header_size = GC_RE2_HEADER_LEN + sizeof(float);

  /* calculate read size for element data on each MPI rank */
  int read_size = nelt * elemDataSize;
  if (rank == 0)
    read_size += header_size;

  char *buf = (char *)calloc(read_size, sizeof(char));
  char *buf0 = buf;
  MPI_Status st;
  int err = MPI_File_read_ordered(file, buf, read_size, MPI_BYTE, &st);
  if (err)
    return 1;

  if (rank == 0)
    buf0 += header_size;

  /* initialize array */
  uint nUnits = nelt * nVertex;
  array_init(struct Point_private, &mesh->elements, nUnits);
  Point ptr = mesh->elements.ptr;

  /* read elements for each rank */
  double x[GC_MAX_VERTICES], y[GC_MAX_VERTICES], z[GC_MAX_VERTICES];
  int i, j, k;
  for (i = 0; i < nelt; i++) {
    // skip group id
    buf0 += sizeof(double);
    READ_T(x, buf0, double, nVertex);
    buf0 += sizeof(double) * nVertex;
    READ_T(y, buf0, double, nVertex);
    buf0 += sizeof(double) * nVertex;
    if (nDim == 3) {
      READ_T(z, buf0, double, nVertex);
      buf0 += sizeof(double) * nVertex;
    }

    for (k = 0; k < nVertex; k++) {
      j = PRE_TO_SYM_VERTEX[k];
      ptr->x[0] = x[j], ptr->x[1] = y[j];
      if (nDim == 3)
        ptr->x[2] = z[j];
      ptr->elementId = start + i;
      ptr->sequenceId = nVertex * (start + i) + k;
      ptr->origin = rank;
      ptr++;
    }
  }
  mesh->elements.n = nUnits;

#if defined(GENMAP_DEBUG)
  printf("io: rank=%d npts=%u\n", rank, nUnits);
#endif

  free(buf);
  return 0;
}

int readRe2Boundaries(Mesh mesh, MPI_File file, struct comm *c) {
  uint rank = c->id;
  uint size = c->np;
  MPI_Comm comm = c->c;

  int nelt = mesh->nelt;
  int nelgt = mesh->nelgt;
  int nDim = mesh->nDim;
  int nVertex = mesh->nVertex;

  int elemDataSize = nVertex * nDim * sizeof(double) + sizeof(double);
  int header_size = GC_RE2_HEADER_LEN + sizeof(float);

  MPI_Status st;
  char bufL[8];

  /* calculate offset for the curve side data */
  MPI_Offset curveOffset = header_size + nelgt * elemDataSize;
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

  int nbcsLocal = nbcs / size, nrem = nbcs - nbcsLocal * size;
  nbcsLocal += (rank >= (size - nrem) ? 1 : 0);

  slong out[2][1], bfr[2][1], in = nbcsLocal;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bfr);
  slong start = out[0][0];

  int offset = boundaryOffset + sizeof(long) + start * 8 * sizeof(long);
  int read_size = nbcsLocal * sizeof(long) * 8;
  char *buf = calloc(read_size, sizeof(char)), *buf0 = buf;
  MPI_File_read_at_all(file, offset, buf, read_size, MPI_BYTE, &st);

  double tmp[5];
  char cbc[4];
  struct Boundary_private boundary;
  sint i;
  for (i = 0; i < nbcsLocal; i++) {
    READ_T(tmp, buf0, long, 1);
    buf0 += sizeof(long);
    boundary.elementId = tmp[0] - 1;

    READ_T(tmp, buf0, long, 1);
    buf0 += sizeof(long);
    boundary.faceId = PRE_TO_SYM_FACE[(long)tmp[0] - 1];

    READ_T(tmp, buf0, long, 5);
    buf0 += 5 * sizeof(long);
    READ_T(cbc, buf0, char, 3);
    buf0 += sizeof(long);
    cbc[3] = '\0';

    if (strcmp(cbc, GC_PERIODIC) == 0) {
      boundary.bc[0] = (long)tmp[0] - 1;
      boundary.bc[1] = PRE_TO_SYM_FACE[(long)tmp[1] - 1];
      array_cat(struct Boundary_private, &mesh->boundary, &boundary, 1);
    }
  }
  free(buf);
}

int read_geometry(Mesh *mesh_, char *fname, struct comm *c) {
  int nelt, nDim, nVertex;
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

  readRe2Header(mesh_, file, c);
  Mesh mesh = *mesh_;
  readRe2Coordinates(mesh, file, c);
  readRe2Boundaries(mesh, file, c);
  transferBoundaryFaces(mesh, c);

  err = MPI_File_close(&file);
  if (err)
    errs++;

  MPI_Barrier(comm);

  return errs;
}

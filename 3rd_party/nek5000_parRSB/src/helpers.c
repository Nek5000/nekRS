#include "parrsb-impl.h"
#include <getopt.h>

#include <sys/resource.h>

#if defined __GLIBC__
#include <execinfo.h>

// Obtain a backtrace and print it to stdout.
void parrsb_print_stack(void) {
  void *bt[50];
  int bt_size = backtrace(bt, 50);
  char **symbols = backtrace_symbols(bt, bt_size);
  printf("backtrace(): obtained %d stack frames.\n", bt_size);
  for (unsigned i = 0; i < bt_size; i++)
    printf("%s\n", symbols[i]);
  free(symbols);
}
#else
void parrsb_print_stack(){};
#endif // defined __GLIBC__

double parrsb_get_max_rss() {
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, &r_usage);
#if defined(__APPLE__) && defined(__MACH__)
  return (double)r_usage.ru_maxrss;
#else
  return (double)(r_usage.ru_maxrss * 1024L);
#endif
}
int log2ll(long long n) {
  int k = 0;
  while (n > 1)
    n /= 2, k++;

  return k;
}

int parrsb_dist_mesh(unsigned int *nelt_, long long **vl_, double **coord_,
                     int *part, int nv, MPI_Comm comm) {
  typedef struct {
    int proc;
    long long vtx[MAXNV];
    double coord[MAXNV * MAXDIM];
  } elem_data;

  uint nelt = *nelt_;
  struct array elements;
  array_init(elem_data, &elements, nelt);

  elem_data data;
  long long *vl = *vl_;
  uint e, n;
  for (e = 0; e < nelt; ++e) {
    data.proc = part[e];
    for (n = 0; n < nv; ++n)
      data.vtx[n] = vl[e * nv + n];
    array_cat(elem_data, &elements, &data, 1);
  }
  assert(elements.n == nelt);

  int ndim = (nv == 8) ? 3 : 2;
  elem_data *ed = elements.ptr;
  double *coord = (coord_ == NULL ? NULL : *coord_);
  if (coord != NULL) {
    for (e = 0; e < nelt; e++)
      for (n = 0; n < ndim * nv; n++)
        ed[e].coord[n] = coord[e * ndim * nv + n];
  }

  struct comm c;
  comm_init(&c, comm);

  struct crystal cr;
  crystal_init(&cr, &c);

  sarray_transfer(elem_data, &elements, proc, 0, &cr);
  *nelt_ = nelt = elements.n;
  ed = elements.ptr;

  vl = *vl_ = (long long *)realloc(*vl_, nv * nelt * sizeof(long long));
  for (e = 0; e < nelt; ++e)
    for (n = 0; n < nv; ++n)
      vl[e * nv + n] = ed[e].vtx[n];

  if (coord != NULL) {
    coord = *coord_ =
        (double *)realloc(*coord_, ndim * nv * nelt * sizeof(double));
    for (e = 0; e < nelt; ++e) {
      for (n = 0; n < ndim * nv; ++n)
        coord[e * ndim * nv + n] = ed[e].coord[n];
    }
  }

  crystal_free(&cr);
  comm_free(&c);
  array_free(&elements);

  return 0;
}

int parrsb_setup_mesh(unsigned *nelt, unsigned *nv, long long **vl,
                      double **coord, struct parrsb_cmd_opts *in,
                      MPI_Comm comm) {
  int id, err;
  MPI_Comm_rank(comm, &id);

  // Read the geometry from the .re2 file
  unsigned int nbcs;
  long long *bcs = NULL;
  err = parrsb_read_mesh(nelt, nv, NULL, coord, &nbcs, &bcs, in->mesh, comm, 1);
  parrsb_check_error(err, comm);

  // Find connectivity
  *vl = (long long *)calloc(*nelt * *nv, sizeof(long long));
  err = (*vl == NULL);
  parrsb_check_error(err, comm);

  int ndim = (*nv == 8 ? 3 : 2);
  err = parrsb_conn_mesh(*vl, *coord, *nelt, ndim, bcs, nbcs, in->tol, comm);
  parrsb_check_error(err, comm);

  // Partition the mesh
  int *part = (int *)calloc(*nelt, sizeof(int));
  err = (part == NULL);
  parrsb_check_error(err, comm);

  parrsb_options paropt = parrsb_default_options;
  err = parrsb_part_mesh(part, NULL, *vl, *coord, *nelt, *nv, paropt, comm);
  parrsb_check_error(err, comm);

  // Redistribute data based on identified partitions
  err = parrsb_dist_mesh(nelt, vl, coord, part, *nv, comm);
  parrsb_check_error(err, comm);

  free(part), free(bcs);

  return 0;
}

void parrsb_get_part_stat(int *nc, int *ns, int *nss, int *nel, long long *vtx,
                          int nelt, int nv, MPI_Comm ce) {
  struct comm comm;
  comm_init(&comm, ce);

  int np = comm.np;
  int id = comm.id;

  if (np == 1)
    return;

  int Npts = nelt * nv;
  int i;
  slong *data = (slong *)malloc((Npts + 1) * sizeof(slong));
  for (i = 0; i < Npts; i++)
    data[i] = vtx[i];
  struct gs_data *gsh = gs_setup(data, Npts, &comm, 0, gs_pairwise, 0);

  int Nmsg;
  pw_data_nmsg(gsh, &Nmsg);

  int *Ncomm = (int *)malloc((Nmsg + 1) * sizeof(int));
  pw_data_size(gsh, Ncomm);

  gs_free(gsh);
  free(data);

  int nelMin, nelMax, nelSum;
  int ncMin, ncMax, ncSum;
  int nsMin, nsMax, nsSum;
  int nssMin, nssMax, nssSum;
  int b;

  ncMax = Nmsg;
  ncMin = Nmsg;
  ncSum = Nmsg;
  comm_allreduce(&comm, gs_int, gs_max, &ncMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &ncMin, 1, &b);
  comm_allreduce(&comm, gs_int, gs_add, &ncSum, 1, &b);

  nsMax = Ncomm[0];
  nsMin = Ncomm[0];
  nsSum = Ncomm[0];
  for (i = 1; i < Nmsg; ++i) {
    nsMax = Ncomm[i] > Ncomm[i - 1] ? Ncomm[i] : Ncomm[i - 1];
    nsMin = Ncomm[i] < Ncomm[i - 1] ? Ncomm[i] : Ncomm[i - 1];
    nsSum += Ncomm[i];
  }
  comm_allreduce(&comm, gs_int, gs_max, &nsMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nsMin, 1, &b);

  nssMin = nsSum;
  nssMax = nsSum;
  nssSum = nsSum;
  comm_allreduce(&comm, gs_int, gs_max, &nssMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nssMin, 1, &b);
  comm_allreduce(&comm, gs_int, gs_add, &nssSum, 1, &b);

  if (Nmsg)
    nsSum = nsSum / Nmsg;
  else
    nsSum = 0;
  comm_allreduce(&comm, gs_int, gs_add, &nsSum, 1, &b);

  nelMax = nelt;
  nelMin = nelt;
  nelSum = nelt;
  comm_allreduce(&comm, gs_int, gs_max, &nelMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nelMin, 1, &b);
  comm_allreduce(&comm, gs_int, gs_add, &nelSum, 1, &b);

  free(Ncomm);
  comm_free(&comm);

  if (nc != NULL) {
    nc[0] = ncMin;
    nc[1] = ncMax;
    nc[2] = ncSum;
  }

  if (ns != NULL) {
    ns[0] = nsMin;
    ns[1] = nsMax;
    ns[2] = nsSum;
  }

  if (nss != NULL) {
    nss[0] = nssMin;
    nss[1] = nssMax;
    nss[2] = nssSum;
  }

  if (nel != NULL) {
    nel[0] = nelMin;
    nel[1] = nelMax;
    nel[2] = nelSum;
  }
}

void parrsb_print_part_stat(long long *vtx, unsigned nelt, unsigned nv,
                            MPI_Comm ce) {
  int id, np;
  MPI_Comm_rank(ce, &id);
  MPI_Comm_size(ce, &np);

  int nc[3], ns[3], nss[3], nel[3];
  parrsb_get_part_stat(&nc[0], &ns[0], &nss[0], &nel[0], vtx, nelt, nv, ce);

  if (id == 0) {
    printf("Min neighbors: %d | Max neighbors: %d | Avg neighbors: %lf\n",
           nc[0], nc[1], (double)nc[2] / np);
    printf("Min nvolume: %d | Max nvolume: %d | Avg nvolume: %lf\n", ns[0],
           ns[1], (double)ns[2] / np);
    printf("Min volume: %d | Max volume: %d | Avg volume: %lf\n", nss[0],
           nss[1], (double)nss[2] / np);
    printf("Min elements: %d | Max elements: %d\n", nel[0], nel[1]);
    fflush(stdout);
  }
}

struct parrsb_cmd_opts *parrsb_parse_cmd_opts(int argc, char *argv[]) {
  struct parrsb_cmd_opts *in = tcalloc(struct parrsb_cmd_opts, 1);

  in->mesh = NULL, in->tol = 2e-1;
  in->test = 0, in->dump = 1, in->verbose = 0, in->nactive = INT_MAX;
  in->ilu_type = 0, in->ilu_tol = 1e-1, in->ilu_pivot = 0;
  in->crs_type = 0, in->crs_tol = 1e-3;

  static struct option long_options[] = {
      {"mesh", required_argument, 0, 0},
      {"tol", optional_argument, 0, 1},
      {"test", optional_argument, 0, 2},
      {"dump", optional_argument, 0, 3},
      {"nactive", optional_argument, 0, 4},
      {"verbose", optional_argument, 0, 5},
      {"ilu_type", optional_argument, 0, 10},
      {"ilu_tol", optional_argument, 0, 11},
      {"ilu_pivot", optional_argument, 0, 12},
      {"crs_type", optional_argument, 0, 20},
      {"crs_tol", optional_argument, 0, 21},
      {0, 0, 0, 0}};

  size_t len;
  for (;;) {
    int c = getopt_long(argc, argv, "", long_options, NULL);
    if (c == -1)
      break;

    switch (c) {
    case 0:
      len = strnlen(optarg, PATH_MAX);
      in->mesh = tcalloc(char, len + 1);
      strncpy(in->mesh, optarg, len);
      break;
    case 1:
      in->tol = atof(optarg);
      break;
    case 2:
      in->test = atoi(optarg);
      break;
    case 3:
      in->dump = atoi(optarg);
      break;
    case 4:
      in->nactive = atoi(optarg);
      break;
    case 5:
      in->verbose = atoi(optarg);
      break;
    case 10:
      in->ilu_type = atoi(optarg);
      break;
    case 11:
      in->ilu_tol = atof(optarg);
      break;
    case 12:
      in->ilu_pivot = atoi(optarg);
      break;
    case 20:
      in->crs_type = atoi(optarg);
      break;
    case 21:
      in->crs_tol = atof(optarg);
      break;
    default:
      exit(1);
    }
  }

  return in;
}

void parrsb_cmd_opts_free(struct parrsb_cmd_opts *opts) {
  if (opts) {
    if (opts->mesh)
      free(opts->mesh);
    free(opts);
  }
}

void parrsb_check_error_(int err, char *file, int line, MPI_Comm comm) {
  int sum;
  MPI_Allreduce(&err, &sum, 1, MPI_INT, MPI_SUM, comm);

  if (sum != 0) {
    int id;
    MPI_Comm_rank(comm, &id);
    if (id == 0) {
      fprintf(stderr, "parrsb_check_error failure in %s:%d\n", file, line);
      fflush(stderr);
    }
    MPI_Finalize();
    exit(1);
  }
}

void parrsb_barrier(struct comm *c) {
#if defined(PARRSB_SYNC_BY_REDUCTION)
  sint dummy = c->id, wrk;
  comm_allreduce(c, gs_int, gs_max, &dummy, 1, &wrk);
#else
  comm_barrier(c);
#endif
}

#define WRITE_T(dest, val, T, nunits)                                          \
  do {                                                                         \
    memcpy(dest, (val), sizeof(T) * nunits);                                   \
    dest += sizeof(T) * nunits;                                                \
  } while (0)

int parrsb_vector_dump(const char *fname, scalar *y, struct rsb_element *elm,
                       uint nelt, unsigned nv, struct comm *c) {
  MPI_File file;
  sint err = MPI_File_open(c->c, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                           MPI_INFO_NULL, &file);
  slong wrk[2][1];
  comm_allreduce(c, gs_int, gs_add, &err, 1, wrk);

  uint rank = c->id;
  if (err) {
    if (rank == 0) {
      fprintf(stderr, "%s:%d Error opening file: %s\n", __FILE__, __LINE__,
              fname);
      fflush(stderr);
    }
    return err;
  }

  slong out[2][1], in = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0], nelgt = out[1][0];

  int ndim = (nv == 8) ? 3 : 2;
  uint write_size = ((ndim + 1) * sizeof(double) + sizeof(slong)) * nelt;
  if (rank == 0)
    write_size += sizeof(long) + sizeof(int); // for nelgt and ndim

  char *bfr, *bfr0;
  bfr = bfr0 = (char *)calloc(write_size, sizeof(char));
  if (rank == 0) {
    WRITE_T(bfr0, &nelgt, slong, 1);
    WRITE_T(bfr0, &ndim, int, 1);
  }

  uint i;
  for (i = 0; i < nelt; i++) {
    WRITE_T(bfr0, &elm[i].globalId, ulong, 1);
    WRITE_T(bfr0, &elm[i].coord[0], double, ndim);
    WRITE_T(bfr0, &y[i], double, 1);
  }

  MPI_Status st;
  err = MPI_File_write_ordered(file, bfr, write_size, MPI_BYTE, &st);
  comm_allreduce(c, gs_int, gs_add, &err, 1, wrk);
  if (err) {
    if (rank == 0) {
      fprintf(stderr, "%s:%d Error writing file: %s.\n", __FILE__, __LINE__,
              fname);
      fflush(stdout);
    }
    return err;
  }

  err = MPI_File_close(&file);
  comm_scan(out, c, gs_int, gs_add, &err, 1, wrk);
  if (err) {
    if (rank == 0) {
      fprintf(stderr, "%s:%d Error closing file: %s.\n", __FILE__, __LINE__,
              fname);
      fflush(stdout);
    }
    return err;
  }

  parrsb_barrier(c);
  free(bfr);
  return err;
}

#undef WRITE_T

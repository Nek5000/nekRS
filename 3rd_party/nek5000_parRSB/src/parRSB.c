#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <genmap-impl.h>
#include <parRSB.h>

parRSB_options parrsb_default_options = {0, -1, 0, 0, 0, 1, 1, 1};

void fparRSB_partMesh(int *part, int *seq, long long *vtx, double *coord,
                      int *nel, int *nv, int *options, int *comm, int *err) {
  *err = 1;
  comm_ext c = MPI_Comm_f2c(*comm);
  // TODO: Convert int options to parRSB_options instead of default options
  parRSB_options opt = parrsb_default_options;
  *err = parRSB_partMesh(part, seq, vtx, coord, *nel, *nv, &opt, c);
}

/*
 * part = [nel], out,
 * seq = [nel], out,
 * vtx = [nel x nv], in,
 * coord = [nel x nv x ndim], in,
 * nel = in,
 * nv = in,
 * options = in/out */
int parRSB_partMesh(int *part, int *seq, long long *vtx, double *coord, int nel,
                    int nv, parRSB_options *options, MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  int rank = c.id;
  int size = c.np;

  if (rank == 0)
    printf("running parRSB ...\n");
  fflush(stdout);

  comm_barrier(&c);
  double time0 = comm_time();

  struct crystal cr;
  crystal_init(&cr, &c);

  buffer bfr;
  buffer_init(&bfr, 1024);

  struct array eList;

  /* Load balance input data */
  genmap_load_balance(&eList, nel, nv, coord, vtx, &cr, &bfr);

  double time1 = comm_time();
  comm_barrier(&c);
  double time2 = comm_time();

  /* Run RSB now */
  comm_ext comm_rsb;
#ifdef MPI
  MPI_Comm_split(c.c, nel > 0, rank, &comm_rsb);
#endif

  // TODO: Move this into another file
  if (nel > 0) {
    metric_init();

    genmap_handle h;
    genmap_init(&h, comm_rsb, options);

    genmap_set_elements(h, &eList);
    genmap_comm_scan(h, genmap_global_comm(h));
    genmap_set_nvertices(h, nv);

    GenmapLong nelg = genmap_get_partition_nel(h);
    GenmapInt id = genmap_comm_rank(genmap_global_comm(h));
    GenmapInt size_ = genmap_comm_size(genmap_global_comm(h));

    if (size_ > nelg) {
      if (id == 0) {
        printf("Total number of elements is smaller than the "
               "number of processors.\n"
               "Run with smaller number of processors.\n");
      }
      // This is wrong
      return 1;
    }

    switch (options->global_partitioner) {
    case 0:
      genmap_rsb(h);
      break;
    case 1:
      genmap_rcb(h);
      break;
    case 2:
      genmap_rib(h);
      break;
    default:
      break;
    }

    genmap_finalize(h);

    if (options->print_timing_info > 0)
      metric_print(&c);
    metric_finalize();
  }

#ifdef MPI
  MPI_Comm_free(&comm_rsb);
#endif

  double time3 = comm_time();
  comm_barrier(&c);
  double time4 = comm_time();

  /* Restore original input */
  sarray_transfer(struct rsb_element, &eList, origin, 1, &cr);
  nel = eList.n;
  sarray_sort(struct rsb_element, eList.ptr, nel, globalId, 1, &bfr);

  struct rsb_element *e_ptr = eList.ptr;
  int e;
  for (e = 0; e < nel; e++)
    part[e] = e_ptr[e].origin;

  if (seq != NULL)
    for (e = 0; e < nel; e++)
      seq[e] = e_ptr[e].seq;

  double time5 = comm_time();
  comm_barrier(&c);
  double time = comm_time() - time0;

  /* Report time and finish */
  if (rank == 0)
    printf(" finished in %g s\n", time);

  if (options->print_timing_info > 0) {
    double min[3], max[3], sum[3], buf[3];
    min[0] = max[0] = sum[0] = time1 - time0;
    min[1] = max[1] = sum[1] = time3 - time2;
    min[2] = max[2] = sum[2] = time5 - time4;
    comm_allreduce(&c, gs_double, gs_min, min, 3, buf); // min
    comm_allreduce(&c, gs_double, gs_max, max, 3, buf); // max
    comm_allreduce(&c, gs_double, gs_add, sum, 3, buf); // sum
    if (rank == 0) {
      printf("LOADBALANCE : %g/%g/%g\n", min[0], max[0], sum[0] / c.np);
      printf("RSB         : %g/%g/%g\n", min[1], max[1], sum[1] / c.np);
      printf("RESTORE     : %g/%g/%g\n", min[2], max[2], sum[2] / c.np);
    }
    fflush(stdout);
  }

  array_free(&eList);
  buffer_free(&bfr);
  crystal_free(&cr);
  comm_free(&c);

  return 0;
}

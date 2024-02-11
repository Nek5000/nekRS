#include "crs_box_impl.hpp"

#define TIME_MAX 32
static int timer_on = 0;
static double timer_st;
static double time_box[TIME_MAX];
static unsigned num_solves = 0;

inline void timer_init() {
  num_solves = 0;
  for (unsigned i = 0; i < TIME_MAX; i++)
    time_box[i] = 0;

  timer_on = 1;
}

inline void timer_tic(const struct comm *c) {
  if (timer_on) {
    comm_barrier(c);
    timer_st = comm_time();
  }
}

inline void timer_toc(BOX_METRIC m) {
  if (timer_on && m != NONE) {
    double time = comm_time() - timer_st;
    time_box[m] += time;
  }
}

inline void timer_dump(struct comm *c, unsigned interval) {
  if (!timer_on || (num_solves != interval))
    return;

  struct btime_t {
    double tbox[TIME_MAX];
    uint p;
  };

  struct array btime;
  array_init(struct btime_t, &btime, 1);

  struct btime_t bt = {.p = 0};
  for (unsigned i = 0; i < TIME_MAX; i++)
    bt.tbox[i] = time_box[i];
  array_cat(struct btime_t, &btime, &bt, 1);

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct btime_t, &btime, p, 1, &cr);
  crystal_free(&cr);

  if (c->id == 0) {
    buffer bfr;
    buffer_init(&bfr, 1024);
    sarray_sort(struct btime_t, btime.ptr, btime.n, p, 0, &bfr);
    buffer_free(&bfr);

    FILE *fp = fopen("cumulative_time.txt", "w+");
    if (fp) {
      fprintf(fp, "copy_rhs,crs_dsavg1,asm1,mult_dsavg,mult_rhs_update,copy_to_"
                  "nek5000,map_vtx_to_box,asm2,map_box_to_vtx,copy_from_"
                  "nek5000,crs_dsavg2,copy_solution\n");
      fprintf(fp, "%d\n", btime.n);

      struct btime_t *ptr = (struct btime_t *)btime.ptr;
      for (uint i = 0; i < btime.n; i++) {
        fprintf(fp, "%d", ptr[i].p);
        for (uint j = 0; j < 12; j++)
          fprintf(fp, " %e", ptr[i].tbox[j] / interval);
        fprintf(fp, "\n");
      }
      fclose(fp);
    }
  }

  array_free(&btime);
}

inline void timer_print(struct comm *c, unsigned interval) {
  if (!timer_on)
    return;

  num_solves++;
  if (0 == num_solves % interval) {
    double max[TIME_MAX], wrk[2 * TIME_MAX];
    for (unsigned i = 0; i < TIME_MAX; i++)
      max[i] = time_box[i];
    comm_allreduce(c, gs_double, gs_max, max, TIME_MAX, wrk);
    if (c->id == 0) {
      printf("box copy_rhs          : %e\n", time_box[0] / num_solves);
      printf("box crs_dsavg1        : %e\n", time_box[1] / num_solves);
      printf("box asm1              : %e\n", time_box[2] / num_solves);
      printf("box crs_dsavg2        : %e\n", time_box[3] / num_solves);
      printf("box mult_rhs_update   : %e\n", time_box[4] / num_solves);
      printf("box copy_to_nek5000   : %e\n", time_box[5] / num_solves);
      printf("box map_vtx_to_box    : %e\n", time_box[6] / num_solves);
      printf("box asm2              : %e\n", time_box[7] / num_solves);
      printf("box map_box_to_vtx    : %e\n", time_box[8] / num_solves);
      printf("box copy_from_nek5000 : %e\n", time_box[9] / num_solves);
      printf("box crs_dsavg3        : %e\n", time_box[10] / num_solves);
      printf("box copy_solution     : %e\n", time_box[11] / num_solves);
      fflush(stdout);
    }
  }
}

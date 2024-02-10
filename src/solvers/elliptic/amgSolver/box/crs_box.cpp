#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// FIXME: Get rid of this once box solver is ported to C.
#include "crs_box_impl.hpp"
#include "nekInterfaceAdapter.hpp"

static void crs_box_dump(uint n, const ulong *id, uint nnz, const uint *Ai,
                         const uint *Aj, const double *A, uint null_space,
                         const struct comm *comm) {
  char file_name[BUFSIZ];
  snprintf(file_name, BUFSIZ, "crs_box_dump_%02d.txt", comm->id);

  FILE *fp = fopen(file_name, "w");
  if (!fp) {
    fprintf(stderr, "Failed to open file %s for writing.\n", file_name);
    fflush(stderr);
    MPI_Abort(comm->c, EXIT_FAILURE);
  }

  fprintf(fp, "%u %u %u\n", n, nnz, null_space);
  for (uint i = 0; i < n; i++)
    fprintf(fp, "%lu\n", id[i]);
  for (uint i = 0; i < nnz; i++)
    fprintf(fp, "%u %u %.*e\n", Ai[i], Aj[i], A[i]);

  fclose(fp);
}

static void dump_asm1(const char *name, const uint n, const ulong *const dofs,
                      const uint nnz, const uint *const Ai,
                      const uint *const Aj, const double *const Av,
                      const struct comm *const comm) {
  char file_name[BUFSIZ];
  snprintf(file_name, BUFSIZ, "%s_%02d.txt", name, comm->id);

  FILE *fp = fopen(file_name, "w");
  if (!fp) {
    fprintf(stderr, "Failed to open file %s for writing.\n", file_name);
    fflush(stderr);
    MPI_Abort(comm->c, EXIT_FAILURE);
  }

  fprintf(fp, "%u %u\n", n, nnz);
  for (uint i = 0; i < nnz; i++)
    fprintf(fp, "%u %u %.*e\n", dofs[Ai[i]], dofs[Aj[i]], Av[i]);

  fclose(fp);
}

static void dump_asm1_solution(const char *name, const uint n,
                               const double *const x, const uint ndim,
                               const double *const xyz,
                               const struct comm *const comm) {
  char file_name[BUFSIZ];
  snprintf(file_name, BUFSIZ, "%s_%02d.txt", name, comm->id);

  FILE *fp = fopen(file_name, "w");
  if (!fp) {
    fprintf(stderr, "Failed to open file %s for writing.\n", file_name);
    fflush(stderr);
    MPI_Abort(comm->c, EXIT_FAILURE);
  }

  const uint nv = (ndim == 3) ? 8 : 4;
  assert(n % nv == 0);

  for (uint e = 0; e < n / nv; e++) {
    const double *xyz_e = &xyz[e * nv * ndim];
    const double *x_e = &x[e * nv];
    fprintf(fp, "%d %e %e %e %e\n", e, x_e[0], xyz_e[0], xyz_e[1], xyz_e[2]);
    fprintf(fp, "%d %e %e %e %e\n", e, x_e[1], xyz_e[3], xyz_e[4], xyz_e[5]);
    fprintf(fp, "%d %e %e %e %e\n", e, x_e[2], xyz_e[9], xyz_e[10], xyz_e[11]);
    fprintf(fp, "%d %e %e %e %e\n", e, x_e[3], xyz_e[6], xyz_e[7], xyz_e[8]);

    if (nv == 4)
      continue;

    xyz_e = &xyz_e[3 * 4];
    x_e = &x_e[4];
    fprintf(fp, "%d %e %e %e %e\n", e, x_e[0], xyz_e[0], xyz_e[1], xyz_e[2]);
    fprintf(fp, "%d %e %e %e %e\n", e, x_e[1], xyz_e[3], xyz_e[4], xyz_e[5]);
    fprintf(fp, "%d %e %e %e %e\n", e, x_e[2], xyz_e[9], xyz_e[10], xyz_e[11]);
    fprintf(fp, "%d %e %e %e %e\n", e, x_e[3], xyz_e[6], xyz_e[7], xyz_e[8]);
  }

  fclose(fp);
}

void box_debug(const int verbose, const char *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  if (verbose > 0)
    vprintf(fmt, args);
  fflush(stdout);
  va_end(args);
}

static const sint *get_u2c(unsigned *cni, const unsigned n,
                           const ulong *const vtx, buffer *bfr) {
  struct vid_t {
    ulong id;
    uint idx;
    sint perm;
  };

  struct array vids;
  array_init(struct vid_t, &vids, n);

  struct vid_t vid;
  for (uint i = 0; i < n; i++) {
    vid.id = vtx[i], vid.idx = i;
    array_cat(struct vid_t, &vids, &vid, 1);
  }
  sarray_sort(struct vid_t, vids.ptr, vids.n, id, 1, bfr);

  struct vid_t *pv = (struct vid_t *)vids.ptr;
  ulong lid = 0;
  sint cn = 0;
  for (uint i = 0; i < vids.n; i++) {
    if (pv[i].id != lid)
      lid = pv[i].id, cn++;
    pv[i].perm = cn - 1;
  }
  *cni = cn;
  sarray_sort(struct vid_t, vids.ptr, vids.n, idx, 0, bfr);

  // Setup u2c -- user vector to compress vector mapping.
  pv = (struct vid_t *)vids.ptr;
  sint *const u2c = tcalloc(sint, n);
  for (uint i = 0; i < n; i++)
    u2c[i] = pv[i].perm;

  array_free(&vids);

  return u2c;
}

static void crs_box_setup_aux(struct box *box, uint ne, const long long *vtx,
                              const double *mask, const int *frontier, uint nw,
                              double tol, const struct comm *comm,
                              const MPI_Comm *inter_comm,
                              const double *const xyz) {
  const unsigned ncr = box->ncr, sn = box->sn, nnz = sn * ncr;
  const double *va = (const double *)nekData.schwz_amat;
  uint *ia = tcalloc(uint, nnz);
  uint *ja = tcalloc(uint, nnz);
  for (unsigned e = 0; e < ne; e++) {
    for (unsigned j = 0; j < ncr; j++) {
      for (unsigned i = 0; i < ncr; i++) {
        ia[e * ncr * ncr + j * ncr + i] = e * ncr + i;
        ja[e * ncr * ncr + j * ncr + i] = e * ncr + j;
      }
    }
  }

  ulong *tmp_vtx = tcalloc(ulong, box->sn);
  double *const tmp_mask = tcalloc(double, box->sn);
  double mask_min = DBL_MAX;
  for (unsigned i = 0; i < box->sn; i++) {
    tmp_vtx[i] = vtx[i];
    tmp_mask[i] = mask[i];
    if (frontier[i] == 1)
      tmp_mask[i] = 0;

    if (tmp_mask[i] < 0.1)
      tmp_vtx[i] = 0;
    if (tmp_mask[i] < mask_min)
      mask_min = tmp_mask[i];
  }
  free(tmp_mask);

  uint null_space = (mask_min < 1e-10) ? 0 : 1;
  assert(null_space == 0);

  box->cn = 0;
  box->u2c = NULL;
  box->ss = NULL;

  box->u2c = (int *)get_u2c(&box->cn, box->sn, tmp_vtx, &box->bfr);
  struct csr *A = csr_setup(nnz, ia, ja, va, box->u2c, tol, &box->bfr);
  asm1_setup(A, null_space, box);
  csr_free(A);

  free(ia), free(ja);

  // Setup the crs_dsavg which basically average the solution of original
  // parRSB domains.
  slong *gs_vtx = tcalloc(slong, box->sn);
  for (uint i = 0; i < box->un; i++)
    gs_vtx[i] = tmp_vtx[i];
  for (uint i = box->un; i < box->sn; i++)
    gs_vtx[i] = -tmp_vtx[i];
  box->gsh = gs_setup((const slong *)gs_vtx, box->sn, comm, 0, gs_auto, 0);
  free(gs_vtx), free(tmp_vtx);

  box->inv_mul = tcalloc(double, box->sn);
  for (uint i = 0; i < box->un; i++)
    box->inv_mul[i] = 1.0;
  gs(box->inv_mul, gs_double, gs_add, 0, box->gsh, &box->bfr);
  for (uint i = 0; i < box->sn; i++)
    box->inv_mul[i] = 1.0 / box->inv_mul[i];

// Allocate work arrays.
#define allocate_work_arrays(T)                                                \
  {                                                                            \
    box->sx = malloc(sizeof(T) * 2 * box->sn);                                 \
    box->srhs = (void *)((T *)box->sx + box->sn);                              \
  }
  BOX_DOMAIN_SWITCH(box->dom, allocate_work_arrays);
#undef allocate_work_arrays
}

struct box *crs_box_setup(uint n, const ulong *id, uint nnz, const uint *Ai,
                          const uint *Aj, const double *A, uint null_space,
                          const struct comm *comm, const MPI_Comm *inter_comm,
                          gs_dom dom) {
  struct box *box = tcalloc(struct box, 1);
  box->un = n;
  box->ncr = nnz / n;
  buffer_init(&box->bfr, 1024);

  const char *tmp = getenv("NEKRS_CRS_TIMER");
  if (tmp && atoi(tmp) > 0)
    timer_init();

  box->dom = dom;
  tmp = getenv("NEKRS_CRS_DOM");
  if (tmp && strncmp(tmp, "gs_double", 32) == 0)
    box->dom = gs_double;
  if (tmp && strncmp(tmp, "gs_float", 32) == 0)
    box->dom = gs_float;

  box->mult = 1;
  tmp = getenv("NEKRS_CRS_MULT");
  if (tmp)
    box->mult = atoi(tmp);

  box->algo = BOX_GPU_BLAS;
  tmp = getenv("NEKRS_CRS_ALGO");
  if (tmp)
    box->algo = atoi(tmp);

  // Copy the global communicator.
  comm_dup(&box->global, comm);

  // Copy the local communicator.
  MPI_Comm local;
  MPI_Comm_split(comm->c, comm->id, 1, &local);
  comm_init(&box->local, local);
  MPI_Comm_free(&local);

  // Fortran setup for ASM2. We should port this to C.
  nek::box_crs_setup();

  // Setup ASM1 solver on C side.
  box->sn = (*nekData.schwz_ne) * box->ncr;
  crs_box_setup_aux(
      box, (*nekData.schwz_ne), (const long long *)nekData.schwz_vtx,
      (const double *)nekData.schwz_mask, (const int *)nekData.schwz_frontier,
      (*nekData.schwz_nw), 1e-12, comm, inter_comm, nekData.schwz_xyz);

  // Print some info.
  if (box->global.id == 0) {
    printf("crs_box_setup: n = %u, nnz = %u, null_space = %u, dom = %s, "
           "mult = %u, algo = %u\n",
           box->un, nnz, null_space,
           (box->dom == gs_double) ? "gs_double" : "gs_float", box->mult,
           box->algo);
    fflush(stdout);
  }

  return box;
}

void crs_box_solve(void *x, struct box *box, const void *rhs) {
  struct comm *c = &box->global;

  // Copy RHS.
  timer_tic(c);
#define copy_rhs(T)                                                            \
  {                                                                            \
    const T *rhsi = (const T *)rhs;                                            \
    T *srhs = (T *)box->srhs;                                                  \
    for (uint i = 0; i < box->un; i++)                                         \
      srhs[i] = rhsi[i];                                                       \
  }
  BOX_DOMAIN_SWITCH(box->dom, copy_rhs);
#undef copy_rhs
  timer_toc(COPY_RHS);

  // crs_dsavg1.
  timer_tic(c);
  gs(box->srhs, box->dom, gs_add, 0, box->gsh, &box->bfr);
#define avg(T)                                                                 \
  {                                                                            \
    T *srhs = (T *)box->srhs;                                                  \
    for (uint i = 0; i < box->sn; i++)                                         \
      srhs[i] = box->inv_mul[i] * srhs[i];                                     \
  }
  BOX_DOMAIN_SWITCH(box->dom, avg);
#undef avg
  timer_toc(CRS_DSAVG1);

  // ASM1.
  timer_tic(c);
  switch (box->algo) {
  case BOX_GPU_BLAS:
    asm1_solve(box->sx, box, box->srhs);
    break;
  default:
    break;
  }
  timer_toc(ASM1);

  // crs_dsavg2.
  timer_tic(c);
  gs(box->sx, box->dom, gs_add, 0, box->gsh, &box->bfr);
#define avg(T)                                                                 \
  {                                                                            \
    T *sx = (T *)box->sx;                                                      \
    for (uint i = 0; i < box->sn; i++)                                         \
      sx[i] = box->inv_mul[i] * sx[i];                                         \
  }
  BOX_DOMAIN_SWITCH(box->dom, avg);
#undef avg
  timer_toc(CRS_DSAVG1);

  // mult_rhs_update.
  if (box->mult) {
    timer_tic(c);
// rhs = rhs - A*sx.
#define update_rhs(T)                                                          \
  {                                                                            \
    const double *A = (const double *)nekData.schwz_amat;                      \
    const T *sx = (T *)box->sx;                                                \
    T *srhs = (T *)box->srhs;                                                  \
    uint ncr = box->ncr, ue = box->un / ncr;                                   \
    for (uint e = 0; e < ue; e++) {                                            \
      for (uint c = 0; c < ncr; c++) {                                         \
        for (uint k = 0; k < ncr; k++) {                                       \
          srhs[k + ncr * e] -=                                                 \
              sx[c + ncr * e] * A[k + c * ncr + ncr * ncr * e];                \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  }
    BOX_DOMAIN_SWITCH(box->dom, update_rhs);
#undef update_rhs
    timer_toc(MULT_RHS_UPDATE);
  }

  // Copy to nek5000 to do the global solve.
  timer_tic(c);
#define copy_to_nek5000(T)                                                     \
  {                                                                            \
    const T *srhs = (T *)box->srhs;                                            \
    for (uint i = 0; i < box->un; i++)                                         \
      nekData.box_r[i] = srhs[i];                                              \
  }
  BOX_DOMAIN_SWITCH(box->dom, copy_to_nek5000);
#undef copy_to_nek5000
  timer_toc(COPY_TO_NEK5000);

  // Solve on nek5000.
  timer_tic(c);
  nek::box_map_vtx_to_box();
  timer_toc(MAP_VTX_TO_BOX);

  timer_tic(c);
  nek::box_crs_solve();
  timer_toc(ASM2);

  timer_tic(c);
  nek::box_map_box_to_vtx();
  timer_toc(MAP_BOX_TO_VTX);

  // Copy from nek5000.
  timer_tic(c);
#define copy_from_nek5000(T)                                                   \
  {                                                                            \
    T *sx = (T *)box->sx;                                                      \
    for (uint i = 0; i < box->un; i++)                                         \
      sx[i] += nekData.box_e[i];                                               \
  }
  BOX_DOMAIN_SWITCH(box->dom, copy_from_nek5000);
#undef copy_from_nek5000
  timer_toc(COPY_FROM_NEK5000);

  // crs_dsavg3.
  timer_tic(c);
  gs(box->sx, box->dom, gs_add, 0, box->gsh, &box->bfr);
#define avg(T)                                                                 \
  {                                                                            \
    T *sx = (T *)box->sx;                                                      \
    for (uint i = 0; i < box->un; i++)                                         \
      sx[i] = box->inv_mul[i] * sx[i];                                         \
  }
  BOX_DOMAIN_SWITCH(box->dom, avg);
#undef avg
  timer_toc(CRS_DSAVG1);

  // Copy solution.
  timer_tic(c);
#define copy_to_x(T)                                                           \
  {                                                                            \
    T *sx = (T *)box->sx, *xi = (T *)x;                                        \
    for (uint i = 0; i < box->un; i++)                                         \
      xi[i] = sx[i];                                                           \
  }
  BOX_DOMAIN_SWITCH(box->dom, copy_to_x);
#undef copy_to_x
  timer_toc(COPY_SOLUTION);

  timer_print(&box->global, 1000);
}

void crs_box_solve_no_gs(occa::memory &o_x, struct box *box,
                         occa::memory &o_rhs) {
  assert(box->algo == BOX_GPU_BLAS);
  assert(box->dom == gs_float);

  struct comm *c = &box->global;
  // ASM1.
  timer_tic(c);
  switch (box->algo) {
  case BOX_GPU_BLAS:
    asm1_solve((float *)box->sx, box, o_rhs);
    break;
  default:
    break;
  }
  timer_toc(ASM1);

  // crs_dsavg1.
  timer_tic(c);
  gs(box->sx, box->dom, gs_add, 0, box->gsh, &box->bfr);
#define avg(T)                                                                 \
  {                                                                            \
    T *sx = (T *)box->sx;                                                      \
    for (uint i = 0; i < box->sn; i++)                                         \
      sx[i] = box->inv_mul[i] * sx[i];                                         \
  }
  BOX_DOMAIN_SWITCH(box->dom, avg);
#undef avg
  timer_toc(CRS_DSAVG1);

  // copy the rhs
  o_rhs.copyTo(box->srhs, box->un * sizeof(float), 0);

  // mult_rhs_update: rhs = rhs - A*sx.
  if (box->mult) {
    timer_tic(c);
#define update_rhs(T)                                                          \
  {                                                                            \
    const double *A = (const double *)nekData.schwz_amat;                      \
    const T *sx = (T *)box->sx;                                                \
    T *srhs = (T *)box->srhs;                                                  \
    uint ncr = box->ncr, ue = box->un / ncr;                                   \
    for (uint e = 0; e < ue; e++) {                                            \
      for (uint c = 0; c < ncr; c++) {                                         \
        for (uint k = 0; k < ncr; k++) {                                       \
          srhs[k + ncr * e] -=                                                 \
              sx[c + ncr * e] * A[k + c * ncr + ncr * ncr * e];                \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  }
    BOX_DOMAIN_SWITCH(box->dom, update_rhs);
#undef update_rhs
    timer_toc(MULT_RHS_UPDATE);
  }

  // Copy to nek5000 to do the global solve.
  timer_tic(c);
#define copy_to_nek5000(T)                                                     \
  {                                                                            \
    const T *srhs = (T *)box->srhs;                                            \
    for (uint i = 0; i < box->un; i++)                                         \
      nekData.box_r[i] = srhs[i];                                              \
  }
  BOX_DOMAIN_SWITCH(box->dom, copy_to_nek5000);
#undef copy_to_nek5000
  timer_toc(COPY_TO_NEK5000);

  // Solve on nek5000.
  timer_tic(c);
  nek::box_map_vtx_to_box();
  timer_toc(MAP_VTX_TO_BOX);

  timer_tic(c);
  nek::box_crs_solve();
  timer_toc(ASM2);

  timer_tic(c);
  nek::box_map_box_to_vtx();
  timer_toc(MAP_BOX_TO_VTX);

  // Copy from nek5000.
  timer_tic(c);
#define copy_from_nek5000(T)                                                   \
  {                                                                            \
    T *sx = (T *)box->sx;                                                      \
    for (uint i = 0; i < box->un; i++)                                         \
      sx[i] += nekData.box_e[i];                                               \
  }
  BOX_DOMAIN_SWITCH(box->dom, copy_from_nek5000);
#undef copy_from_nek5000
  timer_toc(COPY_FROM_NEK5000);

  // crs_dsavg2.
  timer_tic(c);
  gs(box->sx, box->dom, gs_add, 0, box->gsh, &box->bfr);
#define avg(T)                                                                 \
  {                                                                            \
    T *sx = (T *)box->sx;                                                      \
    for (uint i = 0; i < box->un; i++)                                         \
      sx[i] = box->inv_mul[i] * sx[i];                                         \
  }
  BOX_DOMAIN_SWITCH(box->dom, avg);
#undef avg
  timer_toc(CRS_DSAVG1);

  // Copy solution.
  o_x.copyFrom(box->sx, box->un * sizeof(float), 0);

  // Print timing info.
  timer_print(&box->global, 1000);
}

void crs_box_solve2(occa::memory &o_x, struct box *box, occa::memory &o_rhs) {
  if (box->dom != gs_float) {
    fprintf(stderr, "Wrong domain !\n");
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  if (box->algo != BOX_GPU_BLAS) {
    fprintf(stderr, "Wrong solver !\n");
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  struct comm *c = &box->global;

  timer_tic(c);
  o_rhs.copyTo(box->srhs, box->un * sizeof(float), 0);
  timer_toc(COPY_RHS);

  // crs_dsavg1.
  timer_tic(c);
  gs(box->srhs, box->dom, gs_add, 0, box->gsh, &box->bfr);
#define avg(T)                                                                 \
  {                                                                            \
    T *srhs = (T *)box->srhs;                                                  \
    for (uint i = 0; i < box->sn; i++)                                         \
      srhs[i] = box->inv_mul[i] * srhs[i];                                     \
  }
  BOX_DOMAIN_SWITCH(box->dom, avg);
#undef avg
  timer_toc(CRS_DSAVG1);

  // ASM1.
  timer_tic(c);
  asm1_solve(box->sx, box, box->srhs);
  timer_toc(ASM1);

  // crs_dsavg2.
  timer_tic(c);
  gs(box->sx, box->dom, gs_add, 0, box->gsh, &box->bfr);
#define avg(T)                                                                 \
  {                                                                            \
    T *sx = (T *)box->sx;                                                      \
    for (uint i = 0; i < box->sn; i++)                                         \
      sx[i] = box->inv_mul[i] * sx[i];                                         \
  }
  BOX_DOMAIN_SWITCH(box->dom, avg);
#undef avg
  timer_toc(CRS_DSAVG1);

  // mult_rhs_update:  rhs = rhs - A*sx.
  timer_tic(c);
  if (box->mult) {
#define update_rhs(T)                                                          \
  {                                                                            \
    const double *A = (const double *)nekData.schwz_amat;                      \
    const T *sx = (T *)box->sx;                                                \
    T *srhs = (T *)box->srhs;                                                  \
    uint ncr = box->ncr, ue = box->un / ncr;                                   \
    for (uint e = 0; e < ue; e++) {                                            \
      for (uint c = 0; c < ncr; c++) {                                         \
        for (uint k = 0; k < ncr; k++) {                                       \
          srhs[k + ncr * e] -=                                                 \
              sx[c + ncr * e] * A[k + c * ncr + ncr * ncr * e];                \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  }
    BOX_DOMAIN_SWITCH(box->dom, update_rhs);
#undef update_rhs
  }
  timer_toc(MULT_RHS_UPDATE);

  // Copy to nek5000 to do the global solve.
  timer_tic(c);
#define copy_to_nek5000(T)                                                     \
  {                                                                            \
    const T *srhs = (T *)box->srhs;                                            \
    for (uint i = 0; i < box->un; i++)                                         \
      nekData.box_r[i] = srhs[i];                                              \
  }
  BOX_DOMAIN_SWITCH(box->dom, copy_to_nek5000);
#undef copy_to_nek5000
  timer_toc(COPY_TO_NEK5000);

  // Solve on nek5000.
  timer_tic(c);
  nek::box_map_vtx_to_box();
  timer_toc(MAP_VTX_TO_BOX);

  timer_tic(c);
  nek::box_crs_solve();
  timer_toc(ASM2);

  timer_tic(c);
  nek::box_map_box_to_vtx();
  timer_toc(MAP_BOX_TO_VTX);

  // Copy from nek5000.
  timer_tic(c);
#define copy_from_nek5000(T)                                                   \
  {                                                                            \
    T *sx = (T *)box->sx;                                                      \
    for (uint i = 0; i < box->un; i++)                                         \
      sx[i] += nekData.box_e[i];                                               \
  }
  BOX_DOMAIN_SWITCH(box->dom, copy_from_nek5000);
#undef copy_from_nek5000
  timer_toc(COPY_FROM_NEK5000);

  // crs_dsavg3.
  timer_tic(c);
  gs(box->sx, box->dom, gs_add, 0, box->gsh, &box->bfr);
#define avg(T)                                                                 \
  {                                                                            \
    T *sx = (T *)box->sx;                                                      \
    for (uint i = 0; i < box->un; i++)                                         \
      sx[i] = box->inv_mul[i] * sx[i];                                         \
  }
  BOX_DOMAIN_SWITCH(box->dom, avg);
#undef avg
  timer_toc(CRS_DSAVG1);

  // Copy solution.
  timer_tic(c);
  o_x.copyFrom(box->sx, box->un * sizeof(float), 0);
  timer_toc(COPY_SOLUTION);

  timer_print(&box->global, 1000);
}

void crs_box_free(struct box *box) {
  if (!box)
    return;

  switch (box->algo) {
  case BOX_GPU_BLAS:
    asm1_free(box);
    break;
  default:
    break;
  }

  gs_free(box->gsh);
  buffer_free(&box->bfr);
  comm_free(&box->local);
  comm_free(&box->global);
  free(box->u2c);
  free(box->inv_mul);
  free(box->sx);
}

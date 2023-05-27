#ifndef _PARRSB_MAT_H_
#define _PARRSB_MAT_H_

#include <gslib.h>

#ifdef scalar
#undef scalar
#endif
#define scalar double

#ifdef SCALAR_MAX
#undef SCALAR_MAX
#endif
#define SCALAR_MAX DBL_MAX

struct nbr {
  ulong r, c;
  uint proc;
};

struct mij {
  ulong r, c;
  uint idx, p;
  scalar v;
};

struct mat {
  ulong start;
  uint n, *Lp, *Li;
  scalar *L, *D;
};

struct par_mat {
  int type; // CSC or CSR
  uint cn, rn, *adj_off, *adj_idx, *diag_idx;
  scalar *adj_val, *diag_val;
  ulong *cols, *rows; // Unique global column and row ids
};

int IS_CSC(const struct par_mat *A);
int IS_CSR(const struct par_mat *A);
int IS_DIAG(const struct par_mat *A);

// Output array `arr` is an array of type `struct nbr`
void find_nbrs(struct array *arr, const ulong *eid, const slong *vtx,
               const uint nelt, const int nv, struct crystal *cr, buffer *buf);
// Output array `eij` is an array of type `struct mij`, input array `nbr` is
// an array of type `struct nbr`
int compress_nbrs(struct array *eij, struct array *nbr, buffer *bfr);

// Input array `entries` is of type `struct mij`
// TODO: Rename to mat_setup
int csr_setup(struct mat *mat, struct array *entries, int sep, buffer *buf);
int mat_print(struct mat *mat);
void mat_dump(const char *name, struct mat *A, struct crystal *cr, buffer *bfr);
int mat_free(struct mat *mat);

// Input array `entries` is of type `struct mij`
// type = 0 (CSC) or 1 (CSR)
int par_mat_setup(struct par_mat *M, struct array *mijs, const int type,
                  const int sd, buffer *bfr);
int par_csc_setup(struct par_mat *mat, struct array *entries, int sd,
                  buffer *buf);
int par_csr_setup(struct par_mat *mat, struct array *entries, int sd,
                  buffer *buf);
struct par_mat *par_csr_setup_ext(struct array *entries, int sd, buffer *bfr);
struct par_mat *par_csc_setup_ext(struct array *entries, int sd, buffer *bfr);
void par_csr_to_csc(struct par_mat *B, const struct par_mat *A, int diag,
                    struct crystal *cr, buffer *bfr);
void par_csc_to_csr(struct par_mat *B, const struct par_mat *A, int diag,
                    struct crystal *cr, buffer *bfr);
void par_mat_print(struct par_mat *A);
void par_mat_dump(const char *name, struct par_mat *A, struct crystal *cr,
                  buffer *bfr);
int par_mat_free(struct par_mat *A);

// Create a par_mat from connectivity
struct par_mat *par_csr_setup_con(const uint nelt, const ulong *eid,
                                  const slong *vtx, int nv, int sep,
                                  struct comm *c, struct crystal *cr,
                                  buffer *bfr);
// Mat vec routines
struct gs_data *setup_Q(const struct par_mat *M, const struct comm *c,
                        buffer *bfr);
void mat_vec_csr(scalar *y, const scalar *x, const struct par_mat *M,
                 struct gs_data *gsh, scalar *buf, buffer *bfr);

void par_arr_dump(const char *name, struct array *arr, struct crystal *cr,
                  buffer *bfr);
#endif // _MAT_H_

/*
 * Low-Order finite element preconditioner computed with HYPRE's AMG solver
*/

#include <math.h>
#include "gslib.h"
#include "_hypre_utilities.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_parcsr_ls.h"
#include "HYPRE.h"
#include "fem_amg_preco.hpp"

namespace{
typedef double (*Basis)(double*);
typedef void (*DBasis)(double**, double*);

void matrix_distribution();
void fem_assembly();
void quadrature_rule(double***, double**, int, int);
void mesh_connectivity(int***, int***, int, int);
void x_map(double**, double*, double**, int, Basis*);
void J_xr_map(double***, double*, double**, int, DBasis*);


double phi_3D_1(double*);
double phi_3D_2(double*);
double phi_3D_3(double*);
double phi_3D_4(double*);
void dphi_3D_1(double**, double*);
void dphi_3D_2(double**, double*);
void dphi_3D_3(double**, double*);
void dphi_3D_4(double**, double*);

double determinant(double**, int);
void inverse(double***, double**, int);
long long maximum(long long, long long);

int* mem_alloc_1D_int(int);
long long* mem_alloc_1D_long(int);
double* mem_alloc_1D_double(int);
int** mem_alloc_2D_int(int, int);
double** mem_alloc_2D_double(int, int);
void mem_free_1D_int(int**, int);
void mem_free_1D_long(long long**, int);
void mem_free_1D_double(double**, int);
void mem_free_2D_int(int***, int, int);
void mem_free_2D_double(double***, int, int);

static constexpr int n_dim = 3;
static int n_x, n_y, n_z, n_elem;
static int n_xyz, n_xyze;
static double *x_m, *y_m, *z_m;
static double *z_buffer;
static long long *glo_num;
static double *pmask;
static int num_loc_dofs;
static long long *dof_map;
static long long row_start;
static long long row_end;
static HYPRE_IJMatrix A_bc;
}
static struct comm comm;
struct gs_data {
  struct comm comm;
};
static struct gs_data *gsh;


/* Interface definition */
SEMFEMData* fem_amg_setup(const sint *n_x_, const sint *n_y_, const sint *n_z_,
                   const sint *n_elem_, double *x_m_, double *y_m_,
                   double *z_m_, double *pmask_, MPI_Comm mpiComm,
                   long long int *gatherGlobalNodes) {
  n_x = *n_x_;
  n_y = *n_y_;
  n_z = *n_z_;
  n_elem = *n_elem_;
  x_m = x_m_;
  y_m = y_m_;
  z_m = z_m_;
  pmask = pmask_;

  n_xyz = n_x * n_y * n_z;
  n_xyze = n_x * n_y * n_z * n_elem;
  z_buffer = (double *)calloc(n_xyze, sizeof(double));
  int NuniqueBases = n_xyze;

  {
    /* gslib stuff */
    comm_ext world;
    /*  MPI_Comm_dup(MPI_COMM_WORLD, (MPI_Comm*) &world); */
    world = (comm_ext)mpiComm; // MPI_COMM_WORLD;
    comm_init(&comm, world);
    gsh = gs_setup(gatherGlobalNodes, NuniqueBases, &comm, 0, gs_pairwise,
                   0); // gs_auto, gs_crystal_router, gs_pw
  }

  if (comm.id == 0)
    printf("fem_amg_setup ...\n");
  static_assert(sizeof(HYPRE_Int) == sizeof(long long),
                "HYPRE_Int must be long long!\n");

  double time0 = comm_time();
  matrix_distribution();
  fem_assembly();

  const long long numRows = row_end - row_start + 1;
  long long* ownedRows = (long long*) calloc(numRows, sizeof(long long));
  long long ctr = 0;
  for(long long row = row_start; row <= row_end; ++row)
    ownedRows[ctr++] = row;
  
  long long * ncols = (long long*) calloc(numRows, sizeof(long long));
  HYPRE_IJMatrixGetRowCounts(A_bc,
    numRows,
    ownedRows,
    ncols);

  long long nnz = 0;
  for(long long i = 0; i < numRows; ++i)
    nnz += ncols[i];
  
  // construct COO matrix from Hypre matrix
  long long *Ai = (long long*) calloc(nnz, sizeof(long long));
  long long *Aj = (long long*) calloc(nnz, sizeof(long long));
  double *Av = (double*) calloc(nnz, sizeof(double));
  HYPRE_IJMatrixGetValues(A_bc,
    -numRows,
    ncols,
    ownedRows,
    &Aj[0],
    &Av[0]);
  
  ctr = 0;
  for(long long i = 0; i < numRows; ++i){
    long long row = ownedRows[i];
    for(long long col = 0; col < ncols[i]; ++col)
      Ai[ctr++] = row;
  }

  HYPRE_IJMatrixDestroy(A_bc);

  SEMFEMData* data = new SEMFEMData();
  {
    data->Ai = Ai;
    data->Aj = Aj;
    data->Av = Av;
    data->nnz = nnz;
    data->rowStart = row_start;
    data->rowEnd = row_end;
    data->dofMap = dof_map;
  }

  double time1 = comm_time();
  if (comm.id == 0)
    printf("fem_amg_setup: done %fs\n", time1 - time0);
  fflush(stdout);

  return data;
}

namespace{

/* FEM Assembly definition */
void matrix_distribution() {
  /*
   * Ranks the global numbering array after removing the Dirichlet nodes
   * which is then used in the assembly of the matrices to map degrees of
   * freedom to rows of the matrix
   */

  int idx;
  buffer my_buffer;
  long long idx_start = n_xyze;
  long long scan_out[2], scan_buf[2];
  comm_scan(scan_out, &comm, gs_long_long, gs_add, &idx_start, 1, scan_buf);
  idx_start = scan_out[0];

  glo_num = mem_alloc_1D_long(n_xyze);

  for (idx = 0; idx < n_xyze; idx++) {
    if (pmask[idx] > 0.0)
      glo_num[idx] = idx_start + (long long)idx;
    else
      glo_num[idx] = -1;
  }

  gs(glo_num, gs_long_long, gs_min, 0, gsh, 0);

  /* Rank ids */
  long long maximum_value_local = 0;
  long long maximum_value = 0;

  for (idx = 0; idx < n_xyze; idx++) {
    maximum_value_local = (glo_num[idx] > maximum_value_local)
                              ? glo_num[idx]
                              : maximum_value_local;
  }

  comm_allreduce(&comm, gs_long_long, gs_max, &maximum_value_local, 1,
                 &maximum_value);
  const long long nstar = maximum_value / comm.np + 1;

  struct ranking_tuple {
    long long rank;
    unsigned int proc;
    unsigned int idx;
  };

  struct array ranking_transfer;
  array_init(ranking_tuple, &ranking_transfer, n_xyze);
  ranking_transfer.n = n_xyze;
  struct ranking_tuple *ranking_tuple_array =
      (struct ranking_tuple *)ranking_transfer.ptr;

  for (idx = 0; idx < ranking_transfer.n; idx++) {
    ranking_tuple_array[idx].rank = glo_num[idx];
    ranking_tuple_array[idx].proc = glo_num[idx] / nstar;
    ranking_tuple_array[idx].idx = idx;
  }

  struct crystal crystal_router_handle;
  crystal_init(&crystal_router_handle, &comm);
  sarray_transfer(ranking_tuple, &ranking_transfer, proc, 1,
                  &crystal_router_handle);
  ranking_tuple_array = (struct ranking_tuple *)ranking_transfer.ptr;

  buffer_init(&my_buffer, 1);
  sarray_sort(ranking_tuple, ranking_transfer.ptr, ranking_transfer.n, rank, 1,
              &my_buffer);

  long long current_rank = ranking_tuple_array[0].rank;
  long long current_count = 0;
  ranking_tuple_array[0].rank = current_count;

  for (idx = 1; idx < ranking_transfer.n; idx++) {

    if (ranking_tuple_array[idx].rank > current_rank) {
      current_count++;
      current_rank = ranking_tuple_array[idx].rank;
      ranking_tuple_array[idx].rank = current_count;
    } else if (ranking_tuple_array[idx].rank == current_rank) {
      ranking_tuple_array[idx].rank = current_count;
    } else {
      break;
    }
  }

  current_count += 1;

  long long rank_start;
  comm_scan(scan_out, &comm, gs_long_long, gs_add, &current_count, 1, scan_buf);
  rank_start = scan_out[0];

  for (idx = 0; idx < ranking_transfer.n; idx++) {
    ranking_tuple_array[idx].rank += rank_start;
  }

  sarray_transfer(ranking_tuple, &ranking_transfer, proc, 1,
                  &crystal_router_handle);
  ranking_tuple_array = (struct ranking_tuple *)ranking_transfer.ptr;

  buffer_init(&my_buffer, 1);
  sarray_sort(ranking_tuple, ranking_transfer.ptr, ranking_transfer.n, idx, 0,
              &my_buffer);

  for (idx = 0; idx < n_xyze; idx++) {
    glo_num[idx] = ranking_tuple_array[idx].rank;
  }

  array_free(&ranking_transfer);
  crystal_free(&crystal_router_handle);
}

void fem_assembly() {
  /*
   * Assembles the low-order FEM matrices from the spectral element mesh
   *
   * Returns A_fem and B_fem
   */

  /* Variables */
  int i, j, k, e, d, t, q;
  int idx;
  long long row;

  /*
   * Rank and prepare data to be mapped to rows of the matrix so it can
   * be solved with Hypre (Ranking done on the Fortran side since here it fails
   * for some reason I couldn't figure out)
   */
  long long *ranking = mem_alloc_1D_long(n_xyze);

  for (idx = 0; idx < n_xyze; idx++)
    ranking[idx] = glo_num[idx];

  row_start = 0;
  row_end = 0;

  for (idx = 0; idx < n_xyze; idx++)
    if (ranking[idx] >= 0)
      row_end = maximum(row_end, ranking[idx]);

  long long scan_out[2], scan_buf[2];
  comm_scan(scan_out, &comm, gs_long_long, gs_max, &row_end, 1, scan_buf);
  if (comm.id > 0)
    row_start = scan_out[0] + 1;

  num_loc_dofs = row_end - row_start + 1;

  dof_map = mem_alloc_1D_long(num_loc_dofs);

  for (idx = 0; idx < n_xyze; idx++) {
    if ((row_start <= ranking[idx]) && (ranking[idx] <= row_end)) {
      dof_map[ranking[idx] - row_start] = idx;
    }
  }

  /* Assemble FE matrices with boundary conditions applied */
  HYPRE_IJMatrixCreate(comm.c, row_start, row_end, row_start, row_end, &A_bc);
  HYPRE_IJMatrixSetObjectType(A_bc, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(A_bc);

  /* Set quadrature rule */
  int n_quad = 4;
  double **q_r;
  double *q_w;

  quadrature_rule(&q_r, &q_w, n_quad, n_dim);

  /* Mesh connectivity (Can be changed to fill-out or one-per-vertex) */
  int num_fem;
  int **v_coord;
  int **t_map;

  num_fem = 8;

  mesh_connectivity(&v_coord, &t_map, num_fem, n_dim);

  /* Finite element assembly */
  double **A_loc = mem_alloc_2D_double(n_dim + 1, n_dim + 1);
  double **J_xr = mem_alloc_2D_double(n_dim, n_dim);
  double **J_rx = mem_alloc_2D_double(n_dim, n_dim);
  double **x_t = mem_alloc_2D_double(n_dim, n_dim + 1);
  double *q_x = mem_alloc_1D_double(n_dim);
  double *dp = mem_alloc_1D_double(n_dim);

  Basis phi[n_dim + 1];
  DBasis dphi[n_dim + 1];

  phi[0] = phi_3D_1;
  phi[1] = phi_3D_2;
  phi[2] = phi_3D_3;
  phi[3] = phi_3D_4;
  dphi[0] = dphi_3D_1;
  dphi[1] = dphi_3D_2;
  dphi[2] = dphi_3D_3;
  dphi[3] = dphi_3D_4;

  int s_x, s_y, s_z;
  int E_x = n_x - 1;
  int E_y = n_y - 1;
  int E_z = n_z - 1;

  for (e = 0; e < n_elem; e++) {
    /* Cycle through collocated quads/hexes */
    for (s_z = 0; s_z < E_z; s_z++) {
      for (s_y = 0; s_y < E_y; s_y++) {
        for (s_x = 0; s_x < E_x; s_x++) {
          /* Get indices */
          int s[n_dim];

          s[0] = s_x;
          s[1] = s_y;
          s[2] = s_z;

          int idx[(int)(pow(2, n_dim))];

          for (i = 0; i < pow(2, n_dim); i++) {
            idx[i] = 0;

            for (d = 0; d < n_dim; d++) {
              idx[i] += (s[d] + v_coord[i][d]) * pow(n_x, d);
            }
          }

          /* Cycle through collocated triangles/tets */
          for (t = 0; t < num_fem; t++) {
            /* Get vertices */
            for (i = 0; i < n_dim + 1; i++) {
              if (n_dim == 2) {
                x_t[0][i] = x_m[idx[t_map[t][i]] + e * n_xyz];
                x_t[1][i] = y_m[idx[t_map[t][i]] + e * n_xyz];
              } else {
                x_t[0][i] = x_m[idx[t_map[t][i]] + e * n_xyz];
                x_t[1][i] = y_m[idx[t_map[t][i]] + e * n_xyz];
                x_t[2][i] = z_m[idx[t_map[t][i]] + e * n_xyz];
              }
            }

            /* Local FEM matrices */
            /* Reset local stiffness and mass matrices */
            for (i = 0; i < n_dim + 1; i++) {
              for (j = 0; j < n_dim + 1; j++) {
                A_loc[i][j] = 0.0;
              }
            }

            /* Build local stiffness matrices by applying quadrature rules */
            for (q = 0; q < n_quad; q++) {
              /* From r to x */
              x_map(&q_x, q_r[q], x_t, n_dim, phi);
              J_xr_map(&J_xr, q_r[q], x_t, n_dim, dphi);
              inverse(&J_rx, J_xr, n_dim);
              double det_J_xr = determinant(J_xr, n_dim);

              /* Integrand */
              for (i = 0; i < n_dim + 1; i++) {
                for (j = 0; j < n_dim + 1; j++) {
                  int alpha, beta;
                  double func = 0.0;

                  for (alpha = 0; alpha < n_dim; alpha++) {
                    double a = 0.0, b = 0.0;

                    for (beta = 0; beta < n_dim; beta++) {
                      dphi[i](&dp, q_r[q]);
                      a += dp[beta] * J_rx[beta][alpha];

                      dphi[j](&dp, q_r[q]);
                      b += dp[beta] * J_rx[beta][alpha];
                    }

                    func += a * b;
                  }

                  A_loc[i][j] += func * det_J_xr * q_w[q];
                }
              }
            }

            /* Add to global matrix */
            for (i = 0; i < n_dim + 1; i++) {
              for (j = 0; j < n_dim + 1; j++) {
                if ((pmask[idx[t_map[t][i]] + e * n_xyz] > 0.0) &&
                    (pmask[idx[t_map[t][j]] + e * n_xyz] > 0.0)) {
                  long long row = ranking[idx[t_map[t][i]] + e * n_xyz];
                  long long col = ranking[idx[t_map[t][j]] + e * n_xyz];

                  double A_val = A_loc[i][j];
                  long long ncols = 1;
                  int insert_error = 0;

                  if (fabs(A_val) > 1.0e-14) {
                    insert_error = HYPRE_IJMatrixAddToValues(
                        A_bc, 1, &ncols, &row, &col, &A_val);
                  }
                  if (insert_error != 0) {
                    if (comm.id == 0)
                      printf(
                          "There was an error with entry A(%lld, %lld) = %f\n",
                          row, col, A_val);
                    exit(EXIT_FAILURE);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  HYPRE_IJMatrixAssemble(A_bc);
  /* Free memory */
  mem_free_1D_long(&ranking, n_xyze);
  mem_free_1D_long(&glo_num, n_xyze);
  mem_free_2D_double(&q_r, n_quad, n_dim);
  mem_free_1D_double(&q_w, n_quad);
  mem_free_2D_int(&v_coord, pow(n_dim, 2), n_dim);
  mem_free_2D_int(&t_map, num_fem, n_dim + 1);
  mem_free_2D_double(&A_loc, n_dim + 1, n_dim + 1);
  mem_free_2D_double(&J_xr, n_dim, n_dim);
  mem_free_2D_double(&J_rx, n_dim, n_dim);
  mem_free_2D_double(&x_t, n_dim, n_dim + 1);
  mem_free_1D_double(&q_x, n_dim);
  mem_free_1D_double(&dp, n_dim);
}

void quadrature_rule(double ***q_r, double **q_w, int n_quad, int n_dim) {
  (*q_r) = mem_alloc_2D_double(n_quad, n_dim);
  (*q_w) = mem_alloc_1D_double(n_quad);

  if (n_quad == 4) {
    double a = (5.0 + 3.0 * sqrt(5.0)) / 20.0;
    double b = (5.0 - sqrt(5.0)) / 20.0;

    (*q_r)[0][0] = a;
    (*q_r)[0][1] = b;
    (*q_r)[0][2] = b;
    (*q_r)[1][0] = b;
    (*q_r)[1][1] = a;
    (*q_r)[1][2] = b;
    (*q_r)[2][0] = b;
    (*q_r)[2][1] = b;
    (*q_r)[2][2] = a;
    (*q_r)[3][0] = b;
    (*q_r)[3][1] = b;
    (*q_r)[3][2] = b;

    (*q_w)[0] = 1.0 / 24.0;
    (*q_w)[1] = 1.0 / 24.0;
    (*q_w)[2] = 1.0 / 24.0;
    (*q_w)[3] = 1.0 / 24.0;
  } else {
    printf("No quadrature rule for %d points available\n", n_quad);
    exit(EXIT_FAILURE);
  }
}

void mesh_connectivity(int ***v_coord, int ***t_map, int num_fem, int n_dim) {
  (*v_coord) = mem_alloc_2D_int(pow(n_dim, 2), n_dim);
  (*t_map) = mem_alloc_2D_int(num_fem, n_dim + 1);

  (*v_coord)[0][0] = 0;
  (*v_coord)[0][1] = 0;
  (*v_coord)[0][2] = 0;
  (*v_coord)[1][0] = 1;
  (*v_coord)[1][1] = 0;
  (*v_coord)[1][2] = 0;
  (*v_coord)[2][0] = 0;
  (*v_coord)[2][1] = 1;
  (*v_coord)[2][2] = 0;
  (*v_coord)[3][0] = 1;
  (*v_coord)[3][1] = 1;
  (*v_coord)[3][2] = 0;
  (*v_coord)[4][0] = 0;
  (*v_coord)[4][1] = 0;
  (*v_coord)[4][2] = 1;
  (*v_coord)[5][0] = 1;
  (*v_coord)[5][1] = 0;
  (*v_coord)[5][2] = 1;
  (*v_coord)[6][0] = 0;
  (*v_coord)[6][1] = 1;
  (*v_coord)[6][2] = 1;
  (*v_coord)[7][0] = 1;
  (*v_coord)[7][1] = 1;
  (*v_coord)[7][2] = 1;

  (*t_map)[0][0] = 0;
  (*t_map)[0][1] = 2;
  (*t_map)[0][2] = 1;
  (*t_map)[0][3] = 4;
  (*t_map)[1][0] = 1;
  (*t_map)[1][1] = 0;
  (*t_map)[1][2] = 3;
  (*t_map)[1][3] = 5;
  (*t_map)[2][0] = 2;
  (*t_map)[2][1] = 6;
  (*t_map)[2][2] = 3;
  (*t_map)[2][3] = 0;
  (*t_map)[3][0] = 3;
  (*t_map)[3][1] = 2;
  (*t_map)[3][2] = 7;
  (*t_map)[3][3] = 1;
  (*t_map)[4][0] = 4;
  (*t_map)[4][1] = 5;
  (*t_map)[4][2] = 6;
  (*t_map)[4][3] = 0;
  (*t_map)[5][0] = 5;
  (*t_map)[5][1] = 7;
  (*t_map)[5][2] = 4;
  (*t_map)[5][3] = 1;
  (*t_map)[6][0] = 6;
  (*t_map)[6][1] = 7;
  (*t_map)[6][2] = 2;
  (*t_map)[6][3] = 4;
  (*t_map)[7][0] = 7;
  (*t_map)[7][1] = 3;
  (*t_map)[7][2] = 6;
  (*t_map)[7][3] = 5;
}

void x_map(double **x, double *r, double **x_t, int n_dim, Basis *phi) {
  int i, d;

  for (d = 0; d < n_dim; d++) {
    (*x)[d] = 0.0;

    for (i = 0; i < n_dim + 1; i++) {
      (*x)[d] += x_t[d][i] * phi[i](r);
    }
  }
}

void J_xr_map(double ***J_xr, double *r, double **x_t, int n_dim,
              DBasis *dphi) {
  int i, j, k;
  double *deriv = mem_alloc_1D_double(n_dim);

  for (i = 0; i < n_dim; i++) {
    for (j = 0; j < n_dim; j++) {
      (*J_xr)[i][j] = 0.0;

      for (k = 0; k < n_dim + 1; k++) {
        dphi[k](&deriv, r);

        (*J_xr)[i][j] += x_t[i][k] * deriv[j];
      }
    }
  }

  mem_free_1D_double(&deriv, n_dim);
}

/* Basis functions and derivatives in 3D */
double phi_3D_1(double *r) { return r[0]; }
double phi_3D_2(double *r) { return r[1]; }
double phi_3D_3(double *r) { return r[2]; }
double phi_3D_4(double *r) { return 1.0 - r[0] - r[1] - r[2]; }
void dphi_3D_1(double **dp, double *r) {
  (*dp)[0] = 1.0;
  (*dp)[1] = 0.0;
  (*dp)[2] = 0.0;
}
void dphi_3D_2(double **dp, double *r) {
  (*dp)[0] = 0.0;
  (*dp)[1] = 1.0;
  (*dp)[2] = 0.0;
}
void dphi_3D_3(double **dp, double *r) {
  (*dp)[0] = 0.0;
  (*dp)[1] = 0.0;
  (*dp)[2] = 1.0;
}
void dphi_3D_4(double **dp, double *r) {
  (*dp)[0] = -1.0;
  (*dp)[1] = -1.0;
  (*dp)[2] = -1.0;
}

/* Math functions */
long long maximum(long long a, long long b) { return a > b ? a : b; }

double determinant(double **A, int n) {
  /*
   * Computes the determinant of a matrix
   */

  double d_1 = A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2]);
  double d_2 = A[0][1] * (A[1][0] * A[2][2] - A[2][0] * A[1][2]);
  double d_3 = A[0][2] * (A[1][0] * A[2][1] - A[2][0] * A[1][1]);

  return d_1 - d_2 + d_3;
}

void inverse(double ***inv_A, double **A, int n) {
  /*
   * Computes the inverse of a matrix
   */
  double det_A = determinant(A, n);

  (*inv_A)[0][0] = (1.0 / det_A) * (A[1][1] * A[2][2] - A[2][1] * A[1][2]);
  (*inv_A)[0][1] = (1.0 / det_A) * (A[0][2] * A[2][1] - A[2][2] * A[0][1]);
  (*inv_A)[0][2] = (1.0 / det_A) * (A[0][1] * A[1][2] - A[1][1] * A[0][2]);
  (*inv_A)[1][0] = (1.0 / det_A) * (A[1][2] * A[2][0] - A[2][2] * A[1][0]);
  (*inv_A)[1][1] = (1.0 / det_A) * (A[0][0] * A[2][2] - A[2][0] * A[0][2]);
  (*inv_A)[1][2] = (1.0 / det_A) * (A[0][2] * A[1][0] - A[1][2] * A[0][0]);
  (*inv_A)[2][0] = (1.0 / det_A) * (A[1][0] * A[2][1] - A[2][0] * A[1][1]);
  (*inv_A)[2][1] = (1.0 / det_A) * (A[0][1] * A[2][0] - A[2][1] * A[0][0]);
  (*inv_A)[2][2] = (1.0 / det_A) * (A[0][0] * A[1][1] - A[1][0] * A[0][1]);
}

/* Memory management */
int *mem_alloc_1D_int(int n) {
  int *array = (int *)malloc(n * sizeof(int));

  return array;
}

long long *mem_alloc_1D_long(int n) {
  long long *array = (long long *)malloc(n * sizeof(long long));

  return array;
}

double *mem_alloc_1D_double(int n) {
  double *array = (double *)malloc(n * sizeof(double));

  return array;
}

int **mem_alloc_2D_int(int n, int m) {
  int i;
  int **array = (int **)malloc(n * sizeof(int *));

  for (i = 0; i < n; i++)
    array[i] = (int *)malloc(m * sizeof(int));

  return array;
}

double **mem_alloc_2D_double(int n, int m) {
  int i;
  double **array = (double **)malloc(n * sizeof(double *));

  for (i = 0; i < n; i++)
    array[i] = (double *)malloc(m * sizeof(double));

  return array;
}

void mem_free_1D_int(int **array, int n) { free((*array)); }

void mem_free_1D_long(long long **array, int n) { free((*array)); }

void mem_free_1D_double(double **array, int n) { free((*array)); }

void mem_free_2D_int(int ***array, int n, int m) {
  int i;

  for (i = 0; i < n; i++)
    free((*array)[i]);

  free((*array));
}

void mem_free_2D_double(double ***array, int n, int m) {
  int i;

  for (i = 0; i < n; i++)
    free((*array)[i]);

  free((*array));
}
}
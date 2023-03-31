#include <math.h>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

#include "SEMFEMSolver.hpp"
#include "platform.hpp"
#include "gslib.h"

namespace {
void quadrature_rule(double q_r[4][3], double q_w[4])
{
  double a = (5.0 + 3.0 * sqrt(5.0)) / 20.0;
  double b = (5.0 - sqrt(5.0)) / 20.0;

  q_r[0][0] = a;
  q_r[0][1] = b;
  q_r[0][2] = b;
  q_r[1][0] = b;
  q_r[1][1] = a;
  q_r[1][2] = b;
  q_r[2][0] = b;
  q_r[2][1] = b;
  q_r[2][2] = a;
  q_r[3][0] = b;
  q_r[3][1] = b;
  q_r[3][2] = b;

  q_w[0] = 1.0 / 24.0;
  q_w[1] = 1.0 / 24.0;
  q_w[2] = 1.0 / 24.0;
  q_w[3] = 1.0 / 24.0;
}
long long bisection_search_index(const long long *sortedArr, long long value, long long start, long long end)
{
  long long fail = -1;
  long long L = start;
  long long R = end - 1;
  while (L <= R) {
    const long long m = (L + R) / 2;
    if (sortedArr[m] < value) {
      L = m + 1;
    }
    else if (sortedArr[m] > value) {
      R = m - 1;
    }
    else {
      return m;
    }
  }
  return fail;
}

long long linear_search_index(const long long *unsortedArr, long long value, long long start, long long end)
{
  long long fail = -1;
  for (long long idx = start; idx < end; ++idx) {
    if (unsortedArr[idx] == value) {
      return idx;
    }
  }
  return fail;
}

/* Basis functions and derivatives in 3D */
double phi_3D_1(double q_r[4][3], int q) { return q_r[q][0]; }
double phi_3D_2(double q_r[4][3], int q) { return q_r[q][1]; }
double phi_3D_3(double q_r[4][3], int q) { return q_r[q][2]; }
double phi_3D_4(double q_r[4][3], int q) { return 1.0 - q_r[q][0] - q_r[q][1] - q_r[q][2]; }
void dphi(double deriv[3], int q)
{
  if (q == 0) {
    deriv[0] = 1.0;
    deriv[1] = 0.0;
    deriv[2] = 0.0;
  }

  if (q == 1) {
    deriv[0] = 0.0;
    deriv[1] = 1.0;
    deriv[2] = 0.0;
  }

  if (q == 2) {
    deriv[0] = 0.0;
    deriv[1] = 0.0;
    deriv[2] = 1.0;
  }

  if (q == 3) {
    deriv[0] = -1.0;
    deriv[1] = -1.0;
    deriv[2] = -1.0;
  }
}

/* Math functions */
long long maximum(long long a, long long b) { return a > b ? a : b; }

double determinant(double A[3][3])
{
  /*
   * Computes the determinant of a matrix
   */

  double d_1 = A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2]);
  double d_2 = A[0][1] * (A[1][0] * A[2][2] - A[2][0] * A[1][2]);
  double d_3 = A[0][2] * (A[1][0] * A[2][1] - A[2][0] * A[1][1]);

  return d_1 - d_2 + d_3;
}

void inverse(double invA[3][3], double A[3][3])
{
  /*
   * Computes the inverse of a matrix
   */
  double inv_det_A = 1.0 / determinant(A);
  invA[0][0] = inv_det_A * (A[1][1] * A[2][2] - A[2][1] * A[1][2]);
  invA[0][1] = inv_det_A * (A[0][2] * A[2][1] - A[2][2] * A[0][1]);
  invA[0][2] = inv_det_A * (A[0][1] * A[1][2] - A[1][1] * A[0][2]);
  invA[1][0] = inv_det_A * (A[1][2] * A[2][0] - A[2][2] * A[1][0]);
  invA[1][1] = inv_det_A * (A[0][0] * A[2][2] - A[2][0] * A[0][2]);
  invA[1][2] = inv_det_A * (A[0][2] * A[1][0] - A[1][2] * A[0][0]);
  invA[2][0] = inv_det_A * (A[1][0] * A[2][1] - A[2][0] * A[1][1]);
  invA[2][1] = inv_det_A * (A[0][1] * A[2][0] - A[2][1] * A[0][0]);
  invA[2][2] = inv_det_A * (A[0][0] * A[1][1] - A[1][0] * A[0][1]);
}

void x_map(double x[3], double q_r[4][3], double x_t[3][4], int q)
{
  const int n_dim = 3;
  int i, d;

  for (d = 0; d < n_dim; d++) {
    x[d] = x_t[d][0] * phi_3D_1(q_r, q);
    x[d] += x_t[d][1] * phi_3D_2(q_r, q);
    x[d] += x_t[d][2] * phi_3D_3(q_r, q);
    x[d] += x_t[d][3] * phi_3D_4(q_r, q);
  }
}

void J_xr_map(double J_xr[3][3], double q_r[4][3], double x_t[3][4])
{
  const int n_dim = 3;
  int i, j, k;
  double deriv[3];

  for (i = 0; i < n_dim; i++) {
    for (j = 0; j < n_dim; j++) {
      J_xr[i][j] = 0.0;

      for (k = 0; k < n_dim + 1; k++) {
        dphi(deriv, k);

        J_xr[i][j] += x_t[i][k] * deriv[j];
      }
    }
  }
}

occa::memory scratchOrAllocateMemory(int nWords,
                                     int sizeT,
                                     void *src,
                                     size_t &bytesRemaining,
                                     size_t &byteOffset,
                                     size_t &bytesAllocated,
                                     bool &allocated);
static occa::kernel computeStiffnessMatrixKernel;
static occa::memory o_stiffness;
static occa::memory o_x;
static occa::memory o_y;
static occa::memory o_z;
static bool constructOnHost = false;

void load();

void construct_coo_graph();
void fem_assembly_device(hypreWrapper::IJ_t &hypreIJ);
void fem_assembly_host(hypreWrapper::IJ_t &hypreIJ);

void matrix_distribution();
void fem_assembly(hypreWrapper::IJ_t &hypreIJ);
void mesh_connectivity(int[8][3], int[8][4]);
long long maximum(long long, long long);

static constexpr int n_dim = 3;
static int n_x, n_y, n_z, n_elem;
static int n_xyz, n_xyze;
static long long *glo_num;
static double *pmask;
static double lambda0;
static int num_loc_dofs;
static long long *dof_map;
static long long row_start;
static long long row_end;

struct COOGraph {
  int nrows;
  long long nnz;
  long long *rows;
  long long *rowOffsets;
  int *ncols;
  long long *cols;
  float *vals;
};

static COOGraph coo_graph;

} // namespace

static struct comm comm;
struct gs_data {
  struct comm comm;
};
static struct gs_data *gsh;

/* Interface definition */
SEMFEMSolver_t::matrix_t *SEMFEMSolver_t::build(const int N_,
                                const int n_elem_,
                                occa::memory _o_x,
                                occa::memory _o_y,
                                occa::memory _o_z,
                                double *pmask_,
                                double lambda0_,
                                hypreWrapper::IJ_t &hypreIJ,
                                MPI_Comm mpiComm,
                                long long int *gatherGlobalNodes)
{
  n_x = N_;
  n_y = N_;
  n_z = N_;
  o_x = _o_x;
  o_y = _o_y;
  o_z = _o_z;
  n_elem = n_elem_;
  pmask = pmask_;
  lambda0 = lambda0_;

  n_xyz = n_x * n_y * n_z;
  n_xyze = n_x * n_y * n_z * n_elem;
  int NuniqueBases = n_xyze;

  {
    comm_ext world;
    world = (comm_ext)mpiComm;
    comm_init(&comm, world);
    gsh = gs_setup(gatherGlobalNodes,
                   NuniqueBases,
                   &comm,
                   0,
                   gs_pairwise,
                   /* mode */ 0);
  }

  constructOnHost = !platform->device.deviceAtomic;

  if (!constructOnHost)
    load();

  matrix_distribution();

  fem_assembly(hypreIJ);

  auto matrix = new matrix_t();

  {

    const int numRows = row_end - row_start + 1;

    {
      long long numRowsGlobal64;
      long long numRows64 = numRows;
      comm_allreduce(&comm, gs_long_long, gs_add, &numRows64, 1, &numRowsGlobal64);
      nrsCheck(numRowsGlobal64 > std::numeric_limits<int>::max(), comm.c, EXIT_FAILURE,
               "Number of global rows requires BigInt support!", "");
    }

    hypreWrapper::BigInt *ownedRows = (hypreWrapper::BigInt *)calloc(numRows, sizeof(hypreWrapper::BigInt));
    int ctr = 0;
    for (long long row = row_start; row <= row_end; ++row)
      ownedRows[ctr++] = row;

    hypreWrapper::Int *ncols = (hypreWrapper::Int *)calloc(numRows, sizeof(hypreWrapper::Int));
    hypreIJ.MatrixGetRowCounts(numRows, ownedRows, ncols);

    int nnz = 0;
    for (int i = 0; i < numRows; ++i)
      nnz += ncols[i];

    // construct COO matrix from Hypre matrix
    hypreWrapper::BigInt *hAj = (hypreWrapper::BigInt *)calloc(nnz, sizeof(hypreWrapper::BigInt));
    hypreWrapper::Real *hAv = (hypreWrapper::Real *)calloc(nnz, sizeof(hypreWrapper::Real));
    hypreIJ.MatrixGetValues(-numRows, ncols, ownedRows, &hAj[0], &hAv[0]);

    long long *Ai = (long long *)calloc(nnz, sizeof(long long));
    long long *Aj = (long long *)calloc(nnz, sizeof(long long));
    double *Av = (double *)calloc(nnz, sizeof(double));
    for (int n = 0; n < nnz; ++n) {
      Aj[n] = hAj[n];
      Av[n] = hAv[n];
    }
    ctr = 0;
    for (int i = 0; i < numRows; ++i) {
      long long row = ownedRows[i];
      for (int col = 0; col < ncols[i]; ++col)
        Ai[ctr++] = row;
    }

    free(hAj);
    free(hAv);
    free(ownedRows);
    free(ncols);

    double dropTol = 0.0;
    platform->options.getArgs("AMG DROP TOLERANCE", dropTol);

    int nnzTol = 0;
    for (int n = 0; n < nnz; ++n) {
      if (std::abs(Av[n]) > dropTol) {
        nnzTol++;
      }
    }

    long long *AiTol = (long long *)calloc(nnzTol, sizeof(long long));
    long long *AjTol = (long long *)calloc(nnzTol, sizeof(long long));
    double *AvTol = (double *)calloc(nnzTol, sizeof(double));

    int idx = 0;
    for (int n = 0; n < nnz; ++n) {
      if (std::abs(Av[n]) > dropTol) {
        AiTol[idx] = Ai[n];
        AjTol[idx] = Aj[n];
        AvTol[idx] = Av[n];
        idx++;
      }
    }

    free(Ai);
    free(Aj);
    free(Av);

    matrix->Ai = AiTol;
    matrix->Aj = AjTol;
    matrix->Av = AvTol;
    matrix->nnz = nnzTol;
    matrix->rowStart = row_start;
    matrix->rowEnd = row_end;
    matrix->dofMap = dof_map;
  }

  return matrix;
}

namespace {

/* FEM Assembly definition */
void matrix_distribution()
{
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

  glo_num = (long long *)malloc(n_xyze * sizeof(long long));

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
    maximum_value_local = (glo_num[idx] > maximum_value_local) ? glo_num[idx] : maximum_value_local;
  }

  comm_allreduce(&comm, gs_long_long, gs_max, &maximum_value_local, 1, &maximum_value);
  const long long nstar = maximum_value / comm.np + 1;

  struct ranking_tuple {
    long long rank;
    unsigned int proc;
    unsigned int idx;
  };

  struct array ranking_transfer;
  array_init(ranking_tuple, &ranking_transfer, n_xyze);
  ranking_transfer.n = n_xyze;
  struct ranking_tuple *ranking_tuple_array = (struct ranking_tuple *)ranking_transfer.ptr;

  for (idx = 0; idx < ranking_transfer.n; idx++) {
    ranking_tuple_array[idx].rank = glo_num[idx];
    ranking_tuple_array[idx].proc = glo_num[idx] / nstar;
    ranking_tuple_array[idx].idx = idx;
  }

  struct crystal crystal_router_handle;
  crystal_init(&crystal_router_handle, &comm);
  sarray_transfer(ranking_tuple, &ranking_transfer, proc, 1, &crystal_router_handle);
  ranking_tuple_array = (struct ranking_tuple *)ranking_transfer.ptr;

  buffer_init(&my_buffer, 1);
  sarray_sort(ranking_tuple, ranking_transfer.ptr, ranking_transfer.n, rank, 1, &my_buffer);

  long long current_rank = ranking_tuple_array[0].rank;
  long long current_count = 0;
  ranking_tuple_array[0].rank = current_count;

  for (idx = 1; idx < ranking_transfer.n; idx++) {

    if (ranking_tuple_array[idx].rank > current_rank) {
      current_count++;
      current_rank = ranking_tuple_array[idx].rank;
      ranking_tuple_array[idx].rank = current_count;
    }
    else if (ranking_tuple_array[idx].rank == current_rank) {
      ranking_tuple_array[idx].rank = current_count;
    }
    else {
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

  sarray_transfer(ranking_tuple, &ranking_transfer, proc, 1, &crystal_router_handle);
  ranking_tuple_array = (struct ranking_tuple *)ranking_transfer.ptr;

  buffer_init(&my_buffer, 1);
  sarray_sort(ranking_tuple, ranking_transfer.ptr, ranking_transfer.n, idx, 0, &my_buffer);

  for (idx = 0; idx < n_xyze; idx++) {
    glo_num[idx] = ranking_tuple_array[idx].rank;
  }

  array_free(&ranking_transfer);
  crystal_free(&crystal_router_handle);
}
void construct_coo_graph()
{

  /* Mesh connectivity (Can be changed to fill-out or one-per-vertex) */
  constexpr int num_fem = 8;
  int v_coord[8][3];
  int t_map[8][4];

  mesh_connectivity(v_coord, t_map);

  int E_x = n_x - 1;
  int E_y = n_y - 1;
  int E_z = n_z - 1;

  std::unordered_map<long long, std::unordered_set<long long>> graph;
  std::unordered_map<int, std::unordered_set<int>> rowIdxToColIdxMap;
  const int nvert = 8;
  for (int s_z = 0; s_z < E_z; s_z++) {
    for (int s_y = 0; s_y < E_y; s_y++) {
      for (int s_x = 0; s_x < E_x; s_x++) {
        /* Get indices */
        int s[n_dim];

        s[0] = s_x;
        s[1] = s_y;
        s[2] = s_z;

        int idx[nvert];

        for (int i = 0; i < nvert; i++) {
          idx[i] = 0;

          idx[i] += (s[0] + v_coord[i][0]) * 1;
          idx[i] += (s[1] + v_coord[i][1]) * n_x;
          idx[i] += (s[2] + v_coord[i][2]) * n_x * n_x;
        }
        for (int t = 0; t < num_fem; t++) {
          for (int i = 0; i < n_dim + 1; i++) {
            for (int j = 0; j < n_dim + 1; j++) {
              rowIdxToColIdxMap[idx[t_map[t][i]]].insert(idx[t_map[t][j]]);
            }
          }
        }
      }
    }
  }
  for (int e = 0; e < n_elem; ++e) {
    for (auto &&idx_row_lookup_and_col_lookups : rowIdxToColIdxMap) {
      auto row_lookup = idx_row_lookup_and_col_lookups.first;
      for (auto &&col_lookup : idx_row_lookup_and_col_lookups.second) {
        if (pmask[e * n_xyz + row_lookup] > 0.0 && pmask[e * n_xyz + col_lookup] > 0.0) {
          long long row = glo_num[e * n_xyz + row_lookup];
          long long col = glo_num[e * n_xyz + col_lookup];
          graph[row].emplace(col);
        }
      }
    }
  }
  const int nrows = graph.size();
  long long *rows = (long long *)malloc(nrows * sizeof(long long));
  long long *rowOffsets = (long long *)malloc((nrows + 1) * sizeof(long long));
  int *ncols = (int *)malloc(nrows * sizeof(int));
  long long nnz = 0;
  long long ctr = 0;

  for (auto &&row_and_colset : graph) {
    rows[ctr++] = row_and_colset.first;
    nnz += row_and_colset.second.size();
  }
  long long *cols = (long long *)malloc(nnz * sizeof(long long));
  float *vals = (float *)calloc(nnz, sizeof(float));
  std::sort(rows, rows + nrows);
  long long entryCtr = 0;
  rowOffsets[0] = 0;
  for (auto localrow = 0; localrow < nrows; ++localrow) {
    const long long row = rows[localrow];
    const auto &colset = graph[row];
    const auto size = colset.size();
    ncols[localrow] = size;
    rowOffsets[localrow + 1] = rowOffsets[localrow] + size;
    for (auto &&col : colset) {
      cols[entryCtr++] = col;
    }
  }

  coo_graph.nrows = nrows;
  coo_graph.nnz = nnz;
  coo_graph.rows = rows;
  coo_graph.rowOffsets = rowOffsets;
  coo_graph.ncols = ncols;
  coo_graph.vals = vals;
  coo_graph.cols = cols;
}

void fem_assembly_host(hypreWrapper::IJ_t &hypreIJ)
{

  const int nrows = coo_graph.nrows;
  long long *rows = coo_graph.rows;
  long long *rowOffsets = coo_graph.rowOffsets;
  int *ncols = coo_graph.ncols;
  long long nnz = coo_graph.nnz;
  long long *cols = coo_graph.cols;
  float *vals = coo_graph.vals;

  double q_r[4][3];
  double q_w[4];
  int v_coord[8][3];
  int t_map[8][4];
  quadrature_rule(q_r, q_w);
  mesh_connectivity(v_coord, t_map);

  double *x = (double *)malloc(n_xyze * sizeof(double));
  double *y = (double *)malloc(n_xyze * sizeof(double));
  double *z = (double *)malloc(n_xyze * sizeof(double));
  o_x.copyTo(x, n_xyze * sizeof(double));
  o_y.copyTo(y, n_xyze * sizeof(double));
  o_z.copyTo(z, n_xyze * sizeof(double));

  const int n_quad = 4;
  const int num_fem = 8;

  for (int e = 0; e < n_elem; e++) {

    /* Cycle through collocated quads/hexes */
    for (int s_z = 0; s_z < n_x - 1; s_z++) {
      for (int s_y = 0; s_y < n_x - 1; s_y++) {
        for (int s_x = 0; s_x < n_x - 1; s_x++) {
          double A_loc[4][4];
          double J_xr[3][3];
          double J_rx[3][3];
          double x_t[3][4];
          double q_x[3];
          /* Get indices */
          int s[n_dim];

          s[0] = s_x;
          s[1] = s_y;
          s[2] = s_z;

          int idx[8];

          for (int i = 0; i < 8; i++) {
            idx[i] = 0;
            idx[i] += (s[0] + v_coord[i][0]) * 1;
            idx[i] += (s[1] + v_coord[i][1]) * n_x;
            idx[i] += (s[2] + v_coord[i][2]) * n_x * n_x;
          }

          /* Cycle through collocated triangles/tets */
          for (int t = 0; t < num_fem; t++) {
            /* Get vertices */
            for (int i = 0; i < n_dim + 1; i++) {
              x_t[0][i] = x[idx[t_map[t][i]] + e * n_x * n_x * n_x];
              x_t[1][i] = y[idx[t_map[t][i]] + e * n_x * n_x * n_x];
              x_t[2][i] = z[idx[t_map[t][i]] + e * n_x * n_x * n_x];
            }

            /* Local FEM matrices */
            /* Reset local stiffness and mass matrices */
            for (int i = 0; i < n_dim + 1; i++) {
              for (int j = 0; j < n_dim + 1; j++) {
                A_loc[i][j] = 0.0;
              }
            }

            /* Build local stiffness matrices by applying quadrature rules */
            J_xr_map(J_xr, q_r, x_t);
            inverse(J_rx, J_xr);
            const double det_J_xr = determinant(J_xr);
            for (int q = 0; q < n_quad; q++) {
              /* From r to x */
              x_map(q_x, q_r, x_t, q);

              /* Integrand */
              for (int i = 0; i < n_dim + 1; i++) {
                double deriv_i[3];
                dphi(deriv_i, i);
                for (int j = 0; j < n_dim + 1; j++) {
                  double deriv_j[3];
                  dphi(deriv_j, j);
                  int alpha, beta;
                  double func = 0.0;

                  for (alpha = 0; alpha < n_dim; alpha++) {
                    double a = 0.0, b = 0.0;

                    for (beta = 0; beta < n_dim; beta++) {
                      a += deriv_i[beta] * J_rx[beta][alpha];

                      b += deriv_j[beta] * J_rx[beta][alpha];
                    }

                    func += a * b;
                  }

                  A_loc[i][j] += func * det_J_xr * q_w[q];
                }
              }
            }
            for (int i = 0; i < n_dim + 1; i++) {
              for (int j = 0; j < n_dim + 1; j++) {
                if ((pmask[idx[t_map[t][i]] + e * n_x * n_x * n_x] > 0.0) &&
                    (pmask[idx[t_map[t][j]] + e * n_x * n_x * n_x] > 0.0)) {
                  long long row = glo_num[idx[t_map[t][i]] + e * n_x * n_x * n_x];
                  long long col = glo_num[idx[t_map[t][j]] + e * n_x * n_x * n_x];
                  long long local_row_id = bisection_search_index(rows, row, 0, nrows);
                  long long start = rowOffsets[local_row_id];
                  long long end = rowOffsets[local_row_id + 1];

                  long long id = linear_search_index(cols, col, start, end);
                  vals[id] += lambda0 * A_loc[i][j];
                }
              }
            }
          }
        }
      }
    }
  }

  int err = hypreIJ.MatrixAddToValues(nrows, ncols, rows, cols, vals);
  nrsCheck(err != 0, comm.c, EXIT_FAILURE, 
           "hypreWrapper::IJMatrixAddToValues failed!\n", "");

  free(rows);
  free(rowOffsets);
  free(ncols);
  free(cols);
  free(vals);
  free(x);
  free(y);
  free(z);
}
void fem_assembly_device(hypreWrapper::IJ_t &hypreIJ)
{

  const int nrows = coo_graph.nrows;
  long long *rows = coo_graph.rows;
  long long *rowOffsets = coo_graph.rowOffsets;
  int *ncols = coo_graph.ncols;
  long long nnz = coo_graph.nnz;
  long long *cols = coo_graph.cols;
  float *vals = coo_graph.vals;

  struct AllocationTracker {
    bool o_maskAlloc;
    bool o_glo_numAlloc;
    bool o_rowOffsetsAlloc;
    bool o_rowsAlloc;
    bool o_colsAlloc;
    bool o_valsAlloc;
  };
  AllocationTracker allocations;
  size_t bytesRemaining = platform->o_mempool.bytesAllocated;
  size_t byteOffset = 0;
  size_t bytesAllocated = 0;
  occa::memory o_mask = scratchOrAllocateMemory(n_xyze,
                                                sizeof(double),
                                                pmask,
                                                bytesRemaining,
                                                byteOffset,
                                                bytesAllocated,
                                                allocations.o_maskAlloc);
  occa::memory o_glo_num = scratchOrAllocateMemory(n_xyze,
                                                   sizeof(long long),
                                                   glo_num,
                                                   bytesRemaining,
                                                   byteOffset,
                                                   bytesAllocated,
                                                   allocations.o_glo_numAlloc);
  occa::memory o_rows = scratchOrAllocateMemory(nrows,
                                                sizeof(long long),
                                                rows,
                                                bytesRemaining,
                                                byteOffset,
                                                bytesAllocated,
                                                allocations.o_rowsAlloc);
  occa::memory o_rowOffsets = scratchOrAllocateMemory(nrows + 1,
                                                      sizeof(long long),
                                                      rowOffsets,
                                                      bytesRemaining,
                                                      byteOffset,
                                                      bytesAllocated,
                                                      allocations.o_rowOffsetsAlloc);
  occa::memory o_cols = scratchOrAllocateMemory(nnz,
                                                sizeof(long long),
                                                cols,
                                                bytesRemaining,
                                                byteOffset,
                                                bytesAllocated,
                                                allocations.o_colsAlloc);
  occa::memory o_vals = scratchOrAllocateMemory(nnz,
                                                sizeof(float),
                                                vals,
                                                bytesRemaining,
                                                byteOffset,
                                                bytesAllocated,
                                                allocations.o_valsAlloc);

  computeStiffnessMatrixKernel(n_elem,
                               (int)nrows,
                               o_x,
                               o_y,
                               o_z,
                               o_mask,
                               o_glo_num,
                               o_rows,
                               o_rowOffsets,
                               o_cols,
                               o_vals);
  o_vals.copyTo(vals, nnz * sizeof(float));

  if (allocations.o_maskAlloc)
    o_mask.free();
  if (allocations.o_glo_numAlloc)
    o_glo_num.free();
  if (allocations.o_rowOffsetsAlloc)
    o_rowOffsets.free();
  if (allocations.o_rowsAlloc)
    o_rows.free();
  if (allocations.o_colsAlloc)
    o_cols.free();
  if (allocations.o_valsAlloc)
    o_vals.free();

  int err = hypreIJ.MatrixAddToValues(nrows, ncols, rows, cols, vals);
  nrsCheck(err != 0, comm.c, EXIT_FAILURE,
           "hypreWrapper::IJMatrixAddToValues failed!\n", "");

  free(rows);
  free(rowOffsets);
  free(ncols);
  free(cols);
  free(vals);
}

void fem_assembly(hypreWrapper::IJ_t &hypreIJ)
{
  /*
   * Assembles the low-order FEM matrices from the spectral element mesh
   *
   * Returns A_fem and B_fem
   */

  /* Variables */
  double tStart = MPI_Wtime();
  if (comm.id == 0)
    printf("building matrix ... ");
  int i, j, k, e, d, t, q;
  int idx;
  long long row;

  row_start = 0;
  row_end = 0;

  for (idx = 0; idx < n_xyze; idx++)
    if (glo_num[idx] >= 0)
      row_end = maximum(row_end, glo_num[idx]);

  long long scan_out[2], scan_buf[2];
  comm_scan(scan_out, &comm, gs_long_long, gs_max, &row_end, 1, scan_buf);
  if (comm.id > 0)
    row_start = scan_out[0] + 1;

  num_loc_dofs = row_end - row_start + 1;

  dof_map = (long long *)malloc(num_loc_dofs * sizeof(long long));

  for (idx = 0; idx < n_xyze; idx++) {
    if ((row_start <= glo_num[idx]) && (glo_num[idx] <= row_end)) {
      dof_map[glo_num[idx] - row_start] = idx;
    }
  }

  /* Assemble FE matrices with boundary conditions applied */
  hypreIJ.MatrixCreate(comm.c, row_start, row_end, row_start, row_end);
  hypreIJ.MatrixSetObjectType();
  hypreIJ.MatrixInitialize();

  construct_coo_graph();

  if (constructOnHost) {
    fem_assembly_host(hypreIJ);
  }
  else {
    fem_assembly_device(hypreIJ);
  }

  {
    hypreIJ.MatrixAssemble();
  }

  free(glo_num);

  MPI_Barrier(comm.c);
  if (comm.id == 0)
    printf("done (%gs)\n", MPI_Wtime() - tStart);
}

void load() { computeStiffnessMatrixKernel = platform->kernels.get("computeStiffnessMatrix"); }
void mesh_connectivity(int v_coord[8][3], int t_map[8][4])
{

  (v_coord)[0][0] = 0;
  (v_coord)[0][1] = 0;
  (v_coord)[0][2] = 0;
  (v_coord)[1][0] = 1;
  (v_coord)[1][1] = 0;
  (v_coord)[1][2] = 0;
  (v_coord)[2][0] = 0;
  (v_coord)[2][1] = 1;
  (v_coord)[2][2] = 0;
  (v_coord)[3][0] = 1;
  (v_coord)[3][1] = 1;
  (v_coord)[3][2] = 0;
  (v_coord)[4][0] = 0;
  (v_coord)[4][1] = 0;
  (v_coord)[4][2] = 1;
  (v_coord)[5][0] = 1;
  (v_coord)[5][1] = 0;
  (v_coord)[5][2] = 1;
  (v_coord)[6][0] = 0;
  (v_coord)[6][1] = 1;
  (v_coord)[6][2] = 1;
  (v_coord)[7][0] = 1;
  (v_coord)[7][1] = 1;
  (v_coord)[7][2] = 1;

  (t_map)[0][0] = 0;
  (t_map)[0][1] = 2;
  (t_map)[0][2] = 1;
  (t_map)[0][3] = 4;
  (t_map)[1][0] = 1;
  (t_map)[1][1] = 0;
  (t_map)[1][2] = 3;
  (t_map)[1][3] = 5;
  (t_map)[2][0] = 2;
  (t_map)[2][1] = 6;
  (t_map)[2][2] = 3;
  (t_map)[2][3] = 0;
  (t_map)[3][0] = 3;
  (t_map)[3][1] = 2;
  (t_map)[3][2] = 7;
  (t_map)[3][3] = 1;
  (t_map)[4][0] = 4;
  (t_map)[4][1] = 5;
  (t_map)[4][2] = 6;
  (t_map)[4][3] = 0;
  (t_map)[5][0] = 5;
  (t_map)[5][1] = 7;
  (t_map)[5][2] = 4;
  (t_map)[5][3] = 1;
  (t_map)[6][0] = 6;
  (t_map)[6][1] = 7;
  (t_map)[6][2] = 2;
  (t_map)[6][3] = 4;
  (t_map)[7][0] = 7;
  (t_map)[7][1] = 3;
  (t_map)[7][2] = 6;
  (t_map)[7][3] = 5;
}

occa::memory scratchOrAllocateMemory(int nWords,
                                     int sizeT,
                                     void *src,
                                     size_t &bytesRemaining,
                                     size_t &byteOffset,
                                     size_t &bytesAllocated,
                                     bool &allocated)
{
  occa::memory o_mem;
  if (nWords * sizeT < bytesRemaining) {
    o_mem = platform->o_mempool.o_ptr.slice(byteOffset);
    o_mem.copyFrom(src, nWords * sizeT);
    bytesRemaining -= nWords * sizeT;
    byteOffset += nWords * sizeT;
    allocated = false;
  }
  else {
    o_mem = platform->device.malloc(nWords * sizeT, src);
    allocated = true;
    bytesAllocated += nWords * sizeT;
  }
  return o_mem;
}

} // namespace

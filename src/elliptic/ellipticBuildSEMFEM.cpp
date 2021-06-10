/*
 * Low-Order finite element preconditioner computed with HYPRE's AMG solver
*/

#include <platform.hpp>
#include <math.h>
#include <limits>

#include "_hypre_utilities.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_parcsr_ls.h"
#include "HYPRE.h"

#include "gslib.h"
#include "ellipticBuildSEMFEM.hpp"
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

namespace{

occa::memory scratchOrAllocateMemory(int nWords, int sizeT, void* src, long long& bytesRemaining, long long& byteOffset, long long& bytesAllocated, bool& allocated);
static occa::kernel computeStiffnessMatrixKernel;
static occa::memory o_stiffness;
static occa::memory o_x;
static occa::memory o_y;
static occa::memory o_z;

void build_kernel();

void fem_assembly_device();

void matrix_distribution();
void fem_assembly();
void mesh_connectivity(int[8][3], int[8][4]);
long long maximum(long long, long long);

static constexpr int n_dim = 3;
static int n_x, n_y, n_z, n_elem;
static int n_xyz, n_xyze;
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
SEMFEMData* ellipticBuildSEMFEM(const int N_, const int n_elem_,
                          occa::memory _o_x, occa::memory _o_y, occa::memory _o_z,
                          double *pmask_, MPI_Comm mpiComm,
                          long long int *gatherGlobalNodes
                          )
{
  n_x = N_;
  n_y = N_;
  n_z = N_;
  o_x = _o_x;
  o_y = _o_y;
  o_z = _o_z;
  n_elem = n_elem_;
  pmask = pmask_;

  n_xyz = n_x * n_y * n_z;
  n_xyze = n_x * n_y * n_z * n_elem;
  int NuniqueBases = n_xyze;

  {
    comm_ext world;
    world = (comm_ext)mpiComm; // MPI_COMM_WORLD;
    comm_init(&comm, world);
    gsh = gs_setup(gatherGlobalNodes, NuniqueBases, &comm, 0, gs_pairwise,
                   /* mode */ 0);
  }

  build_kernel();

  if(platform->options.compareArgs("BUILD ONLY", "TRUE")) return NULL;

  matrix_distribution();

  fem_assembly();

  SEMFEMData* data;

  {

    const int numRows = row_end - row_start + 1;

    {
      long long numRowsGlobal64;
      long long numRows64 = numRows;
      comm_allreduce(&comm, gs_long_long, gs_add, &numRows64, 1, &numRowsGlobal64);
      if(numRowsGlobal64 > std::numeric_limits<int>::max()) { 
        if(comm.id == 0) printf("Number of global rows requires BigInt support!");
        MPI_Abort(comm.c, EXIT_FAILURE);  
      }
    }

    HYPRE_BigInt *ownedRows = (HYPRE_BigInt*) calloc(numRows, sizeof(HYPRE_BigInt));
    int ctr = 0;
    for(long long row = row_start; row <= row_end; ++row)
      ownedRows[ctr++] = row;
  
    HYPRE_Int *ncols = (HYPRE_Int*) calloc(numRows, sizeof(HYPRE_Int));
    HYPRE_IJMatrixGetRowCounts(A_bc,
      numRows,
      ownedRows,
      ncols);

    int nnz = 0;
    for(int i = 0; i < numRows; ++i)
      nnz += ncols[i];
  
    // construct COO matrix from Hypre matrix
    HYPRE_BigInt *hAj = (HYPRE_BigInt*) calloc(nnz, sizeof(HYPRE_BigInt));
    HYPRE_Real   *hAv = (HYPRE_Real*) calloc(nnz, sizeof(HYPRE_Real));
    HYPRE_IJMatrixGetValues(A_bc,
      -numRows,
      ncols,
      ownedRows,
      &hAj[0],
      &hAv[0]);

    long long *Ai = (long long*) calloc(nnz, sizeof(long long));
    long long *Aj = (long long*) calloc(nnz, sizeof(long long));
    double    *Av = (double*) calloc(nnz, sizeof(double));
    for(int n = 0; n < nnz; ++n) {
       Aj[n] = hAj[n];
       Av[n] = hAv[n];
    } 
    ctr = 0;
    for(int i = 0; i < numRows; ++i){
      long long row = ownedRows[i];
      for(int col = 0; col < ncols[i]; ++col)
        Ai[ctr++] = row;
    }

    free(hAj);
    free(hAv);
    free(ownedRows);
    free(ncols);
    HYPRE_IJMatrixDestroy(A_bc);

    data = (SEMFEMData*) malloc(sizeof(SEMFEMData));
    data->Ai = Ai;
    data->Aj = Aj;
    data->Av = Av;
    data->nnz = nnz;
    data->rowStart = row_start;
    data->rowEnd = row_end;
    data->dofMap = dof_map;

  }

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

  glo_num = (long long*) malloc(n_xyze * sizeof(long long));

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
void fem_assembly_device() {

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

          for (int d = 0; d < n_dim; d++) {
            idx[i] += (s[d] + v_coord[i][d]) * pow(n_x, d);
          }
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
  for(int e = 0; e < n_elem; ++e){
    for(auto&& idx_row_lookup_and_col_lookups : rowIdxToColIdxMap){
      auto row_lookup = idx_row_lookup_and_col_lookups.first;
      for(auto && col_lookup : idx_row_lookup_and_col_lookups.second)
      {
        if(pmask[e * n_xyz + row_lookup] > 0.0 && pmask[e * n_xyz + col_lookup] > 0.0)
        {
          long long row = glo_num[e*n_xyz + row_lookup];
          long long col = glo_num[e*n_xyz + col_lookup];
          graph[row].emplace(col);
        }
      }
    }
  }
  const long long nrows = graph.size();
  long long * rows = (long long*) malloc(nrows * sizeof(long long));
  long long * rowOffsets = (long long*) malloc((nrows+1) * sizeof(long long));
  long long * ncols = (long long*) malloc(nrows * sizeof(long long));
  long long nnz = 0;
  long long ctr = 0;

  for(auto&& row_and_colset : graph){
    rows[ctr++] = row_and_colset.first;
    nnz += row_and_colset.second.size();
  }
  long long * cols = (long long*) malloc(nnz * sizeof(long long));
  double* vals = (double*) calloc(nnz,sizeof(double));
  std::sort(rows, rows + nrows);
  long long entryCtr = 0;
  rowOffsets[0] = 0;
  for(long long localrow = 0; localrow < nrows; ++localrow){
    const long long row = rows[localrow];
    const auto& colset = graph[row];
    const int size = colset.size();
    ncols[localrow] = size;
    rowOffsets[localrow+1] = rowOffsets[localrow] + size;
    for(auto&& col : colset){
      cols[entryCtr++] = col;
    }
  }

  struct AllocationTracker{
    bool o_maskAlloc;
    bool o_glo_numAlloc;
    bool o_rowOffsetsAlloc;
    bool o_rowsAlloc;
    bool o_colsAlloc;
    bool o_valsAlloc;
  };
  AllocationTracker allocations;
  long long bytesRemaining = platform->o_mempool.bytesAllocated;
  long long byteOffset = 0;
  long long bytesAllocated = 0;
  occa::memory o_mask = scratchOrAllocateMemory(
    n_xyze,
    sizeof(double),
    pmask,
    bytesRemaining,
    byteOffset,
    bytesAllocated,
    allocations.o_maskAlloc
  );
  occa::memory o_glo_num = scratchOrAllocateMemory(
    n_xyze,
    sizeof(long long),
    glo_num,
    bytesRemaining,
    byteOffset,
    bytesAllocated,
    allocations.o_glo_numAlloc
  );
  occa::memory o_rows = scratchOrAllocateMemory(
    nrows,
    sizeof(long long),
    rows,
    bytesRemaining,
    byteOffset,
    bytesAllocated,
    allocations.o_rowsAlloc
  );
  occa::memory o_rowOffsets = scratchOrAllocateMemory(
    nrows+1,
    sizeof(long long),
    rowOffsets,
    bytesRemaining,
    byteOffset,
    bytesAllocated,
    allocations.o_rowOffsetsAlloc
  );
  occa::memory o_cols = scratchOrAllocateMemory(
    nnz,
    sizeof(long long),
    cols,
    bytesRemaining,
    byteOffset,
    bytesAllocated,
    allocations.o_colsAlloc
  );
  occa::memory o_vals = scratchOrAllocateMemory(
    nnz,
    sizeof(double),
    vals,
    bytesRemaining,
    byteOffset,
    bytesAllocated,
    allocations.o_valsAlloc
  );

  MPI_Allreduce(MPI_IN_PLACE, &bytesAllocated, 1, MPI_LONG_LONG, MPI_SUM, platform->comm.mpiComm);
  double bytesTotal = (double) bytesAllocated / 1e9;
  double bytesPerProc = bytesTotal / platform->comm.mpiCommSize;

  computeStiffnessMatrixKernel(
    n_elem,
    (int)nrows,
    o_x,
    o_y,
    o_z,
    o_mask,
    o_glo_num,
    o_rows,
    o_rowOffsets,
    o_cols,
    o_vals
  );
  o_vals.copyTo(vals, nnz * sizeof(double));

  if(allocations.o_maskAlloc) o_mask.free();
  if(allocations.o_glo_numAlloc) o_glo_num.free();
  if(allocations.o_rowOffsetsAlloc) o_rowOffsets.free();
  if(allocations.o_rowsAlloc) o_rows.free();
  if(allocations.o_colsAlloc) o_cols.free();
  if(allocations.o_valsAlloc) o_vals.free();

  int err = HYPRE_IJMatrixAddToValues(A_bc, nrows, ncols, rows, cols, vals);
  if (err != 0) {
    if (comm.id == 0)
      printf("err!\n");
    exit(EXIT_FAILURE);
  }

  free(rows);
  free(rowOffsets);
  free(ncols);
  free(cols);
  free(vals);
}

void fem_assembly() {
  /*
   * Assembles the low-order FEM matrices from the spectral element mesh
   *
   * Returns A_fem and B_fem
   */

  /* Variables */
  double tStart = MPI_Wtime();
  if(comm.id == 0) printf("building matrix ... ");
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

  dof_map = (long long *) malloc(num_loc_dofs * sizeof(long long));

  for (idx = 0; idx < n_xyze; idx++) {
    if ((row_start <= glo_num[idx]) && (glo_num[idx] <= row_end)) {
      dof_map[glo_num[idx] - row_start] = idx;
    }
  }

  /* Assemble FE matrices with boundary conditions applied */
  HYPRE_IJMatrixCreate(comm.c, row_start, row_end, row_start, row_end, &A_bc);
  HYPRE_IJMatrixSetObjectType(A_bc, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(A_bc);

  {
    fem_assembly_device();
  }

  {
    HYPRE_IJMatrixAssemble(A_bc);
  }

  free(glo_num);

  MPI_Barrier(comm.c);
  if(comm.id == 0) printf("done (%gs)\n", MPI_Wtime() - tStart);
}

void build_kernel(){
  std::string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  std::string oklpath = install_dir + "/okl/";
  occa::properties stiffnessKernelInfo = platform->kernelInfo;
  std::string filename = oklpath + "elliptic/ellipticSEMFEMStiffness.okl";
  stiffnessKernelInfo["defines/" "p_Nq"] = n_x;
  stiffnessKernelInfo["defines/" "p_Np"] = n_x * n_x * n_x;
  stiffnessKernelInfo["defines/" "p_rows_sorted"] = 1;
  stiffnessKernelInfo["defines/" "p_cols_sorted"] = 0;

  computeStiffnessMatrixKernel = platform->device.buildKernel(
    filename,
    "computeStiffnessMatrix",
    stiffnessKernelInfo
  );
}
void mesh_connectivity(int v_coord[8][3], int t_map[8][4]) {

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

occa::memory scratchOrAllocateMemory(int nWords, int sizeT, void* src, long long& bytesRemaining, long long& byteOffset, long long& bytesAllocated, bool& allocated)
{
  occa::memory o_mem;
  if(nWords * sizeT < bytesRemaining){
    o_mem = platform->o_mempool.o_ptr.slice(byteOffset);
    o_mem.copyFrom(src, nWords * sizeT);
    bytesRemaining -= nWords * sizeT;
    byteOffset += nWords * sizeT;
    allocated = false;
  } else {
    o_mem = platform->device.malloc(nWords * sizeT, src);
    allocated = true;
    bytesAllocated += nWords * sizeT;
  }
  return o_mem;
}
long long maximum(long long a, long long b) { return a > b ? a : b; }

}

#include <stdio.h>
#include <string.h>
#include <mpi.h>

#include "AMGX.hpp"

#ifdef ENABLE_AMGX

static int AMGXinit = 0;

void printHandler(const char *msg, int length)
{
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if(myid == 0) printf("%s", msg);
}

AMGX_t::AMGX_t(const int nLocalRows_, const int nnz,
               const long long *rows, const long long *cols, const double *values, /* COO */ 
               const int nullspace, const MPI_Comm comm_, int deviceID,
               int useFP32, int MPI_DIRECT, const char* cfgFile)
{
  MPI_Comm_dup(comm_,&comm);
  nLocalRows = nLocalRows_;
 
  int myid, commSize;
  MPI_Comm_rank(comm, &myid);
  MPI_Comm_size(comm, &commSize);

  const int useDevice = (deviceID < 0) ? 0 : 1; 

  AMGX_Mode mode;
  if(useDevice && !useFP32)
    mode = AMGX_mode_dDDI;
  else if (useDevice && useFP32)
    mode = AMGX_mode_dFFI;
  else if (!useDevice && !useFP32)
    mode = AMGX_mode_hDDI;
  else if (useDevice && useFP32)
    mode = AMGX_mode_hFFI;

  if(!AMGXinit) {
    AMGX_SAFE_CALL(AMGX_initialize());
    AMGX_SAFE_CALL(AMGX_initialize_plugins());
    AMGX_SAFE_CALL(AMGX_register_print_callback(&printHandler));
  }

  if (cfgFile) { 
    AMGX_SAFE_CALL(AMGX_config_create_from_file(&cfg, cfgFile));
  } else {
    char solverSettings[] = 
      "\"solver\": {\
            \"scope\":\"main\",\
            \"solver\":\"AMG\",\
            \"algorithm\":\"CLASSICAL\",\
            \"strength_threshold\":0.25,\
            \"max_row_sum\":0.9,\
            \"interpolator\":\"D2\",\
            \"interp_max_elements\":4,\
            \"max_levels\":20,\
            \"print_config\":0,\
            \"print_grid_stats\":1,\
            \"max_iters\":1,\
            \"cycle\":\"V\",\
            \"presweeps\":1,\
            \"postsweeps\":1,\
            \"coarsest_sweeps\":3,\
            \"use_sum_stopping_criteria\":1,\
            \"coarse_solver\": \"NOSOLVER\"\
      }";

    char cfgStr[4096] = "";
    strcat(cfgStr, "{ \"config_version\": 2,");
    strcat(cfgStr, (MPI_DIRECT) ? "\"communicator\":\"MPI_DIRECT\"," : "\"communicator\":\"MPI\",");
    strcat(cfgStr, solverSettings);
    strcat(cfgStr, "}");
    //printf("cfgStr: %s\n", cfgStr); fflush(stdout);
    AMGX_SAFE_CALL(AMGX_config_create(&cfg, cfgStr));
  }

  AMGX_resources_create(&rsrc, cfg, &comm, 1, &deviceID);

  AMGX_solver_create(&solver, rsrc, mode, cfg);

  AMGX_vector_create(&AmgXP, rsrc, mode);
  AMGX_vector_create(&AmgXRHS, rsrc, mode);
  AMGX_matrix_create(&AmgXA, rsrc, mode);

  AMGX_distribution_handle dist;
  AMGX_distribution_create(&dist, cfg);

  long long *partitionOffsets = (long long*) calloc(commSize + 1, sizeof(long long));
  long long n64 = nLocalRows;
  MPI_Allgather(&n64, 1, MPI_LONG_LONG, &partitionOffsets[1], 1, MPI_LONG_LONG, comm);
  for (int r=2; r<commSize+1; r++) partitionOffsets[r] += partitionOffsets[r-1]; 
  const long long nGlobalRows = partitionOffsets[commSize]; 
  AMGX_distribution_set_partition_data(dist, AMGX_DIST_PARTITION_OFFSETS, partitionOffsets);

  int *csrRows = (int*) calloc(nLocalRows + 1, sizeof(int));
  for (int i = 0; i < nnz; i++) {
    int powPtr = rows[i] - partitionOffsets[myid]; 
    csrRows[powPtr + 1]++;
  }
  for (int i = 0; i < nLocalRows; i++)
    csrRows[i + 1] += csrRows[i]; 

  if(useFP32) {
    float *csrValues = (float*) malloc(nnz * sizeof(float));
    for (int i = 0; i < nnz; i++) csrValues[i] = values[i]; 
    AMGX_matrix_upload_distributed(
      AmgXA, (int) nGlobalRows, nLocalRows, nnz, 1, 1,
      csrRows, cols, csrValues, NULL, dist);
    free(csrValues);
  } else {
    double *csrValues = (double*) malloc(nnz * sizeof(double));
    for (int i = 0; i < nnz; i++) csrValues[i] = values[i]; 
    AMGX_matrix_upload_distributed(
      AmgXA, (int) nGlobalRows, nLocalRows, nnz, 1, 1,
      csrRows, cols, csrValues, NULL, dist);
    free(csrValues);
  }
  AMGX_distribution_destroy(dist);
  free(csrRows);
  free(partitionOffsets);

  AMGX_solver_setup(solver, AmgXA);
  AMGX_vector_bind(AmgXP, AmgXA);
  AMGX_vector_bind(AmgXRHS, AmgXA);

  AMGXinit = 1;
}

int AMGX_t::solve(void *rhs, void *x)
{
  AMGX_vector_upload(AmgXP, nLocalRows, 1, x);
  AMGX_vector_upload(AmgXRHS, nLocalRows, 1, rhs);

  AMGX_solver_solve(solver, AmgXRHS, AmgXP);
  AMGX_SOLVE_STATUS status; 
  AMGX_solver_get_status(solver, &status);
  if (status != AMGX_SOLVE_SUCCESS) return status; 

  AMGX_vector_download(AmgXP, x);

  return 0;
}

void AMGXfinalize() 
{
  if(!AMGXinit) return;

  /* destroy config (need to use AMGX_SAFE_CALL after this point) */
  AMGX_SAFE_CALL(AMGX_finalize_plugins());
  AMGX_SAFE_CALL(AMGX_finalize());
}

AMGX_t::~AMGX_t()
{
  AMGX_vector_destroy(AmgXP);
  AMGX_vector_destroy(AmgXRHS);
  AMGX_matrix_destroy(AmgXA);
  AMGX_solver_destroy(solver);
  AMGX_resources_destroy(rsrc);
  AMGX_SAFE_CALL(AMGX_config_destroy(cfg));
  MPI_Comm_free(&comm);
}

int AMGXenabled()
{
  return 1;
}

#else
AMGX_t::AMGX_t(const int nLocalRows, const int nnz,
               const long long *rows, const long long *cols, const double *values, /* COO */
               const int nullspace, const MPI_Comm comm, int deviceID,
               int useFP32, int MPI_DIRECT, const char* cfgFile)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);  
  if(rank == 0) printf("ERROR: Recompile with AMGX support!\n");
}

int AMGX_t::solve(void *x, void *rhs)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);  
  if(rank == 0) printf("ERROR: Recompile with AMGX support!\n");
  return 1;
}

void AMGXfinalize()
{
}

AMGX_t::~AMGX_t()
{
}

int AMGXenabled()
{
  return 0;
}
#endif

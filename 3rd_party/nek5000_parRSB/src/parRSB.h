#ifndef _PARRSB_H_
#define _PARRSB_H_

#include "gslib.h"

#if !defined(MPI)
#error "gslib needs to be compiled with MPI"
#endif

#if !defined(GLOBAL_LONG_LONG)
#error "gslib needs to be compiled with GLOBAL_LONG_LONG"
#endif

#ifdef __cplusplus
extern "C" {
#endif

//==============================================================================
// Partitioning
//
typedef struct {
  // General options
  int partitioner;   // Partition algo: 0 - RSB, 1 - RCB, 2 - RIB (Default: 0)
  int verbose_level; // Verbose level: 0, 1, 2, .. etc (Default: 1)
  int profile_level; // Profile level: 0, 1, 2, .. etc (Default: 1)
  int two_level;     // Enable two level partitioning (Default: 0)
  int repair; // Repair disconnected components: 0 - No, 1 - Yes (Default: 0)
  // RSB common (Lanczos + MG) options
  int rsb_algo; // RSB algo: 0 - Lanczos, 1 - MG (Default: 0)
  int rsb_pre;  // RSB pre-partition : 0 - None, 1 - RCB , 2 - RIB (Default: 1)
  int rsb_max_iter;   // Max iterations in Lanczos / MG (Default: 50)
  int rsb_max_passes; // Max Lanczos restarts / Inverse iterations (Default: 50)
  double rsb_tol;     // Tolerance for Lanczos or RQI (Default: 1e-5)
  // RSB MG specific options
  int rsb_mg_grammian; // MG Grammian: 0 or 1 (Default: 0)
  int rsb_mg_factor;   // MG Coarsening factor (Default: 2, should be > 1)
  int rsb_mg_sagg;     // MG smooth aggregation: 0 or 1 (Default: 0)
} parrsb_options;

extern parrsb_options parrsb_default_options;

int parrsb_part_mesh(int *part, int *seq, long long *vtx, double *coord,
                     int nel, int nv, parrsb_options options, MPI_Comm comm);

#define fparrsb_part_mesh FORTRAN_UNPREFIXED(fparrsb_partmesh, FPARRSB_PARTMESH)
void fparrsb_part_mesh(int *part, int *seq, long long *vtx, double *coord,
                       int *nel, int *nve, int *options, int *comm, int *err);

//==============================================================================
// Connectivity
//
int parrsb_conn_mesh(long long *vtx, double *coord, int nel, int nDim,
                     long long *periodicInfo, int nPeriodicFaces, double tol,
                     MPI_Comm comm);

#define fparrsb_conn_mesh                                                      \
  FORTRAN_UNPREFIXED(fparrsb_conn_mesh, FPARRSB_CONN_MESH)
void fparrsb_conn_mesh(long long *vtx, double *coord, int *nel, int *nDim,
                       long long *periodicInfo, int *nPeriodicFaces,
                       double *tol, MPI_Fint *fcomm, int *err);

//==============================================================================
// I/O routines
//
int parrsb_read_mesh(unsigned *nel, unsigned *nv, long long **vl,
                     double **coord, unsigned *nbcs, long long **bcs,
                     char *name, MPI_Comm comm, int read);

int parrsb_dump_con(char *name, unsigned nelt, unsigned nv, long long *vl,
                    MPI_Comm comm);

int parrsb_dump_map(char *name, unsigned nelt, unsigned nv, long long *vl,
                    MPI_Comm comm);

int parrsb_dump_part(char *name, unsigned nelt, unsigned nv, double *coord,
                     int gid, MPI_Comm comm);

//==============================================================================
// Auxiliary functions
//
struct parrsb_cmd_opts {
  char *mesh;  // Mesh name, required.
  double tol;  // gencon tolerance, default: 0.2
  int test;    // run tests, default: 0
  int dump;    // dump the connectivity or map file, default: 1
  int nactive; // # of active MPI ranks, default: INT_MAX
  int verbose; // Verbosity, default: 0

  int ilu_type;   // ILU type, default: 0
  double ilu_tol; // ILU tolerance, default: 0.1
  int ilu_pivot;  // Pivoting for ILU: default: 0

  int crs_type;   // Coarse solver type, default: 0
  double crs_tol; // Coarse tolerance, default: 1e-3
};

struct parrsb_cmd_opts *parrsb_parse_cmd_opts(int argc, char *argv[]);
void parrsb_cmd_opts_free(struct parrsb_cmd_opts *opts);

int parrsb_dist_mesh(unsigned *nelt, long long **vl, double **coord, int *part,
                     int nv, MPI_Comm comm);

int parrsb_setup_mesh(unsigned *nelt, unsigned *nv, long long **vl,
                      double **coord, struct parrsb_cmd_opts *opts,
                      MPI_Comm comm);

void parrsb_print_part_stat(long long *vtx, unsigned nelt, unsigned nv,
                            MPI_Comm comm);

void parrsb_get_part_stat(int *nc, int *ns, int *nss, int *nel, long long *vtx,
                          int nelt, int nv, MPI_Comm comm);

void parrsb_check_error_(int err, char *file, int line, MPI_Comm comm);
#define parrsb_check_error(err, comm)                                          \
  parrsb_check_error_(err, __FILE__, __LINE__, comm)

#ifdef __cplusplus
}
#endif

#endif

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
  int partitioner; // Partition algo: 0 - RSB, 1 - RCB, 2 - RIB (Default: 0)
  int tagged;      // Tagged partitioning: 0 - No, 1 - Yes (Default: 0)
  int levels;      // Number of levels: 1, or 2 (Default: 2)
  int find_disconnected_comps; // Find number of components: 0 - No, 1 - Yes
                               // (Default: 1)
  int repair; // Repair disconnected components: 0 - No, 1 - Yes (Default: 0)
  int verbose_level; // Verbose level: 0, 1, 2, .. etc (Default: 1)
  int profile_level; // Profile level: 0, 1, 2, .. etc (Default: 0)
  // RSB common (Lanczos and MG) options
  int rsb_algo; // RSB algo: 0 - Lanczos, 1 - MG (Default: 0)
  int rsb_pre;  // RSB pre-partition : 0 - None, 1 - RCB , 2 - RIB (Default: 1)
  int rsb_max_iter;   // Max iterations in Lanczos / MG (Default: 50)
  int rsb_max_passes; // Max Lanczos restarts / Inverse iterations (Default: 50)
  double rsb_tol;     // Tolerance for Lanczos or RQI (Default: 1e-5)
  int rsb_dump_stats; // Dump partition statistics to a text file.
  // RSB MG specific options
  int rsb_mg_grammian; // MG Grammian: 0 or 1 (Default: 0)
  int rsb_mg_factor;   // MG Coarsening factor (Default: 2, should be > 1)
} parrsb_options;

extern parrsb_options parrsb_default_options;

int parrsb_part_mesh(int *part, const long long *const vtx,
                     const double *const xyz, const int *const tag,
                     const int nel, const int nv, parrsb_options *const options,
                     MPI_Comm comm);

void parrsb_part_solid(int *part, const long long *vtx2, unsigned nel2,
                       const long long *vtx1, unsigned nel1, unsigned nv,
                       MPI_Comm comm);

void parrsb_check_tagged_partitions(const long long *const eids,
                                    const long long *const vtx, const uint nel,
                                    const unsigned nv, const uint ntags,
                                    const struct comm *const c,
                                    const int verbose);
//==============================================================================
// Connectivity
//
int parrsb_conn_mesh(long long *vtx, double *coord, uint nel, unsigned nDim,
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

//==============================================================================
// Auxiliary functions
//
typedef struct {
  char *mesh;  // Mesh name, required.
  double tol;  // gencon tolerance, default: 0.2
  int test;    // run tests, default: 0
  int dump;    // dump the connectivity or map file, default: 1
  int nactive; // # of active MPI ranks, default: INT_MAX
  int verbose; // Verbosity, default: 0
} parrsb_cmd_line_opts;

parrsb_cmd_line_opts *parrsb_parse_cmd_opts(int argc, char *argv[]);

void parrsb_cmd_opts_free(parrsb_cmd_line_opts *opts);

int parrsb_dist_mesh(unsigned *nelt, long long **vl, double **coord, int *part,
                     unsigned nv, MPI_Comm comm);

int parrsb_setup_mesh(unsigned *nelt, unsigned *nv, long long **vl,
                      double **coord, parrsb_cmd_line_opts *opts,
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

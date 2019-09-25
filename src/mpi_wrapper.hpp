#if !defined(nekrs_mpi_wrapper_hpp_)
#define nekrs_mpi_wrapper_hpp_

#include "mpi.h"

namespace nekrs::mpi{
  static MPI_Comm comm_;
  static int root_;

  int mpiInit(MPI_Comm comm, int root=0);
  int mpiInit(int *argc,char **argv[],root=0);
  int Rank();
  int Size();
  int Barrier();
  double Wtime();
  MPI_Comm Comm();
  int mpiFinalize();
}

#endif

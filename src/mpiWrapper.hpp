#if !defined(nekrs_mpi_wrapper_hpp_)
#define nekrs_mpi_wrapper_hpp_

#include "mpi.h"

namespace nekrs::mpi{
  static MPI_Comm comm_;
  static int root_;

  int mpiInit(MPI_Comm comm, int root=0);
  int mpiInit(int *argc,char **argv[],int root=0);
  int Rank();
  int Size();
  int Barrier();
  double Wtime();
  MPI_Comm Comm();
  int Bcast(void *buffer,int size,int root);
  int Bcast(void *buffer,int size);
  int Allreduce(const void *sendbuf,void *recvbuf,int count,MPI_Op op);
  int mpiFinalize();
}

#endif

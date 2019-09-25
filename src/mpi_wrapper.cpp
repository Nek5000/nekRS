#include "mpi_wrapper.hpp"

int nekrs::mpi::mpiInit(MPI_Comm comm, int root){
  comm_=comm;
  root_=root;
}

int nekrs::mpi::mpiInit(int *argc,char **argv[],int root){
  int ret=MPI_Init(argc,argv);
  mpiInit(MPI_COMM_WORLD,root);
  return ret;
}

int nekrs::mpi::Rank(){
  int rank;
  MPI_Comm_rank(comm_,&rank);
  return rank;
}

int nekrs::mpi::Size(){
  int size;
  MPI_Comm_size(comm_,&size);
  return size;
}

int nekrs::mpi::Barrier(){
  MPI_Barrier(comm_);
}

double nekrs::mpi::Wtime(){
  return MPI_Wtime();
}

MPI_Comm nekrs::mpi::Comm(){
  return comm_;
}

int nekrs::mpi::mpiFinalize(){
  return MPI_Finalize();
}

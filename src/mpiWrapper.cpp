#include "mpiWrapper.hpp"

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

int nekrs::mpi::Bcast(void *buffer,int size,int root){
  return MPI_Bcast(buffer,size,MPI_BYTE,root,comm_);
}

int nekrs::mpi::Bcast(void *buffer,int size){
  return Bcast(buffer,size,root_);
}

int nekrs::mpi::Allreduce(const void *sendbuf,void *recvbuf,int count,MPI_Op op){
  MPI_Allreduce(sendbuf,recvbuf,count,MPI_BYTE,op,comm_);
}

int nekrs::mpi::mpiFinalize(){
  return MPI_Finalize();
}

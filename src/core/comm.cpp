#include "platform.hpp"
#include "comm.hpp"
#include "nrs.hpp"

comm_t::comm_t(MPI_Comm _commg, MPI_Comm _comm)
{

  mpiCommParent = _commg;
  mpiComm = _comm;
  MPI_Comm_rank(_comm, &mpiRank);
  MPI_Comm_size(_comm, &mpiCommSize);

  MPI_Comm_split_type(_comm, MPI_COMM_TYPE_SHARED, mpiRank, MPI_INFO_NULL, &mpiCommLocal);
  MPI_Comm_rank(mpiCommLocal, &localRank);
  MPI_Comm_size(mpiCommLocal, &mpiCommLocalSize);

  useGPUAware = true;
  const char * gpuMPIEnv = getenv("NEKRS_GPU_MPI");
  if(gpuMPIEnv){
    if(std::stoi(gpuMPIEnv)){
      useGPUAware = true;
    } else {
      useGPUAware = false;
    }
  }

}

MPI_Datatype comm_t::toMPI_Datatype(comm_t::type t) const
{
  switch(t)
  {
    case comm_t::type::dfloat:
      return MPI_DFLOAT;
    case comm_t::type::dlong:
      return MPI_DLONG;
    case comm_t::type::hlong:
      return MPI_HLONG;
  }
}

MPI_Op comm_t::toMPI_Op(comm_t::op o)const
{
  switch(o)
  {
    case comm_t::op::sum:
      return MPI_SUM;
    case comm_t::op::max:
      return MPI_MAX;
    case comm_t::op::min:
      return MPI_MIN;
  }
}

void comm_t::reallocScratch(size_t Nbytes) const
{
  if(useGPUAware) return;

  if(h_recvBuf.size() < Nbytes){
    if(h_recvBuf.size()) h_recvBuf.free();
    if(h_sendBuf.size()) h_sendBuf.free();
    h_recvBuf = platform->device.mallocHost(Nbytes);
    h_sendBuf = platform->device.mallocHost(Nbytes);
    recv = (void*) h_recvBuf.ptr();
    send = (void*) h_sendBuf.ptr();
  }

};

int comm_t::allreduce(const void *sendbuf, void *recvbuf, int count,
                comm_t::type datatype, comm_t::op op, MPI_Comm comm) const
{
  auto mpiDataType = toMPI_Datatype(datatype);
  auto mpiOp = toMPI_Op(op);

  return MPI_Allreduce(sendbuf, recvbuf, count, mpiDataType, mpiOp, comm);
}

int comm_t::allreduce(occa::memory sendbuf, occa::memory recvbuf, int count,
                comm_t::type datatype, comm_t::op op, MPI_Comm comm) const
{
  auto mpiDataType = toMPI_Datatype(datatype);
  auto mpiOp = toMPI_Op(op);

  int sizeBytes;
  MPI_Type_size(mpiDataType, &sizeBytes);

  const size_t Nbytes = sizeBytes * count;

  reallocScratch(Nbytes);

  if(useGPUAware || platform->serial){
    return MPI_Allreduce((void*) sendbuf.ptr(), (void*) recvbuf.ptr(), count, mpiDataType, mpiOp, comm);
  } else {
    int retVal = 0;

    sendbuf.copyTo(send, Nbytes);
    retVal = MPI_Allreduce(send, recv, count, mpiDataType, mpiOp, comm);
    recvbuf.copyFrom(recv, Nbytes);

    return retVal;
  }
}

// in place
int comm_t::allreduce(occa::memory recvbuf, int count,
                comm_t::type datatype, comm_t::op op, MPI_Comm comm) const
{
  auto mpiDataType = toMPI_Datatype(datatype);
  auto mpiOp = toMPI_Op(op);

  int sizeBytes;
  MPI_Type_size(mpiDataType, &sizeBytes);

  const size_t Nbytes = sizeBytes * count;

  reallocScratch(Nbytes);

  if(useGPUAware || platform->serial){
    return MPI_Allreduce(MPI_IN_PLACE, (void*) recvbuf.ptr(), count, mpiDataType, mpiOp, comm);
  } else {
    int retVal = 0;

    recvbuf.copyTo(recv, Nbytes);
    retVal = MPI_Allreduce(MPI_IN_PLACE, recv, count, mpiDataType, mpiOp, comm);
    recvbuf.copyFrom(recv, Nbytes);

    return retVal;
  }
}

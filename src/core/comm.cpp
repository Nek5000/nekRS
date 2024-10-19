#include "platform.hpp"
#include "comm.hpp"

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
  const char *gpuMPIEnv = getenv("NEKRS_GPU_MPI");
  if (gpuMPIEnv) {
    if (std::stoi(gpuMPIEnv)) {
      useGPUAware = true;
    } else {
      useGPUAware = false;
    }
  }
}

MPI_Datatype comm_t::toMPI_Datatype(const occa::memory& t) const
{
  const auto type = t.dtype(); 

  if (type == occa::dtype::get<double>()) 
   return MPI_DOUBLE;
  else if (type == occa::dtype::get<float>()) 
   return MPI_FLOAT;
  else if (type == occa::dtype::get<int>()) 
   return MPI_INT;
  else if (type == occa::dtype::get<long long int>()) 
   return MPI_LONG_LONG_INT;
  else
    nekrsAbort(MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "Unkown datatype!");

  return 0;
}

MPI_Op comm_t::toMPI_Op(comm_t::op o) const
{
  switch (o) {
  case comm_t::op::sum:
    return MPI_SUM;
  case comm_t::op::max:
    return MPI_MAX;
  case comm_t::op::min:
    return MPI_MIN;
  default:
    nekrsAbort(MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "Unknown operation!");
  }

  return 0;
}

void comm_t::reallocScratch(size_t Nbytes) const
{
  if (useGPUAware) {
    return;
  }

  if (h_recvBuf.size() < Nbytes) {
    if (h_recvBuf.size()) {
      h_recvBuf.free();
    }
    if (h_sendBuf.size()) {
      h_sendBuf.free();
    }
    h_recvBuf = platform->device.mallocHost(Nbytes);
    h_sendBuf = platform->device.mallocHost(Nbytes);
    recv = (void *)h_recvBuf.ptr();
    send = (void *)h_sendBuf.ptr();
  }
};

int comm_t::allreduce(occa::memory recvbuf,
                      int count,
                      comm_t::op op,
                      MPI_Comm comm) const
{
  auto mpiDataType = toMPI_Datatype(recvbuf);
  auto mpiOp = toMPI_Op(op);

  int sizeBytes;
  MPI_Type_size(mpiDataType, &sizeBytes);

  const size_t Nbytes = sizeBytes * count;

  reallocScratch(Nbytes);

  if (useGPUAware || platform->serial) {
    platform->device.finish();
    return MPI_Allreduce(MPI_IN_PLACE, (void *)recvbuf.ptr(), count, mpiDataType, mpiOp, comm);
  } else {
    recvbuf.copyTo(recv, count);
    int retVal = MPI_Allreduce(MPI_IN_PLACE, recv, count, mpiDataType, mpiOp, comm);
    recvbuf.copyFrom(recv, count);
    return retVal;
  }
}

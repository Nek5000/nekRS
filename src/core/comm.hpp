#ifndef comm_hpp_
#define comm_hpp_
class comm_t{
public:
  comm_t(MPI_Comm, MPI_Comm);
  MPI_Comm mpiCommParent;
  MPI_Comm mpiComm;
  int mpiRank;
  int mpiCommSize;

  MPI_Comm mpiCommLocal;
  int mpiCommLocalSize;
  int localRank;

  enum class type{
    dfloat,
    dlong,
    hlong,
  };

  enum class op{
    sum,
    max,
    min,
  };

  std::string to_string() const {
    std::ostringstream ss;
    ss << "mpiRank = " << mpiRank << std::endl;
    ss << "mpiCommSize = " << mpiCommSize << std::endl;
    ss << "mpiCommLocalSize = " << mpiCommLocalSize << std::endl;
    ss << "localRank = " << localRank << std::endl;
    return ss.str();
  }

  int allreduce(const void *sendbuf, void *recvbuf, int count,
                  type datatype, op op, MPI_Comm comm) const;
  int allreduce(occa::memory sendbuf, occa::memory recvbuf, int count,
                  type datatype, op op, MPI_Comm comm) const;
  
  // in place
  int allreduce(occa::memory recvbuf, int count,
                  type datatype, op op, MPI_Comm comm) const;
  
private:

  MPI_Datatype toMPI_Datatype(type t) const;
  MPI_Op toMPI_Op(op t) const;

  void reallocScratch(size_t Nbytes) const;
  bool useGPUAware;

  mutable occa::memory h_recvBuf;
  mutable occa::memory h_sendBuf;
  mutable void* recv;
  mutable void* send;

};

#endif
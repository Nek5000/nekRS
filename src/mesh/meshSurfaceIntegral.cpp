#include <mesh.h>
#include "platform.hpp"

namespace {
  dfloat *sum;
  dfloat *sumFace;
  occa::memory o_sumFace;
  occa::memory h_sumFace;
}

std::vector<dfloat> mesh_t::surfaceIntegral(int nbID, const occa::memory& o_bID, const occa::memory& o_fld)
{
  return surfaceIntegral(1, 0, nbID, o_bID, o_fld);
}

std::vector<dfloat> mesh_t::surfaceIntegral(int Nfields, int fieldOffset, int nbID,
                                            const occa::memory o_bID, const occa::memory& o_fld)
{
  const auto Nbytes = (Nfields * sizeof(dfloat)) * Nelements;

  if (o_sumFace.size() < Nbytes) {
    if (o_sumFace.size()) o_sumFace.free();
    o_sumFace = platform->device.malloc(Nbytes);

    if (h_sumFace.size()) h_sumFace.free();
    h_sumFace = platform->device.mallocHost(Nbytes);
    sumFace = (dfloat *) h_sumFace.ptr();

    if (sum) free(sum);
    sum = (dfloat *) std::malloc(Nbytes);
  }

  surfaceIntegralKernel(Nelements, 
                        Nfields,
                        fieldOffset, 
                        nbID,
                        o_bID,
                        o_sgeo, 
                        o_vmapM, 
                        o_EToB,  
                        o_fld,
                        o_sumFace);

  o_sumFace.copyTo(sumFace, Nbytes);

  for (int j = 0; j < Nfields + 1; ++j) {
    sum[j] = 0;
    for (int i = 0; i < Nelements; ++i) {
      sum[j] += sumFace[i + j * Nelements];
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, sum, Nfields, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);

  std::vector<dfloat> out; 
  for (int i = 0; i < Nfields; ++i)
    out.push_back(sum[i]);
  
  return out;
}

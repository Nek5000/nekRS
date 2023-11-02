#include <mesh.h>
#include "platform.hpp"

namespace {

dfloat *sum;
dfloat *sumFace;
occa::memory o_sumFace;
occa::memory h_sumFace;

std::vector<dfloat> integral(mesh_t *mesh, int Nfields, int fieldOffset, bool vector, int nbID,
                             const occa::memory o_bID, const occa::memory& o_fld)
{
  if (o_sumFace.length() < Nfields * mesh->Nelements) {
    if (o_sumFace.size()) o_sumFace.free();
    o_sumFace = platform->device.malloc<dfloat>(Nfields * mesh->Nelements);
    if (h_sumFace.length()) h_sumFace.free();
    h_sumFace = platform->device.mallocHost<dfloat>(Nfields * mesh->Nelements);
    sumFace = (dfloat *) h_sumFace.ptr();

    if (sum) free(sum);
    sum = (dfloat *) calloc(Nfields * mesh->Nelements, sizeof(dfloat));
  }

  if (vector)
    mesh->surfaceIntegralVectorKernel(mesh->Nelements, 
                                Nfields,
                                fieldOffset, 
                                nbID,
                                o_bID,
                                mesh->o_sgeo, 
                                mesh->o_vmapM, 
                                mesh->o_EToB,  
                                o_fld,
                                o_sumFace);
  else 
    mesh->surfaceIntegralKernel(mesh->Nelements, 
                          Nfields,
                          fieldOffset, 
                          nbID,
                          o_bID,
                          mesh->o_sgeo, 
                          mesh->o_vmapM, 
                          mesh->o_EToB,  
                          o_fld,
                          o_sumFace);

  o_sumFace.copyTo(sumFace, Nfields * mesh->Nelements);

  for (int j = 0; j < Nfields; ++j) {
    sum[j] = 0;
    for (int i = 0; i < mesh->Nelements; ++i) {
      sum[j] += sumFace[i + j * mesh->Nelements];
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, sum, Nfields, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);

  std::vector<dfloat> out; 
  for (int i = 0; i < Nfields; ++i)
    out.push_back(sum[i]);
  
  return out;
}

} // private namespace

std::vector<dfloat> mesh_t::surfaceIntegral(int nbID, const occa::memory& o_bID, const occa::memory& o_fld)
{
  return integral(this, 1, static_cast<dlong>(0), false, nbID, o_bID, o_fld);
}

std::vector<dfloat> mesh_t::surfaceIntegralVector(dlong fieldOffset, int nbID, const occa::memory& o_bID, const occa::memory& o_fld)
{
  return integral(this, 1, fieldOffset, true, nbID, o_bID, o_fld);
}

std::vector<dfloat> mesh_t::surfaceIntegralMany(int Nfields, dlong fieldOffset, int nbID,
                                                const occa::memory& o_bID, const occa::memory& o_fld)
{
  return integral(this, Nfields, fieldOffset, false, nbID, o_bID, o_fld);
}

#include <mesh.h>
#include "platform.hpp"

static dfloat *sum;
static dfloat *sumFace;
static occa::memory o_sumFace;
static occa::memory h_sumFace;

dfloat mesh_t::avgBoundaryValue(int BID, occa::memory o_fld)
{
  dfloat avg = 0.0;
  avgBoundaryValue(BID, 1, fieldOffset, o_fld, &avg);
  return avg;
}

void mesh_t::avgBoundaryValue(int BID, int Nfields, int offsetFld, occa::memory o_fld, dfloat *avgs)
{
  const auto offset = Nfaces * Nelements;
  const auto Nbytes = ((Nfields + 1) * sizeof(dfloat)) * offset;

  if (o_sumFace.size() < Nbytes) {
    if (o_sumFace.size())
      o_sumFace.free();
    if (h_sumFace.size())
      h_sumFace.free();

    // pinned scratch buffer
    {
      h_sumFace = platform->device.mallocHost(Nbytes);
      sumFace = (dfloat *)h_sumFace.ptr();
    }

    o_sumFace = platform->device.malloc(Nbytes);

    if (sum)
      free(sum);
    sum = (dfloat *)calloc(Nfields + 1, sizeof(dfloat));
  }

  avgBIDValueKernel(Nelements, BID, Nfields, offsetFld, offset, o_sgeo, o_EToB, o_vmapM, o_fld, o_sumFace);

  o_sumFace.copyTo(sumFace, Nbytes);

  for (int j = 0; j < Nfields + 1; ++j) {
    sum[j] = 0;
    for (int i = 0; i < offset; ++i) {
      sum[j] += sumFace[i + j * offset];
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, sum, Nfields + 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);

  const auto invArea = 1 / sum[Nfields];
  for (int i = 0; i < Nfields; ++i)
    avgs[i] = sum[i] * invArea;
}
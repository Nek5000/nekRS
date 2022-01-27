#include <mesh.h>
#include "platform.hpp"

dfloat mesh_t::avgBoundaryValue(int BID, occa::memory fld)
{

  if (!o_sum.isInitialized()) {
    o_sum = platform->device.malloc(Nfaces * Nelements * sizeof(dfloat));
    sum = (dfloat *)calloc(Nfaces * Nelements, sizeof(dfloat));
  }
  if (!o_area.isInitialized()) {
    o_area = platform->device.malloc(Nfaces * Nelements * sizeof(dfloat));
    area = (dfloat *)calloc(Nfaces * Nelements, sizeof(dfloat));
  }

  avgBIDValueKernel(Nelements, BID, o_sgeo, o_EToB, o_vmapM, fld, o_sum, o_area);

  o_sum.copyTo(sum, Nfaces * Nelements * sizeof(dfloat));
  o_area.copyTo(area, Nfaces * Nelements * sizeof(dfloat));

  dfloat localSum = 0.0;
  dfloat localArea = 0.0;
  for (int face = 0; face < Nfaces * Nelements; ++face) {
    localSum += sum[face];
    localArea += area[face];
  }

  dfloat values[2] = {localSum, localArea};
  MPI_Allreduce(MPI_IN_PLACE, values, 2, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);

  return values[0] / values[1];
}
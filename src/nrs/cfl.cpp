#include "nrs.hpp"

namespace
{
bool firstTime = true;
static occa::memory h_scratch;
static occa::memory o_dx;

void setup(nrs_t *nrs)
{
  auto mesh = nrs->mesh;
  h_scratch = platform->device.mallocHost<dfloat>(mesh->Nelements);

  if (nrs->elementType == QUADRILATERALS || nrs->elementType == HEXAHEDRA) {
    std::vector<dfloat> dx(mesh->N + 1);

    for (int n = 0; n < (mesh->N + 1); n++) {
      if (n == 0) {
        dx[n] = mesh->gllz[n + 1] - mesh->gllz[n];
      } else if (n == mesh->N) {
        dx[n] = mesh->gllz[n] - mesh->gllz[n - 1];
      } else {
        dx[n] = 0.5 * (mesh->gllz[n + 1] - mesh->gllz[n - 1]);
      }

      dx[n] = 1.0 / dx[n];
    }

    o_dx = platform->device.malloc<dfloat>(mesh->N + 1);
    o_dx.copyFrom(dx.data());
  }
  firstTime = false;
}

} // namespace

dfloat nrs_t::computeCFL()
{
  return computeCFL(this->dt[0]);
}

dfloat nrs_t::computeCFL(dfloat dt)
{
  if (firstTime) {
    setup(this);
  }

  auto o_cfl = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nelements);

  this->cflKernel(mesh->Nelements, dt, mesh->o_vgeo, o_dx, this->fieldOffset, this->o_U, mesh->o_U, o_cfl);

  auto scratch = (dfloat *)h_scratch.ptr();
  o_cfl.copyTo(scratch);

  dfloat cfl = 0;
  for (dlong n = 0; n < mesh->Nelements; ++n) {
    cfl = std::max(cfl, scratch[n]);
  }

  MPI_Allreduce(MPI_IN_PLACE, &cfl, 1, MPI_DFLOAT, MPI_MAX, platform->comm.mpiComm);

  return cfl;
}

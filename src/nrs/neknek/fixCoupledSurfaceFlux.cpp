#include "bcMap.hpp"
#include "neknek.hpp"
#include "nrs.hpp"
#include <array>

static bool findOutlet(mesh_t *mesh)
{
  dlong numOutlet = 0;
  for (dlong e = 0; e < mesh->Nelements; e++) {
    for (dlong f = 0; f < mesh->Nfaces; f++) {
      auto bID = mesh->EToB[f + mesh->Nfaces * e];
      auto bcType = bcMap::id(bID, "velocity");
      if (bcType == bcMap::bcTypeONX || bcType == bcMap::bcTypeONY || bcType == bcMap::bcTypeONZ ||
          bcType == bcMap::bcTypeON || bcType == bcMap::bcTypeO) {
        numOutlet++;
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &numOutlet, 1, MPI_DLONG, MPI_SUM, platform->comm.mpiComm);

  return (numOutlet == 0) ? false : true;
}

void neknek_t::fixCoupledSurfaceFlux(occa::memory o_U)
{
  if (!nrs->flow) {
    return;
  }

  auto mesh = nrs->mesh;

  static bool isCalled = false;
  static bool hasOutlet;
  if (!isCalled) {
    hasOutlet = findOutlet(mesh);
    isCalled = true;
  }
  if (hasOutlet) {
    return;
  }

  constexpr int nReduction = 2;
  auto o_reduction = platform->deviceMemoryPool.reserve<dfloat>(nReduction * mesh->Nelements);

  this->computeFluxKernel(mesh->Nelements,
                          nrs->fieldOffset,
                          mesh->o_sgeo,
                          mesh->o_vmapM,
                          nrs->o_EToB,
                          o_U,
                          o_reduction);

  std::vector<dfloat> reduction(nReduction * mesh->Nelements);
  o_reduction.copyTo(reduction.data());

  std::array<dfloat, nReduction> res;
  for (int fld = 0; fld < nReduction; fld++) {
    res[fld] = 0;
    for (int e = 0; e < mesh->Nelements; e++) {
      res[fld] += reduction[e + fld * mesh->Nelements];
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, res.data(), nReduction, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);

  auto [flux, area] = res;

  dfloat gamma = 0.0;
  if (area > 0.0) {
    gamma = -1.0 * flux / area;
  }

  if (platform->verbose && platform->comm.mpiRank == 0) {
    printf("neknek::fixCoupledSurfaceFlux flux = %11.4e, area = %11.4e, gamma = %11.4e\n", flux, area, gamma);
  }

  this->fixSurfaceFluxKernel(mesh->Nelements,
                             nrs->fieldOffset,
                             mesh->o_sgeo,
                             mesh->o_vmapM,
                             nrs->o_EToB,
                             gamma,
                             o_U);
}

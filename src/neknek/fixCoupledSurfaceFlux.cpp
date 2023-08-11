#include "neknek.hpp"
#include "nrs.hpp"
#include <array>

void fixCoupledSurfaceFlux(nrs_t *nrs, occa::memory o_U)
{
  if (nrs->neknek == nullptr || !nrs->flow) {
    return;
  }

  auto neknek = nrs->neknek;
  auto mesh = nrs->meshV;

  constexpr int nReduction = 2;
  auto o_reduction = platform->o_memPool.reserve<dfloat>(nReduction * mesh->Nelements);

  neknek->computeFluxKernel(mesh->Nelements,
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

  neknek->fixSurfaceFluxKernel(mesh->Nelements,
                               nrs->fieldOffset,
                               mesh->o_sgeo,
                               mesh->o_vmapM,
                               nrs->o_EToB,
                               gamma,
                               o_U);
}
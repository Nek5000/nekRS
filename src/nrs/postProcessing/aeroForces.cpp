#include "nrs.hpp"
#include "aeroForce.hpp"

namespace
{
static std::vector<dfloat> tmp;
} // namespace

AeroForce *nrs_t::aeroForces(int nbID, const occa::memory &o_bID, const occa::memory &o_Sij_)
{
  occa::memory o_Sij = o_Sij_;
  if (!o_Sij.isInitialized()) {
    o_Sij = this->strainRate();
  }

  auto af = new AeroForce();

  auto o_rho = af->rho();
  if (!o_rho.isInitialized()) {
    o_rho = this->o_rho;
  }

  auto o_forces = platform->deviceMemoryPool.reserve<dfloat>(2 * mesh->dim * mesh->Nelements);
  static occa::kernel kernel;
  if (!kernel.isInitialized()) {
    kernel = platform->kernelRequests.load("nrs-aeroForces");
  }
  kernel(mesh->Nelements,
         this->fieldOffset,
         nbID,
         o_bID,
         mesh->o_sgeo,
         mesh->o_vmapM,
         mesh->o_EToB,
         o_rho,
         this->o_mue,
         this->o_P,
         o_Sij,
         o_forces);

  if (tmp.size() < o_forces.size()) {
    tmp.resize(o_forces.size());
  }
  o_forces.copyTo(tmp.data());

  std::vector<dfloat> sum{0, 0, 0, 0, 0, 0};
  for (dlong i = 0; i < mesh->Nelements; i++) {
    sum[0] += tmp[i + 0 * mesh->Nelements];
    sum[1] += tmp[i + 1 * mesh->Nelements];
    sum[2] += tmp[i + 2 * mesh->Nelements];

    sum[3] += tmp[i + 3 * mesh->Nelements];
    sum[4] += tmp[i + 4 * mesh->Nelements];
    sum[5] += tmp[i + 5 * mesh->Nelements];
  }

  MPI_Allreduce(MPI_IN_PLACE, sum.data(), sum.size(), MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);

  af->forceViscous({sum[0], sum[1], sum[2]});
  af->forcePressure({sum[3], sum[4], sum[5]});

  return af;
}

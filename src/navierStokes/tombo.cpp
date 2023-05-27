#include "nrs.hpp"
#include "udf.hpp"
#include "linAlg.hpp"
#include "neknek.hpp"
#include "bdry.hpp"
#include "bcMap.hpp"
#include <limits>

namespace tombo
{
occa::memory pressureSolve(nrs_t* nrs, dfloat time, int stage)
{
  platform->timer.tic("pressure rhs", 1);
  double flopCount = 0.0;
  mesh_t* mesh = nrs->meshV;

  nrs->curlKernel(mesh->Nelements,
	                1,
                  mesh->o_vgeo,
                  mesh->o_D,
                  nrs->fieldOffset,
                  nrs->o_Ue,
                  platform->o_mempool.slice0);
  flopCount += static_cast<double>(mesh->Nelements) * (18 * mesh->Np * mesh->Nq + 36 * mesh->Np);

  oogs::startFinish(platform->o_mempool.slice0, nrs->NVfields, nrs->fieldOffset,ogsDfloat, ogsAdd, nrs->gsh);

  platform->linAlg->axmyVector(
    mesh->Nlocal,
    nrs->fieldOffset,
    0,
    1.0,
    nrs->meshV->o_invLMM,
    platform->o_mempool.slice0
  );
  flopCount += mesh->Nlocal;

  nrs->curlKernel(
    mesh->Nelements,
    1,
    mesh->o_vgeo,
    mesh->o_D,
    nrs->fieldOffset,
    platform->o_mempool.slice0,
    platform->o_mempool.slice3);
  flopCount += static_cast<double>(mesh->Nelements) * (18 * mesh->Np * mesh->Nq + 36 * mesh->Np);

  nrs->gradientVolumeKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    mesh->o_D,
    nrs->fieldOffset,
    nrs->o_div,
    platform->o_mempool.slice0);
  flopCount += static_cast<double>(mesh->Nelements) * (6 * mesh->Np * mesh->Nq + 18 * mesh->Np);

  if (platform->options.compareArgs("VELOCITY STRESSFORMULATION", "TRUE")) {
    nrs->pressureStressKernel(
         mesh->Nelements,
         mesh->o_vgeo,
         mesh->o_D,
         nrs->fieldOffset,
         nrs->o_mue,
         nrs->o_Ue,
         nrs->o_div,
         platform->o_mempool.slice3);
    flopCount += static_cast<double>(mesh->Nelements) * (18 * mesh->Nq * mesh->Np + 100 * mesh->Np);
  }

  const int viscContribution = (nrs->nBDF > 1) ? 1 : 0; 
  occa::memory o_irho = nrs->o_ellipticCoeff;
  nrs->pressureRhsKernel(
    mesh->Nelements * mesh->Np,
    nrs->fieldOffset,
    viscContribution,
    nrs->o_mue,
    o_irho,
    nrs->o_BF,
    platform->o_mempool.slice3,
    platform->o_mempool.slice0,
    platform->o_mempool.slice6);
  flopCount += 12 * static_cast<double>(mesh->Nlocal);

  oogs::startFinish(platform->o_mempool.slice6, nrs->NVfields, nrs->fieldOffset,ogsDfloat, ogsAdd, nrs->gsh);

  platform->linAlg->axmyVector(
    mesh->Nlocal,
    nrs->fieldOffset,
    0,
    1.0,
    nrs->meshV->o_invLMM,
    platform->o_mempool.slice6
  );

  nrs->wDivergenceVolumeKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    mesh->o_D,
    nrs->fieldOffset,
    platform->o_mempool.slice6,
    platform->o_mempool.slice3);
  flopCount += static_cast<double>(mesh->Nelements) * (6 * mesh->Np * mesh->Nq + 18 * mesh->Np);

  nrs->pressureAddQtlKernel(
    mesh->Nlocal,
    mesh->o_LMM,
    nrs->g0 * nrs->idt,
    nrs->o_div,
    platform->o_mempool.slice3);
  flopCount += 3 * mesh->Nlocal;

  nrs->divergenceSurfaceKernel(
    mesh->Nelements,
    mesh->o_sgeo,
    mesh->o_vmapM,
    nrs->o_EToB,
    nrs->g0 * nrs->idt,
    nrs->fieldOffset,
    platform->o_mempool.slice6,
    nrs->o_U,
    platform->o_mempool.slice3);
  flopCount += 25 * static_cast<double>(mesh->Nelements) * mesh->Nq * mesh->Nq;

  platform->timer.toc("pressure rhs");

  platform->o_mempool.slice1.copyFrom(nrs->o_P, mesh->Nlocal * sizeof(dfloat));
  ellipticSolve(nrs->pSolver, platform->o_mempool.slice3, platform->o_mempool.slice1);

  platform->flopCounter->add("pressure RHS", flopCount);

  return platform->o_mempool.slice1;
}

occa::memory velocitySolve(nrs_t* nrs, dfloat time, int stage)
{
  platform->timer.tic("velocity rhs", 1);
  double flopCount = 0.0;
  mesh_t* mesh = nrs->meshV;

  platform->linAlg->axmyz(mesh->Nlocal,
                          (platform->options.compareArgs("VELOCITY STRESSFORMULATION", "TRUE")) ? -2. / 3
                                                                                                : 1. / 3,
                          nrs->o_mue,
                          nrs->o_div,
                          platform->o_mempool.slice3);
  nrs->gradientVolumeKernel(mesh->Nelements,
                            mesh->o_vgeo,
                            mesh->o_D,
                            nrs->fieldOffset,
                            platform->o_mempool.slice3,
                            platform->o_mempool.slice0);
  flopCount += static_cast<double>(mesh->Nelements) * (6 * mesh->Np * mesh->Nq + 18 * mesh->Np);

  bool weakPressure = true;

  if (weakPressure) {
    nrs->wgradientVolumeKernel(mesh->Nelements,
                               mesh->o_vgeo,
                               mesh->o_D,
                               nrs->fieldOffset,
                               nrs->o_P,
                               platform->o_mempool.slice3);

    platform->linAlg->axpby(nrs->NVfields * nrs->fieldOffset,
                            1.0,
                            platform->o_mempool.slice3,
                            1.0,
                            platform->o_mempool.slice0);
  }
  else {
    nrs->gradientVolumeKernel(mesh->Nelements,
                              mesh->o_vgeo,
                              mesh->o_D,
                              nrs->fieldOffset,
                              nrs->o_P,
                              platform->o_mempool.slice3);

    platform->linAlg->axpby(nrs->NVfields * nrs->fieldOffset,
                            -1.0,
                            platform->o_mempool.slice3,
                            1.0,
                            platform->o_mempool.slice0);
  }
  flopCount += static_cast<double>(mesh->Nelements) * 18 * (mesh->Np * mesh->Nq + mesh->Np);

  nrs->velocityNeumannBCKernel(mesh->Nelements,
                               nrs->fieldOffset,
                               mesh->o_sgeo,
                               mesh->o_vmapM,
                               mesh->o_EToB,
                               nrs->o_EToB,
                               time,
                               mesh->o_x,
                               mesh->o_y,
                               mesh->o_z,
                               nrs->o_rho,
                               nrs->o_mue,
                               nrs->o_usrwrk,
                               nrs->o_Ue,
                               platform->o_mempool.slice0);

  flopCount += static_cast<double>(mesh->Nelements) * (3 * mesh->Np + 36 * mesh->Nq * mesh->Nq);

  nrs->velocityRhsKernel(
    mesh->Nlocal,
    nrs->fieldOffset,
    nrs->o_BF,
    platform->o_mempool.slice0,
    nrs->o_rho,
    platform->o_mempool.slice3);

  flopCount += 6 * mesh->Nlocal;

  platform->timer.toc("velocity rhs");
  platform->o_mempool.slice0.copyFrom(nrs->o_U, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));

  occa::memory o_U0;
  o_U0 = platform->options.compareArgs("VELOCITY INITIAL GUESS", "EXTRAPOLATION") && stage == 1 ? nrs->o_Ue
                                                                                                : nrs->o_U;

  platform->o_mempool.slice0.copyFrom(o_U0, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));

  if(nrs->uvwSolver) {
    ellipticSolve(nrs->uvwSolver, platform->o_mempool.slice3, platform->o_mempool.slice0);
  } else {
    ellipticSolve(nrs->uSolver, platform->o_mempool.slice3, platform->o_mempool.slice0);
    ellipticSolve(nrs->vSolver, platform->o_mempool.slice4, platform->o_mempool.slice1);
    ellipticSolve(nrs->wSolver, platform->o_mempool.slice5, platform->o_mempool.slice2);
  }

  platform->flopCounter->add("velocity RHS", flopCount);

  return platform->o_mempool.slice0;
}

} // namespace

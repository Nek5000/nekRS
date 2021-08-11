#include "nrs.hpp"
#include "udf.hpp"
#include "linAlg.hpp"
#include "neknek.hpp"

namespace tombo
{
occa::memory pressureSolve(nrs_t* nrs, dfloat time, int stage)
{
  mesh_t* mesh = nrs->meshV;


  //enforce Dirichlet BCs
  platform->linAlg->fill((1+nrs->NVfields)*nrs->fieldOffset, std::numeric_limits<dfloat>::min(), platform->o_mempool.slice6);
  for (int sweep = 0; sweep < 2; sweep++) {
    nrs->pressureDirichletBCKernel(mesh->Nelements,
                                   time,
                                   nrs->fieldOffset,
                                   mesh->o_sgeo,
                                   mesh->o_x,
                                   mesh->o_y,
                                   mesh->o_z,
                                   mesh->o_vmapM,
                                   mesh->o_EToB,
                                   nrs->o_EToB,
                                   nrs->o_usrwrk,
                                   nrs->o_U,
                                   nrs->o_P,
                                   platform->o_mempool.slice6);

    nrs->velocityDirichletBCKernel(mesh->Nelements,
                                   nrs->fieldOffset,
                                   time,
                                   mesh->o_sgeo,
                                   mesh->o_x,
                                   mesh->o_y,
                                   mesh->o_z,
                                   mesh->o_vmapM,
                                   mesh->o_EToB,
                                   nrs->o_EToB,
                                   nrs->o_VmapB,
                                   nrs->neknek->o_pointMap,
                                   nrs->neknek->o_valInterp,
                                   nrs->o_usrwrk,
                                   nrs->o_U,
                                   platform->o_mempool.slice7);

    //take care of Neumann-Dirichlet shared edges across elements
    if (sweep == 0) oogs::startFinish(platform->o_mempool.slice6, 1+nrs->NVfields, nrs->fieldOffset, ogsDfloat, ogsMax, nrs->gsh);
    if (sweep == 1) oogs::startFinish(platform->o_mempool.slice6, 1+nrs->NVfields, nrs->fieldOffset, ogsDfloat, ogsMin, nrs->gsh);
  }

  if (nrs->pSolver->Nmasked) nrs->maskCopyKernel(nrs->pSolver->Nmasked, 0, nrs->pSolver->o_maskIds,
                                                 platform->o_mempool.slice6, nrs->o_P);

  if (nrs->uvwSolver) {
    if (nrs->uvwSolver->Nmasked) nrs->maskCopyKernel(nrs->uvwSolver->Nmasked, 0*nrs->fieldOffset, nrs->uvwSolver->o_maskIds,
                                                     platform->o_mempool.slice7, nrs->o_U);
  } else {
    if (nrs->uSolver->Nmasked) nrs->maskCopyKernel(nrs->uSolver->Nmasked, 0*nrs->fieldOffset, nrs->uSolver->o_maskIds,
                                                   platform->o_mempool.slice7, nrs->o_U);
    if (nrs->vSolver->Nmasked) nrs->maskCopyKernel(nrs->vSolver->Nmasked, 1*nrs->fieldOffset, nrs->vSolver->o_maskIds,
                                                   platform->o_mempool.slice7, nrs->o_U);
    if (nrs->wSolver->Nmasked) nrs->maskCopyKernel(nrs->wSolver->Nmasked, 2*nrs->fieldOffset, nrs->wSolver->o_maskIds,
                                                   platform->o_mempool.slice7, nrs->o_U);
  }

  nrs->curlKernel(mesh->Nelements,
                  mesh->o_vgeo,
                  mesh->o_D,
                  nrs->fieldOffset,
                  nrs->o_Ue,
                  platform->o_mempool.slice0);

  oogs::startFinish(platform->o_mempool.slice0, nrs->NVfields, nrs->fieldOffset,ogsDfloat, ogsAdd, nrs->gsh);

  platform->linAlg->axmyVector(
    mesh->Nlocal,
    nrs->fieldOffset,
    0,
    1.0,
    nrs->meshV->o_invLMM,
    platform->o_mempool.slice0
  );

  nrs->curlKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    mesh->o_D,
    nrs->fieldOffset,
    platform->o_mempool.slice0,
    platform->o_mempool.slice3);

  nrs->gradientVolumeKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    mesh->o_D,
    nrs->fieldOffset,
    nrs->o_div,
    platform->o_mempool.slice0);

  if(platform->options.compareArgs("STRESSFORMULATION", "TRUE"))
    nrs->pressureStressKernel(
         mesh->Nelements,
         mesh->o_vgeo,
         mesh->o_D,
         nrs->fieldOffset,
         nrs->o_mue,
         nrs->o_Ue,
         nrs->o_div,
         platform->o_mempool.slice3);

  occa::memory o_irho = nrs->o_ellipticCoeff;
  nrs->pressureRhsKernel(
    mesh->Nelements * mesh->Np,
    nrs->fieldOffset,
    nrs->o_mue,
    o_irho,
    nrs->o_BF,
    platform->o_mempool.slice3,
    platform->o_mempool.slice0,
    platform->o_mempool.slice6);


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


  nrs->pressureAddQtlKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    nrs->g0 * nrs->idt,
    nrs->o_div,
    platform->o_mempool.slice3);

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

  platform->o_mempool.slice1.copyFrom(nrs->o_P, mesh->Nlocal * sizeof(dfloat));
  ellipticSolve(nrs->pSolver, platform->o_mempool.slice3, platform->o_mempool.slice1);

  return platform->o_mempool.slice1;
}

occa::memory velocitySolve(nrs_t* nrs, dfloat time, int stage)
{
  mesh_t* mesh = nrs->meshV;

  dfloat scale = -1./3;
  if(platform->options.compareArgs("STRESSFORMULATION", "TRUE")) scale = 2./3;

  nrs->mueDivKernel(
       mesh->Nelements*mesh->Np,
       scale,
       nrs->o_mue,
       nrs->o_div,
       platform->o_mempool.slice3);

  nrs->gradientVolumeKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    mesh->o_D,
    nrs->fieldOffset,
    platform->o_mempool.slice3,
    platform->o_mempool.slice0);

  nrs->wgradientVolumeKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    mesh->o_D,
    nrs->fieldOffset,
    nrs->o_P,
    platform->o_mempool.slice3);

  platform->linAlg->axpby(
    nrs->NVfields*nrs->fieldOffset,
    1.0,
    platform->o_mempool.slice3,
    -1.0,
    platform->o_mempool.slice0);

  nrs->velocityNeumannBCKernel(
       mesh->Nelements,
       nrs->fieldOffset,
       mesh->o_sgeo,
       mesh->o_vmapM,
       mesh->o_EToB,
       nrs->o_EToB,
       time,
       mesh->o_x,
       mesh->o_y,
       mesh->o_z,
       nrs->o_usrwrk,
       nrs->o_U,
       platform->o_mempool.slice0);

  nrs->velocityRhsKernel(
    mesh->Nelements,
    nrs->fieldOffset,
    nrs->o_BF,
    platform->o_mempool.slice0,
    nrs->o_rho,
    platform->o_mempool.slice3);

  if(platform->options.compareArgs("VELOCITY INITIAL GUESS DEFAULT", "EXTRAPOLATION") && stage == 1) {
    platform->o_mempool.slice0.copyFrom(nrs->o_Ue, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));
    if (nrs->uvwSolver) {
      if (nrs->uvwSolver->Nmasked) nrs->maskCopyKernel(nrs->uvwSolver->Nmasked, 0*nrs->fieldOffset, nrs->uvwSolver->o_maskIds,
                                                       nrs->o_U, platform->o_mempool.slice0);
    } else {
      if (nrs->uSolver->Nmasked) nrs->maskCopyKernel(nrs->uSolver->Nmasked, 0*nrs->fieldOffset, nrs->uSolver->o_maskIds,
                                                     nrs->o_U, platform->o_mempool.slice0);
      if (nrs->vSolver->Nmasked) nrs->maskCopyKernel(nrs->vSolver->Nmasked, 1*nrs->fieldOffset, nrs->vSolver->o_maskIds,
                                                     nrs->o_U, platform->o_mempool.slice0);
      if (nrs->wSolver->Nmasked) nrs->maskCopyKernel(nrs->wSolver->Nmasked, 2*nrs->fieldOffset, nrs->wSolver->o_maskIds,
                                                     nrs->o_U, platform->o_mempool.slice0);
    }
  } else {
    platform->o_mempool.slice0.copyFrom(nrs->o_U, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));
  }

  if(nrs->uvwSolver) {
    ellipticSolve(nrs->uvwSolver, platform->o_mempool.slice3, platform->o_mempool.slice0);
  } else {
    ellipticSolve(nrs->uSolver, platform->o_mempool.slice3, platform->o_mempool.slice0);
    ellipticSolve(nrs->vSolver, platform->o_mempool.slice4, platform->o_mempool.slice1);
    ellipticSolve(nrs->wSolver, platform->o_mempool.slice5, platform->o_mempool.slice2);
  }

  return platform->o_mempool.slice0;
}

} // namespace

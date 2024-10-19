#include "nrs.hpp"
#include "bcMap.hpp"

// lower than any other possible Dirichlet value
static constexpr dfloat TINY = -1e30;

void createZeroNormalMask(nrs_t *nrs,
                          mesh_t *mesh,
                          const occa::memory &o_EToB,
                          const occa::memory &o_EToBV,
                          occa::memory &o_mask)
{
  nrs->initializeZeroNormalMaskKernel(mesh->Nlocal, nrs->fieldOffset, o_EToBV, o_mask);

  // normal xyz + count
  occa::memory o_avgNormal =
      platform->deviceMemoryPool.reserve<dfloat>((nrs->NVfields + 1) * nrs->fieldOffset);

  int bcType = ellipticBcType::ZERO_NORMAL;
  nrs->averageNormalBcTypeKernel(mesh->Nelements,
                                 nrs->fieldOffset,
                                 bcType,
                                 mesh->o_sgeo,
                                 mesh->o_vmapM,
                                 o_EToB,
                                 o_avgNormal);

  oogs::startFinish(o_avgNormal, nrs->NVfields + 1, nrs->fieldOffset, ogsDfloat, ogsAdd, mesh->oogs);

  nrs->fixZeroNormalMaskKernel(mesh->Nelements,
                               nrs->fieldOffset,
                               mesh->o_sgeo,
                               mesh->o_vmapM,
                               o_EToB,
                               o_avgNormal,
                               o_mask);

  oogs::startFinish(o_mask, nrs->NVfields, nrs->fieldOffset, ogsDfloat, ogsMin, mesh->oogs);
}

void applyZeroNormalMask(nrs_t *nrs,
                         mesh_t *mesh,
                         dlong Nelements,
                         const occa::memory &o_elementList,
                         const occa::memory &o_EToB,
                         const occa::memory &o_mask,
                         occa::memory &o_x)
{
  if (Nelements == 0) {
    return;
  }

  nrs->applyZeroNormalMaskKernel(Nelements,
                                 nrs->fieldOffset,
                                 o_elementList,
                                 mesh->o_sgeo,
                                 o_mask,
                                 mesh->o_vmapM,
                                 o_EToB,
                                 o_x);
}

void applyZeroNormalMask(nrs_t *nrs,
                         mesh_t *mesh,
                         const occa::memory &o_EToB,
                         const occa::memory &o_mask,
                         occa::memory &o_x)
{
  nrs->applyZeroNormalMaskKernel(mesh->Nelements,
                                 nrs->fieldOffset,
                                 mesh->o_elementList,
                                 mesh->o_sgeo,
                                 o_mask,
                                 mesh->o_vmapM,
                                 o_EToB,
                                 o_x);
}

void applyDirichletVelocity(nrs_t *nrs, double time, occa::memory &o_U, occa::memory &o_Ue, occa::memory &o_P)
{
  if (bcMap::unalignedMixedBoundary("velocity")) {
    applyZeroNormalMask(nrs, nrs->mesh, nrs->uvwSolver->o_EToB(), nrs->o_zeroNormalMaskVelocity, o_U);
    applyZeroNormalMask(nrs, nrs->mesh, nrs->uvwSolver->o_EToB(), nrs->o_zeroNormalMaskVelocity, o_Ue);
  }

  const auto neknekFieldOffset = nrs->neknek ? nrs->neknek->fieldOffset() : 0;

  auto mesh = nrs->mesh;

  occa::memory o_tmp = platform->deviceMemoryPool.reserve<dfloat>((nrs->NVfields + 1) * nrs->fieldOffset);
  platform->linAlg->fill((1 + nrs->NVfields) * nrs->fieldOffset, TINY, o_tmp);

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
                                   nrs->o_rho,
                                   nrs->o_mue,
                                   nrs->o_usrwrk,
                                   o_Ue,
                                   o_tmp);

    nrs->velocityDirichletBCKernel(mesh->Nelements,
                                   nrs->fieldOffset,
                                   time,
                                   mesh->o_sgeo,
                                   nrs->o_zeroNormalMaskVelocity,
                                   mesh->o_x,
                                   mesh->o_y,
                                   mesh->o_z,
                                   mesh->o_vmapM,
                                   mesh->o_EToB,
                                   nrs->o_EToB,
                                   nrs->o_rho,
                                   nrs->o_mue,
                                   neknekFieldOffset,
                                   nrs->neknek ? nrs->neknek->o_pointMap() : o_NULL,
                                   nrs->neknek ? nrs->neknek->o_U() : o_NULL,
                                   nrs->o_usrwrk,
                                   o_U,
                                   o_tmp.slice(nrs->fieldOffset));

    oogs::startFinish(o_tmp,
                      1 + nrs->NVfields,
                      nrs->fieldOffset,
                      ogsDfloat,
                      (sweep == 0) ? ogsMax : ogsMin,
                      nrs->gsh);
  }

  if (nrs->pSolver->Nmasked()) {
    auto o_dirichlet = o_tmp.slice(0, nrs->fieldOffset);
    nrs->maskCopyKernel(nrs->pSolver->Nmasked(), 0, 0, nrs->pSolver->o_maskIds(), o_dirichlet, o_P);
  }

  auto o_UDirichlet = o_tmp.slice(nrs->fieldOffset);
  if (nrs->uvwSolver) {
    if (nrs->uvwSolver->Nmasked()) {
      nrs->maskCopy2Kernel(nrs->uvwSolver->Nmasked(),
                           0 * nrs->fieldOffset,
                           0 * nrs->fieldOffset,
                           nrs->uvwSolver->o_maskIds(),
                           o_UDirichlet,
                           o_U,
                           o_Ue);
    }
  } else {
    int cnt = 0;
    for (auto &solver : {nrs->uSolver, nrs->vSolver, nrs->wSolver}) {
      if (solver->Nmasked()) {
        nrs->maskCopy2Kernel(solver->Nmasked(),
                             cnt * nrs->fieldOffset,
                             cnt * nrs->fieldOffset,
                             solver->o_maskIds(),
                             o_UDirichlet,
                             o_U,
                             o_Ue);
      }
      cnt++;
    }
  }
}

void applyDirichletScalars(nrs_t *nrs, double time, occa::memory &o_S, occa::memory &o_Se)
{
  cds_t *cds = nrs->cds;

  const auto neknekFieldOffset = nrs->neknek ? nrs->neknek->fieldOffset() : 0;
  for (int is = 0; is < cds->NSfields; is++) {
    if (!cds->compute[is]) {
      continue;
    }
    if (cds->cvodeSolve[is]) {
      continue;
    }
    mesh_t *mesh = cds->mesh[0];
    oogs_t *gsh = cds->gshT;
    if (is) {
      mesh = cds->meshV;
      gsh = cds->gsh;
    }

    auto o_diff_i = cds->o_diff + cds->fieldOffsetScan[is];
    auto o_rho_i = cds->o_rho + cds->fieldOffsetScan[is];

    occa::memory o_SiDirichlet = platform->deviceMemoryPool.reserve<dfloat>(cds->fieldOffset[is]);
    platform->linAlg->fill(cds->fieldOffset[is], TINY, o_SiDirichlet);

    for (int sweep = 0; sweep < 2; sweep++) {
      cds->dirichletBCKernel(mesh->Nelements,
                             cds->fieldOffset[is],
                             is,
                             time,
                             mesh->o_sgeo,
                             mesh->o_x,
                             mesh->o_y,
                             mesh->o_z,
                             mesh->o_vmapM,
                             mesh->o_EToB,
                             cds->o_EToB + is * cds->EToBOffset,
                             cds->o_Ue,
                             o_diff_i,
                             o_rho_i,
                             neknekFieldOffset,
                             nrs->neknek ? nrs->neknek->o_pointMap() : o_NULL,
                             nrs->neknek ? static_cast<int>(nrs->neknek->o_U().isInitialized()) : 0,
                             nrs->neknek ? nrs->neknek->o_U() : o_NULL,
                             nrs->neknek ? nrs->neknek->o_S() : o_NULL,
                             nrs->neknek ? nrs->neknek->o_scalarIndices() : o_NULL,
                             *(cds->o_usrwrk),
                             o_SiDirichlet);

      oogs::startFinish(o_SiDirichlet,
                        1,
                        cds->fieldOffset[is],
                        ogsDfloat,
                        (sweep == 0) ? ogsMax : ogsMin,
                        gsh);
    }
    occa::memory o_Si = o_S.slice(cds->fieldOffsetScan[is], cds->fieldOffset[is]);

    if (o_Se.isInitialized()) {
      occa::memory o_Si_e = o_Se.slice(cds->fieldOffsetScan[is], cds->fieldOffset[is]);

      if (cds->solver[is]->Nmasked()) {
        cds->maskCopy2Kernel(cds->solver[is]->Nmasked(),
                             0,
                             0,
                             cds->solver[is]->o_maskIds(),
                             o_SiDirichlet,
                             o_Si,
                             o_Si_e);
      }
    } else {
      if (cds->solver[is]->Nmasked()) {
        cds->maskCopyKernel(cds->solver[is]->Nmasked(),
                            0,
                            0,
                            cds->solver[is]->o_maskIds(),
                            o_SiDirichlet,
                            o_Si);
      }
    }
  }
}

void applyDirichletMesh(nrs_t *nrs, double time, occa::memory &o_UM, occa::memory &o_UMe, occa::memory &o_U)
{
  auto mesh = (nrs->cht) ? nrs->cds->mesh[0] : nrs->mesh;
  if (bcMap::unalignedMixedBoundary("mesh")) {
    applyZeroNormalMask(nrs, mesh, nrs->meshSolver->o_EToB(), nrs->o_zeroNormalMaskMeshVelocity, o_UM);
    applyZeroNormalMask(nrs, mesh, nrs->meshSolver->o_EToB(), nrs->o_zeroNormalMaskMeshVelocity, o_UMe);
  }

  occa::memory o_UDirichlet = platform->deviceMemoryPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);
  platform->linAlg->fill(nrs->NVfields * nrs->fieldOffset, TINY, o_UDirichlet);

  for (int sweep = 0; sweep < 2; sweep++) {
    mesh->velocityDirichletKernel(mesh->Nelements,
                                  nrs->fieldOffset,
                                  time,
                                  (int)bcMap::useDerivedMeshBoundaryConditions(),
                                  mesh->o_sgeo,
                                  nrs->o_zeroNormalMaskMeshVelocity,
                                  mesh->o_x,
                                  mesh->o_y,
                                  mesh->o_z,
                                  mesh->o_vmapM,
                                  mesh->o_EToB,
                                  nrs->o_EToBMeshVelocity,
                                  nrs->o_meshRho,
                                  nrs->o_meshMue,
                                  nrs->o_usrwrk,
                                  o_U,
                                  o_UDirichlet);

    oogs::startFinish(o_UDirichlet,
                      nrs->NVfields,
                      nrs->fieldOffset,
                      ogsDfloat,
                      (sweep == 0) ? ogsMax : ogsMin,
                      nrs->gshMesh);
  }

  if (nrs->meshSolver->Nmasked()) {
    nrs->maskCopy2Kernel(nrs->meshSolver->Nmasked(),
                         0 * nrs->fieldOffset,
                         0 * nrs->fieldOffset,
                         nrs->meshSolver->o_maskIds(),
                         o_UDirichlet,
                         o_UM,
                         o_UMe);
  }
}

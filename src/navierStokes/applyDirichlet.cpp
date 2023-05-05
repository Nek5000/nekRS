#include "nrs.hpp"
#include "bcMap.hpp"

// lower than any other possible Dirichlet value
static constexpr dfloat TINY = -1e30;

void createZeroNormalMask(nrs_t *nrs, mesh_t *mesh, occa::memory &o_EToB, occa::memory& o_EToBV, occa::memory &o_mask)
{
  nrs->initializeZeroNormalMaskKernel(mesh->Nlocal, nrs->fieldOffset, o_EToBV, o_mask);

  // normal + count (4 fields)
  auto o_avgNormal = platform->o_mempool.slice0;
  int bcType = ZERO_NORMAL;
  nrs->averageNormalBcTypeKernel(mesh->Nelements,
                                 nrs->fieldOffset,
                                 bcType,
                                 mesh->o_sgeo,
                                 mesh->o_vmapM,
                                 o_EToB,
                                 o_avgNormal);

  oogs::startFinish(o_avgNormal, nrs->NVfields+1, nrs->fieldOffset, ogsDfloat, ogsAdd, mesh->oogs);

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
                         occa::memory &o_elementList,
                         occa::memory &o_EToB,
                         occa::memory &o_mask,
                         occa::memory &o_x)
{
  if (Nelements == 0)
    return;

  nrs->applyZeroNormalMaskKernel(Nelements,
                                 nrs->fieldOffset,
                                 o_elementList,
                                 mesh->o_sgeo,
                                 o_mask,
                                 mesh->o_vmapM,
                                 o_EToB,
                                 o_x);
}

void applyZeroNormalMask(nrs_t *nrs, mesh_t *mesh, occa::memory &o_EToB, occa::memory &o_mask, occa::memory &o_x)
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

void applyDirichletVelocity(nrs_t *nrs, double time, occa::memory& o_U,occa::memory& o_Ue,occa::memory& o_P)
{
  if (bcMap::unalignedMixedBoundary("velocity")) {
    applyZeroNormalMask(nrs, nrs->meshV, nrs->uvwSolver->o_EToB, nrs->o_zeroNormalMaskVelocity, o_U);
    applyZeroNormalMask(nrs, nrs->meshV, nrs->uvwSolver->o_EToB, nrs->o_zeroNormalMaskVelocity, o_Ue);
  }

  mesh_t *mesh = nrs->meshV;

  platform->linAlg->fill((1 + nrs->NVfields) * nrs->fieldOffset, TINY, platform->o_mempool.slice0);

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
                                   platform->o_mempool.slice0);

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
                                   nrs->neknek ? nrs->neknek->o_pointMap : o_NULL,
                                   nrs->neknek ? nrs->neknek->o_U : o_NULL,
                                   nrs->o_usrwrk,
                                   o_U,
                                   platform->o_mempool.slice1);

    if (sweep == 0)
      oogs::startFinish(platform->o_mempool.slice0,
                        1 + nrs->NVfields,
                        nrs->fieldOffset,
                        ogsDfloat,
                        ogsMax,
                        nrs->gsh);
    if (sweep == 1)
      oogs::startFinish(platform->o_mempool.slice0,
                        1 + nrs->NVfields,
                        nrs->fieldOffset,
                        ogsDfloat,
                        ogsMin,
                        nrs->gsh);
  }

  if (nrs->pSolver->Nmasked)
    nrs->maskCopyKernel(nrs->pSolver->Nmasked,
                        0,
                        0,
                        nrs->pSolver->o_maskIds,
                        platform->o_mempool.slice0,
                        o_P);

  if (nrs->uvwSolver) {
    if (nrs->uvwSolver->Nmasked) {
      nrs->maskCopy2Kernel(nrs->uvwSolver->Nmasked,
                          0 * nrs->fieldOffset,
                          0 * nrs->fieldOffset,
                          nrs->uvwSolver->o_maskIds,
                          platform->o_mempool.slice1,
                          o_U, o_Ue);
    }
  } else {
    if (nrs->uSolver->Nmasked) {
      nrs->maskCopy2Kernel(nrs->uSolver->Nmasked,
                          0 * nrs->fieldOffset,
                          0 * nrs->fieldOffset,
                          nrs->uSolver->o_maskIds,
                          platform->o_mempool.slice1,
                          o_U, o_Ue);
    }
    if (nrs->vSolver->Nmasked) {
      nrs->maskCopy2Kernel(nrs->vSolver->Nmasked,
                          1 * nrs->fieldOffset,
                          1 * nrs->fieldOffset,
                          nrs->vSolver->o_maskIds,
                          platform->o_mempool.slice1,
                          o_U, o_Ue);
    }
    if (nrs->wSolver->Nmasked) {
      nrs->maskCopy2Kernel(nrs->wSolver->Nmasked,
                          2 * nrs->fieldOffset,
                          2 * nrs->fieldOffset,
                          nrs->wSolver->o_maskIds,
                          platform->o_mempool.slice1,
                          o_U, o_Ue);
    }
  }
}

void applyDirichletScalars(nrs_t *nrs, double time, occa::memory& o_S, occa::memory& o_Se)
{
  cds_t *cds = nrs->cds;
  for (int is = 0; is < cds->NSfields; is++) {
    if (!cds->compute[is])
      continue;
    if (cds->cvodeSolve[is])
      continue;
    mesh_t *mesh = cds->mesh[0];
    oogs_t *gsh = cds->gshT;
    if (is) {
      mesh = cds->meshV;
      gsh = cds->gsh;
    }

    auto o_diff_i = cds->o_diff + cds->fieldOffsetScan[is] * sizeof(dfloat);
    auto o_rho_i = cds->o_rho + cds->fieldOffsetScan[is] * sizeof(dfloat);

    platform->linAlg->fill(cds->fieldOffset[is], TINY, platform->o_mempool.slice0);

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
                             cds->o_EToB + is * cds->EToBOffset * sizeof(int),
                             cds->o_Ue,
                             o_diff_i,
                             o_rho_i,
                             cds->neknek ? cds->neknek->o_pointMap : o_NULL,
                             cds->neknek ? cds->neknek->o_U : o_NULL,
                             cds->neknek ? cds->neknek->o_S : o_NULL,
                             *(cds->o_usrwrk),
                             platform->o_mempool.slice0);

      if (sweep == 0)
        oogs::startFinish(platform->o_mempool.slice0, 1, cds->fieldOffset[is], ogsDfloat, ogsMax, gsh);
      if (sweep == 1)
        oogs::startFinish(platform->o_mempool.slice0, 1, cds->fieldOffset[is], ogsDfloat, ogsMin, gsh);
    }
    occa::memory o_Si =
        o_S.slice(cds->fieldOffsetScan[is] * sizeof(dfloat), cds->fieldOffset[is] * sizeof(dfloat));
    
    if(o_Se.isInitialized()){
      occa::memory o_Si_e =
          o_Se.slice(cds->fieldOffsetScan[is] * sizeof(dfloat), cds->fieldOffset[is] * sizeof(dfloat));

      if (cds->solver[is]->Nmasked)
        cds->maskCopy2Kernel(cds->solver[is]->Nmasked,
                            0,
                            0,
                            cds->solver[is]->o_maskIds,
                            platform->o_mempool.slice0,
                            o_Si, o_Si_e);
    } else {
      if (cds->solver[is]->Nmasked)
        cds->maskCopyKernel(cds->solver[is]->Nmasked,
                            0,
                            0,
                            cds->solver[is]->o_maskIds,
                            platform->o_mempool.slice0,
                            o_Si);
    }
  }
}

void applyDirichletMesh(nrs_t *nrs, double time, occa::memory& o_UM, occa::memory& o_UMe, occa::memory& o_U)
{
  mesh_t *mesh = nrs->_mesh;
  if (bcMap::unalignedMixedBoundary("mesh")) {
    applyZeroNormalMask(nrs, mesh, nrs->meshSolver->o_EToB, nrs->o_zeroNormalMaskMeshVelocity, o_UM);
    applyZeroNormalMask(nrs, mesh, nrs->meshSolver->o_EToB, nrs->o_zeroNormalMaskMeshVelocity, o_UMe);
  }

  platform->linAlg->fill(nrs->NVfields * nrs->fieldOffset, TINY, platform->o_mempool.slice0);

  for (int sweep = 0; sweep < 2; sweep++) {
    mesh->velocityDirichletKernel(mesh->Nelements,
                                  nrs->fieldOffset,
                                  time,
                                  (int) bcMap::useDerivedMeshBoundaryConditions(),
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
                                  platform->o_mempool.slice0);

    if (sweep == 0)
      oogs::startFinish(platform->o_mempool.slice0,
                        nrs->NVfields,
                        nrs->fieldOffset,
                        ogsDfloat,
                        ogsMax,
                        nrs->gshMesh);
    if (sweep == 1)
      oogs::startFinish(platform->o_mempool.slice0,
                        nrs->NVfields,
                        nrs->fieldOffset,
                        ogsDfloat,
                        ogsMin,
                        nrs->gshMesh);
  }

  if (nrs->meshSolver->Nmasked)
    nrs->maskCopy2Kernel(nrs->meshSolver->Nmasked,
                        0 * nrs->fieldOffset,
                        0 * nrs->fieldOffset,
                        nrs->meshSolver->o_maskIds,
                        platform->o_mempool.slice0,
                        o_UM, o_UMe);
}

void applyDirichlet(nrs_t *nrs, double time)
{
  if (nrs->flow)
    applyDirichletVelocity(nrs, time, nrs->o_U, nrs->o_Ue, nrs->o_P);  

  if (nrs->Nscalar) 
    applyDirichletScalars(nrs, time, nrs->cds->o_S, nrs->cds->o_Se); 

  if (!platform->options.compareArgs("MESH SOLVER", "NONE"))
    applyDirichletMesh(nrs, time, nrs->_mesh->o_U, nrs->_mesh->o_Ue, nrs->o_U);
}

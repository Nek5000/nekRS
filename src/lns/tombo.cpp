#include "nrs.hpp"
#include "udf.hpp"

namespace tombo
{
occa::memory pressureSolve(nrs_t* nrs, dfloat time)
{
  mesh_t* mesh = nrs->mesh;

  //enforce Dirichlet BCs
  nrs->fillKernel((1+nrs->NVfields)*nrs->fieldOffset, std::numeric_limits<dfloat>::min(), nrs->o_wrk6);
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
                                   nrs->o_wrk6);

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
                                   nrs->o_usrwrk,
                                   nrs->o_U,
                                   nrs->o_wrk7);

    //take care of Neumann-Dirichlet shared edges across elements
    if (sweep == 0) oogs::startFinish(nrs->o_wrk6, 1+nrs->NVfields, nrs->fieldOffset, ogsDfloat, ogsMax, nrs->gsh);
    if (sweep == 1) oogs::startFinish(nrs->o_wrk6, 1+nrs->NVfields, nrs->fieldOffset, ogsDfloat, ogsMin, nrs->gsh);
  }

  if (nrs->pSolver->Nmasked) nrs->maskCopyKernel(nrs->pSolver->Nmasked, 0, nrs->pSolver->o_maskIds,
                                                 nrs->o_wrk6, nrs->o_P); 

  if (nrs->uvwSolver) {
    if (nrs->uvwSolver->Nmasked) nrs->maskCopyKernel(nrs->uvwSolver->Nmasked, 0*nrs->fieldOffset, nrs->uvwSolver->o_maskIds,
                                                     nrs->o_wrk7, nrs->o_U);
  } else {
    if (nrs->uSolver->Nmasked) nrs->maskCopyKernel(nrs->uSolver->Nmasked, 0*nrs->fieldOffset, nrs->uSolver->o_maskIds, 
                                                   nrs->o_wrk7, nrs->o_U);
    if (nrs->vSolver->Nmasked) nrs->maskCopyKernel(nrs->vSolver->Nmasked, 1*nrs->fieldOffset, nrs->vSolver->o_maskIds, 
                                                   nrs->o_wrk7, nrs->o_U);
    if (nrs->wSolver->Nmasked) nrs->maskCopyKernel(nrs->wSolver->Nmasked, 2*nrs->fieldOffset, nrs->wSolver->o_maskIds, 
                                                   nrs->o_wrk7, nrs->o_U);
  }

  nrs->curlKernel(mesh->Nelements,
                  mesh->o_vgeo,
                  mesh->o_Dmatrices,
                  nrs->fieldOffset,
                  nrs->o_Ue,
                  nrs->o_wrk0);

  oogs::startFinish(nrs->o_wrk0, nrs->NVfields, nrs->fieldOffset,ogsDfloat, ogsAdd, nrs->gsh);

  nrs->invMassMatrixKernel(
    mesh->Nelements,
    nrs->fieldOffset,
    nrs->NVfields,
    mesh->o_vgeo,
    nrs->mesh->o_invLMM,
    nrs->o_wrk0);

  nrs->curlKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    mesh->o_Dmatrices,
    nrs->fieldOffset,
    nrs->o_wrk0,
    nrs->o_wrk3);

  nrs->gradientVolumeKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    mesh->o_Dmatrices,
    nrs->fieldOffset,
    nrs->o_div,
    nrs->o_wrk0);

  if(nrs->options.compareArgs("STRESSFORMULATION", "TRUE"))
    nrs->pressureStressKernel(
         mesh->Nelements,
         mesh->o_vgeo,
         mesh->o_Dmatrices,
         nrs->fieldOffset,
         nrs->o_mue,
         nrs->o_Ue,
         nrs->o_div,
         nrs->o_wrk3);

  occa::memory o_irho = nrs->o_ellipticCoeff;
  nrs->pressureRhsKernel(
    mesh->Nelements * mesh->Np,
    nrs->fieldOffset,
    nrs->o_mue,
    o_irho,
    nrs->o_BF,
    nrs->o_wrk3,
    nrs->o_wrk0,
    nrs->o_wrk6);

  oogs::startFinish(nrs->o_wrk6, nrs->NVfields, nrs->fieldOffset,ogsDfloat, ogsAdd, nrs->gsh);

  nrs->invMassMatrixKernel(
    mesh->Nelements,
    nrs->fieldOffset,
    nrs->NVfields,
    mesh->o_vgeo,
    nrs->mesh->o_invLMM,
    nrs->o_wrk6);

  nrs->divergenceVolumeKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    mesh->o_Dmatrices,
    nrs->fieldOffset,
    nrs->o_wrk6,
    nrs->o_wrk3);

  nrs->pressureAddQtlKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    nrs->g0 * nrs->idt,
    nrs->o_div,
    nrs->o_wrk3);

  nrs->divergenceSurfaceKernel(
    mesh->Nelements,
    mesh->o_sgeo,
    mesh->o_vmapM,
    nrs->o_EToB,
    nrs->g0 * nrs->idt,
    nrs->fieldOffset,
    nrs->o_wrk6,
    nrs->o_U,
    nrs->o_wrk3);

  nrs->o_wrk1.copyFrom(nrs->o_P, nrs->Ntotal * sizeof(dfloat));
  nrs->NiterP = ellipticSolve(nrs->pSolver, nrs->o_wrk3, nrs->o_wrk1);

  return nrs->o_wrk1;
}

occa::memory velocitySolve(nrs_t* nrs, dfloat time)
{
  mesh_t* mesh = nrs->mesh;

  dfloat scale = -1./3;
  if(nrs->options.compareArgs("STRESSFORMULATION", "TRUE")) scale = 2./3;

#if 0
  nrs->PQKernel(
       mesh->Nelements*mesh->Np,
       -scale,
       nrs->o_mue,
       nrs->o_div,
       nrs->o_P,
       nrs->o_wrk3); 

  nrs->gradientVolumeKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    mesh->o_Dmatrices,
    nrs->fieldOffset,
    nrs->o_wrk3,
    nrs->o_wrk0);
#else
  nrs->mueDivKernel(
       mesh->Nelements*mesh->Np,
       scale,
       nrs->o_mue,
       nrs->o_div,
       nrs->o_wrk3); 

  nrs->gradientVolumeKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    mesh->o_Dmatrices,
    nrs->fieldOffset,
    nrs->o_wrk3,
    nrs->o_wrk0);

  nrs->wgradientVolumeKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    mesh->o_Dmatrices,
    nrs->fieldOffset,
    nrs->o_P,
    nrs->o_wrk3); 

  nrs->scaledAddKernel(
    nrs->NVfields*nrs->fieldOffset,
    1.0,
    0*nrs->fieldOffset,
    nrs->o_wrk3,
    -1.0,
    0*nrs->fieldOffset,
    nrs->o_wrk0);
#endif

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
       nrs->o_wrk0); 

  nrs->velocityRhsKernel(
    mesh->Nelements,
    nrs->fieldOffset,
    nrs->o_BF,
    nrs->o_wrk0,
    nrs->o_rho,
    nrs->o_wrk3);

  if(nrs->options.compareArgs("VELOCITY INITIAL GUESS DEFAULT", "EXTRAPOLATION")) { 
    nrs->o_wrk0.copyFrom(nrs->o_Ue, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));
    if (nrs->uvwSolver) {
      if (nrs->uvwSolver->Nmasked) nrs->maskCopyKernel(nrs->uvwSolver->Nmasked, 0*nrs->fieldOffset, nrs->uvwSolver->o_maskIds,
                                                       nrs->o_U, nrs->o_wrk0);
    } else {
      if (nrs->uSolver->Nmasked) nrs->maskCopyKernel(nrs->uSolver->Nmasked, 0*nrs->fieldOffset, nrs->uSolver->o_maskIds,
                                                     nrs->o_U, nrs->o_wrk0);
      if (nrs->vSolver->Nmasked) nrs->maskCopyKernel(nrs->vSolver->Nmasked, 1*nrs->fieldOffset, nrs->vSolver->o_maskIds,
                                                     nrs->o_U, nrs->o_wrk0);
      if (nrs->wSolver->Nmasked) nrs->maskCopyKernel(nrs->wSolver->Nmasked, 2*nrs->fieldOffset, nrs->wSolver->o_maskIds,
                                                     nrs->o_U, nrs->o_wrk0);
    }
  } else {
    nrs->o_wrk0.copyFrom(nrs->o_U, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));
  }

  if(nrs->uvwSolver) {
    nrs->NiterU = ellipticSolve(nrs->uvwSolver, nrs->o_wrk3, nrs->o_wrk0);
  } else {
    nrs->NiterU = ellipticSolve(nrs->uSolver, nrs->o_wrk3, nrs->o_wrk0);
    nrs->NiterV = ellipticSolve(nrs->vSolver, nrs->o_wrk4, nrs->o_wrk1);
    nrs->NiterW = ellipticSolve(nrs->wSolver, nrs->o_wrk5, nrs->o_wrk2);
  }

  return nrs->o_wrk0;
}

} // namespace

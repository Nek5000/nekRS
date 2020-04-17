#include "nrs.hpp"
#include "udf.hpp"
#include "nekInterfaceAdapter.hpp"

#define NEKDSSUM

int firstStep = 1;

namespace tombo {

void nek_ogsGatherScatterMany(ins_t *ins, occa::memory o_u, int nfld)
{
   occa::memory o_vx = o_u + 0*ins->fieldOffset*sizeof(dfloat);
   occa::memory o_vy = o_u + 1*ins->fieldOffset*sizeof(dfloat);
   occa::memory o_vz = o_u + 2*ins->fieldOffset*sizeof(dfloat);

   dfloat *vx = ins->U + 0*ins->fieldOffset;
   dfloat *vy = ins->U + 1*ins->fieldOffset;
   dfloat *vz = ins->U + 2*ins->fieldOffset;

   *(nekData.istep) = 1;
   o_vx.copyTo(vx, ins->fieldOffset*sizeof(dfloat));
   memcpy(nekData.vx, vx, sizeof(dfloat)*ins->Nlocal);
   if(nfld==3) {
     o_vy.copyTo(vy, ins->fieldOffset*sizeof(dfloat));
     memcpy(nekData.vy, vy, sizeof(dfloat)*ins->Nlocal);

     o_vz.copyTo(vz, ins->fieldOffset*sizeof(dfloat));
     memcpy(nekData.vz, vz, sizeof(dfloat)*ins->Nlocal);
   }

   nek_userchk(); // call dssum for vx, vy, vz

   memcpy(vx, nekData.vx, sizeof(dfloat)*ins->Nlocal);
   o_vx.copyFrom(vx, ins->fieldOffset*sizeof(dfloat));
   if(nfld==3) {
     memcpy(vy, nekData.vy, sizeof(dfloat)*ins->Nlocal);
     o_vy.copyFrom(vy, ins->fieldOffset*sizeof(dfloat));
     memcpy(vz, nekData.vz, sizeof(dfloat)*ins->Nlocal);
     o_vz.copyFrom(vz, ins->fieldOffset*sizeof(dfloat));
   }
}

occa::memory pressureSolve(ins_t *ins, dfloat time)
{
  mesh_t *mesh = ins->mesh;

  ins->curlKernel(mesh->Nelements,
                  mesh->o_vgeo,
                  mesh->o_Dmatrices,
                  mesh->o_MM, 
                  ins->fieldOffset,
                  ins->o_Ue,
                  ins->o_wrk0);

#ifdef NEKDSSUM
  ins->o_wrk0.copyTo(ins->U, ins->NVfields*ins->fieldOffset*sizeof(dfloat));
  nek_dssum(ins->U + 0*ins->fieldOffset);
  nek_dssum(ins->U + 1*ins->fieldOffset);
  nek_dssum(ins->U + 2*ins->fieldOffset);
  ins->o_wrk0.copyFrom(ins->U, ins->NVfields*ins->fieldOffset*sizeof(dfloat));
#else
  ogsGatherScatterMany(ins->o_wrk0, ins->NVfields, ins->fieldOffset,
                       ogsDfloat, ogsAdd, mesh->ogs);
#endif

  ins->invMassMatrixKernel(
       mesh->Nelements,
       ins->fieldOffset,
       ins->NVfields,
       mesh->o_vgeo,
       ins->o_InvM,
       ins->o_wrk0);

  ins->curlKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_Dmatrices,
       mesh->o_MM, 
       ins->fieldOffset,
       ins->o_wrk0,
       ins->o_wrk3);

  ins->gradientVolumeKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_Dmatrices,
       ins->fieldOffset,
       ins->o_div,
       ins->o_wrk0);
  occa::memory o_irho = ins->o_ellipticCoeff;
  ins->ncKernel(
       mesh->Np*mesh->Nelements,
       ins->fieldOffset,
       ins->o_mue,
       o_irho,
       ins->o_wrk0,
       ins->o_wrk3);

  ins->pressureRhsKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_MM,
       ins->idt,
       ins->g0,
       ins->o_extbdfA,
       ins->o_extbdfB,
       ins->fieldOffset,
       ins->o_U,
       ins->o_BF,
       ins->o_wrk3,
       ins->o_FU,
       ins->o_wrk0);

#ifdef NEKDSSUM 
  ins->o_wrk0.copyTo(ins->U, ins->NVfields*ins->fieldOffset*sizeof(dfloat));
  nek_dssum(ins->U + 0*ins->fieldOffset);
  nek_dssum(ins->U + 1*ins->fieldOffset);
  nek_dssum(ins->U + 2*ins->fieldOffset);
  ins->o_wrk0.copyFrom(ins->U, ins->NVfields*ins->fieldOffset*sizeof(dfloat));
#else
  ogsGatherScatterMany(ins->o_wrk0, ins->NVfields, ins->fieldOffset,
                       ogsDfloat, ogsAdd, mesh->ogs);
#endif

  ins->invMassMatrixKernel(
       mesh->Nelements,
       ins->fieldOffset,
       ins->NVfields,
       mesh->o_vgeo,
       ins->o_InvM,
       ins->o_wrk0);

  ins->divergenceVolumeKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_Dmatrices,
       ins->fieldOffset,
       ins->o_wrk0,
       ins->o_wrk3);

  const dfloat lambda = ins->g0*ins->idt;
  ins->divergenceSurfaceKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_sgeo,
       mesh->o_LIFTT,
       mesh->o_vmapM,
       mesh->o_vmapP,
       mesh->o_EToB,
       ins->o_EToB,
       time,
       -lambda, 
       mesh->o_x,
       mesh->o_y,
       mesh->o_z,
       ins->fieldOffset,
       ins->o_usrwrk,
       ins->o_wrk0,
       ins->o_wrk3);

  ins->AxKernel(
       mesh->Nelements,
       ins->fieldOffset,
       mesh->o_ggeo,
       mesh->o_Dmatrices,
       mesh->o_Smatrices,
       mesh->o_MM,
       ins->o_P,
       ins->o_ellipticCoeff,
       ins->o_wrk3);

  ins->pressureAddQtlKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       lambda,
       ins->o_div,
       ins->o_wrk3);

  elliptic_t *solver = ins->pSolver;

#ifdef NEKDSSUM 
  ins->o_wrk3.copyTo(ins->U, ins->fieldOffset*sizeof(dfloat));
  nek_dssum(ins->U + 0*ins->fieldOffset);
  ins->o_wrk3.copyFrom(ins->U, ins->fieldOffset*sizeof(dfloat));
#else
  ogsGatherScatter(ins->o_wrk3, ogsDfloat, ogsAdd, mesh->ogs);
#endif

  if(firstStep){
    ins->o_wrk3.copyTo(ins->U, ins->fieldOffset*sizeof(dfloat)); //dumping respr (no masked)
    nek_copyFrom(time, 0);
    nek_outfld();
    firstStep = 0;
  }

  if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, ins->o_wrk3);

  ins->setScalarKernel(ins->Ntotal, 0.0, ins->o_PI);
  //if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, ins->o_PI);

  ins->NiterP = ellipticSolve(solver, ins->presTOL, ins->o_wrk3, ins->o_PI);

  ins->pressureAddBCKernel(mesh->Nelements,
                           time,
                           ins->dt,
                           ins->fieldOffset,
                           mesh->o_sgeo,
                           mesh->o_x,
                           mesh->o_y,
                           mesh->o_z,
                           mesh->o_vmapM,
                           mesh->o_EToB,
                           ins->o_EToB,
                           ins->o_usrwrk,
                           ins->o_U,
                           ins->o_P,
                           ins->o_PI);

  // update (increment) all points but not Dirichlet
  ins->pressureUpdateKernel(mesh->Nelements,
                            ins->fieldOffset,
                            solver->o_mapB,
                            ins->o_PI,
                            ins->o_P,
                            ins->o_wrk0);
  return ins->o_wrk0;
}


occa::memory velocitySolve(ins_t *ins, dfloat time) 
{
  mesh_t *mesh = ins->mesh;

  ins->pqKernel(
       mesh->Nelements*mesh->Np,
       ins->o_mue,
       ins->o_div,
       ins->o_P,
       ins->o_wrk3); 
  ins->gradientVolumeKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_Dmatrices,
       ins->fieldOffset,
       ins->o_wrk3,
       ins->o_wrk0);

  ins->velocityRhsKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_MM,
       ins->idt,
       ins->o_extbdfA,
       ins->o_extbdfB,
       ins->fieldOffset,
       ins->o_U,
       ins->o_BF,
       ins->o_FU,
       ins->o_wrk0,
       ins->o_rho,
       ins->o_wrk3);

  ins->velocityRhsBCKernel(
       mesh->Nelements,
       ins->fieldOffset,
       mesh->o_ggeo,
       mesh->o_sgeo,
       mesh->o_Dmatrices,
       mesh->o_Smatrices,
       mesh->o_MM,
       mesh->o_vmapM,
       mesh->o_EToB,
       mesh->o_sMT,
       time,
       mesh->o_x,
       mesh->o_y,
       mesh->o_z,
       ins->o_VmapB,
       ins->o_EToB,
       ins->o_usrwrk,
       ins->o_U,
       ins->o_ellipticCoeff,
       ins->o_wrk3);

  ogsGatherScatterMany(ins->o_wrk3, ins->NVfields, ins->fieldOffset,
                       ogsDfloat, ogsAdd, mesh->ogs);


  // Use old velocity as initial condition
  ins->o_wrk0.copyFrom(ins->o_U, ins->NVfields*ins->fieldOffset*sizeof(dfloat));

  if(ins->uvwSolver){
    if (ins->uvwSolver->Nmasked) mesh->maskKernel(ins->uvwSolver->Nmasked, ins->uvwSolver->o_maskIds, ins->o_wrk0);
    if (ins->uvwSolver->Nmasked) mesh->maskKernel(ins->uvwSolver->Nmasked, ins->uvwSolver->o_maskIds, ins->o_wrk3);
    ins->NiterU = ellipticSolve(ins->uvwSolver, ins->velTOL, ins->o_wrk3, ins->o_wrk0);
  } else {
    if (ins->uSolver->Nmasked) mesh->maskKernel(ins->uSolver->Nmasked, ins->uSolver->o_maskIds, ins->o_wrk0);
    if (ins->vSolver->Nmasked) mesh->maskKernel(ins->vSolver->Nmasked, ins->vSolver->o_maskIds, ins->o_wrk1);
    if (ins->wSolver->Nmasked) mesh->maskKernel(ins->wSolver->Nmasked, ins->wSolver->o_maskIds, ins->o_wrk2);
    if (ins->uSolver->Nmasked) mesh->maskKernel(ins->uSolver->Nmasked, ins->uSolver->o_maskIds, ins->o_wrk3);
    if (ins->vSolver->Nmasked) mesh->maskKernel(ins->vSolver->Nmasked, ins->vSolver->o_maskIds, ins->o_wrk4);
    if (ins->wSolver->Nmasked) mesh->maskKernel(ins->wSolver->Nmasked, ins->wSolver->o_maskIds, ins->o_wrk5);
    ins->NiterU = ellipticSolve(ins->uSolver, ins->velTOL, ins->o_wrk3, ins->o_wrk0);
    ins->NiterV = ellipticSolve(ins->vSolver, ins->velTOL, ins->o_wrk4, ins->o_wrk1);
    ins->NiterW = ellipticSolve(ins->wSolver, ins->velTOL, ins->o_wrk5, ins->o_wrk2);
  } 
 
  ins->velocityAddBCKernel(mesh->Nelements,
                           ins->fieldOffset,
                           time,
                           mesh->o_sgeo,
                           mesh->o_x,
                           mesh->o_y,
                           mesh->o_z,
                           mesh->o_vmapM,
                           mesh->o_EToB,
                           ins->o_EToB,
                           ins->o_usrwrk,
                           ins->o_U,
                           ins->o_wrk0);

  return ins->o_wrk0;
}

} // namespace

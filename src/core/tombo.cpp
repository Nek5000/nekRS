#include "nrs.hpp"
#include "udf.hpp"

namespace tombo {

occa::memory pressureSolve(ins_t *ins, dfloat time)
{
  mesh_t *mesh = ins->mesh;

  ins->curlKernel(mesh->Nelements,
                  mesh->o_vgeo,
                  mesh->o_Dmatrices,
                  ins->fieldOffset,
                  ins->o_Ue,
                  ins->o_wrk0);

  ogsGatherScatterMany(ins->o_wrk0, ins->NVfields, ins->fieldOffset,
                       ogsDfloat, ogsAdd, mesh->ogs);

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
       ins->fieldOffset,
       ins->o_wrk0,
       ins->o_wrk3);

  if(ins->options.compareArgs("VARIABLE VISCOSITY", "TRUE"))
    ins->pressureStressKernel( 
         mesh->Nelements,
         mesh->o_vgeo,
         mesh->o_Dmatrices,
         ins->fieldOffset,
         ins->o_mue,
         ins->o_Ue,
         ins->o_div,
         ins->o_wrk3);

  ins->gradientVolumeKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_Dmatrices,
       ins->fieldOffset,
       ins->o_div,
       ins->o_wrk0);

  occa::memory o_irho = ins->o_ellipticCoeff;
  ins->pressureRhsKernel(
       mesh->Nelements*mesh->Np,
       ins->fieldOffset,
       ins->o_mue,
       o_irho,
       ins->o_BF,
       ins->o_wrk3,
       ins->o_wrk0,
       ins->o_wrk6);

  ogsGatherScatterMany(ins->o_wrk6, ins->NVfields, ins->fieldOffset,
                       ogsDfloat, ogsAdd, mesh->ogs);

  ins->invMassMatrixKernel(
       mesh->Nelements,
       ins->fieldOffset,
       ins->NVfields,
       mesh->o_vgeo,
       ins->o_InvM,
       ins->o_wrk6);

  ins->divergenceVolumeKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_Dmatrices,
       ins->fieldOffset,
       ins->o_wrk6,
       ins->o_wrk3);

  const dfloat lambda = ins->g0*ins->idt;
  ins->divergenceSurfaceKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_sgeo,
       mesh->o_vmapM,
       mesh->o_EToB,
       ins->o_EToB,
       time,
       -lambda, 
       mesh->o_x,
       mesh->o_y,
       mesh->o_z,
       ins->fieldOffset,
       ins->o_usrwrk,
       ins->o_U,
       ins->o_wrk3);

  ins->pressureAxKernel(
       mesh->Nelements,
       ins->fieldOffset,
       mesh->o_ggeo,
       mesh->o_Dmatrices,
       mesh->o_Smatrices,
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

  ogsGatherScatter(ins->o_wrk3, ogsDfloat, ogsAdd, mesh->ogs);
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

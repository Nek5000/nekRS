#include "nrs.hpp"
#include "udf.hpp"

void curlCurl(ins_t *ins, occa::memory o_wrk, occa::memory o_U, 
              occa::memory o_NC);

namespace tombo {

void pressureSolve(ins_t *ins, dfloat time, occa::memory o_wrk, occa::memory o_Pnew)
{
  mesh_t *mesh = ins->mesh;
  occa::memory o_wrk1 = o_wrk;
  occa::memory o_wrk2 = o_wrk.slice(3*ins->fieldOffset*sizeof(dfloat));

  // o_wrk2 = nu*( JW*curl(curl(v)) - 4/3*JW*grad(qtl) ) 
  curlCurl(ins, o_wrk1, ins->o_Ue, o_wrk2);
  ins->gradientVolumeKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_Dmatrices,
       ins->fieldOffset,
       ins->o_qtl,
       o_wrk1);
  occa::memory o_irho = ins->o_ellipticCoeff.slice(0*ins->fieldOffset*sizeof(dfloat));
  ins->ncKernel(
       mesh->Np*mesh->Nelements,
       ins->fieldOffset,
       ins->o_mue,
       o_irho,
       o_wrk1,
       o_wrk2);

  // o_wrk1 = JW*o_BF - o_wrk2
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
       o_wrk2,
       ins->o_FU,
       o_wrk1);

  ogsGatherScatterMany(o_wrk1, ins->NVfields, ins->fieldOffset,
                       ogsDfloat, ogsAdd, mesh->ogs);

  ins->invMassMatrixKernel(
       mesh->Nelements,
       ins->fieldOffset,
       ins->NVfields,
       mesh->o_vgeo,
       ins->o_InvM,
       o_wrk1);

  // o_wrk2 = div(o_wrk)
  ins->divergenceVolumeKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_Dmatrices,
       ins->fieldOffset,
       o_wrk1,
       o_wrk2);

  // o_wrk2 += surface term
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
       o_wrk1,
       o_wrk2);


  // we solve for delta p => o_wrk2 += -(grad P, grad phi)
  ins->AxKernel(
       mesh->Nelements,
       ins->fieldOffset,
       mesh->o_ggeo,
       mesh->o_Dmatrices,
       mesh->o_Smatrices,
       mesh->o_MM,
       ins->o_P,
       ins->o_ellipticCoeff,
       o_wrk2);

  // o_wrk += g0/dt*qtl
  ins->pressureAddQtlKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       lambda,
       ins->o_qtl,
       o_wrk2);

  elliptic_t *solver = ins->pSolver;

  ogsGatherScatter(o_wrk2, ogsDfloat, ogsAdd, mesh->ogs);
  if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, o_wrk2);

  ins->setScalarKernel(ins->Ntotal, 0.0, ins->o_PI);
  if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, ins->o_PI);

  ins->NiterP = ellipticSolve(solver, 0.0, ins->presTOL, o_wrk2, ins->o_PI);

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
                            o_Pnew);
}


void velocitySolve(ins_t *ins, dfloat time, occa::memory o_wrk1, occa::memory o_Unew) 
{
  mesh_t *mesh = ins->mesh;

  occa::memory o_wrk2 = o_Unew;

  // o_wrk2 = grad(p - 1/3*mue*qtl)
  ins->pqKernel(
       mesh->Nelements*mesh->Np,
       ins->o_mue,
       ins->o_qtl,
       ins->o_P,
       o_wrk1); 
  ins->gradientVolumeKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_Dmatrices,
       ins->fieldOffset,
       o_wrk1,
       o_wrk2);

  // o_wrk1 = o_rho*o_BF - o_wrk2 
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
       o_wrk2,
       ins->o_rho,
       o_wrk1);

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
       o_wrk1);

  ogsGatherScatterMany(o_wrk1, ins->NVfields, ins->fieldOffset,
                       ogsDfloat, ogsAdd, mesh->ogs);
  occa::memory o_wrk1u = o_wrk1.slice(0*ins->fieldOffset*sizeof(dfloat));
  occa::memory o_wrk1v = o_wrk1.slice(1*ins->fieldOffset*sizeof(dfloat));
  occa::memory o_wrk1w = o_wrk1.slice(2*ins->fieldOffset*sizeof(dfloat));
 
  // Use old velocity for velocity solver initial condition
  o_Unew.copyFrom(ins->o_U,ins->NVfields*ins->fieldOffset*sizeof(dfloat));
  occa::memory o_u = o_Unew.slice(0*ins->fieldOffset*sizeof(dfloat));
  occa::memory o_v = o_Unew.slice(1*ins->fieldOffset*sizeof(dfloat));
  occa::memory o_w = o_Unew.slice(2*ins->fieldOffset*sizeof(dfloat));

  if (ins->uSolver->Nmasked) mesh->maskKernel(ins->uSolver->Nmasked, ins->uSolver->o_maskIds, o_u);
  if (ins->vSolver->Nmasked) mesh->maskKernel(ins->vSolver->Nmasked, ins->vSolver->o_maskIds, o_v);
  if (ins->wSolver->Nmasked) mesh->maskKernel(ins->wSolver->Nmasked, ins->wSolver->o_maskIds, o_w);
  
  if (ins->uSolver->Nmasked) mesh->maskKernel(ins->uSolver->Nmasked, ins->uSolver->o_maskIds, o_wrk1u);
  if (ins->vSolver->Nmasked) mesh->maskKernel(ins->vSolver->Nmasked, ins->vSolver->o_maskIds, o_wrk1v);
  if (ins->wSolver->Nmasked) mesh->maskKernel(ins->wSolver->Nmasked, ins->wSolver->o_maskIds, o_wrk1w);
   
  const dfloat lambda = 1; // dummy
  ins->NiterU = ellipticSolve(ins->uSolver, lambda, ins->velTOL, o_wrk1u, o_u);
  ins->NiterV = ellipticSolve(ins->vSolver, lambda, ins->velTOL, o_wrk1v, o_v);
  ins->NiterW = ellipticSolve(ins->wSolver, lambda, ins->velTOL, o_wrk1w, o_w);
  
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
                           o_Unew);
}


} // namespace

void curlCurl(ins_t *ins, occa::memory o_wrk, occa::memory o_U, 
              occa::memory o_NC)
{
  mesh_t *mesh = ins->mesh;

  ins->curlKernel(mesh->Nelements,
                  mesh->o_vgeo,
                  mesh->o_Dmatrices,
                  mesh->o_MM, 
                  ins->fieldOffset,
                  o_U,
                  o_wrk);

  ogsGatherScatterMany(o_wrk, ins->NVfields, ins->fieldOffset,
                       ogsDfloat, ogsAdd, mesh->ogs);

  ins->invMassMatrixKernel(
       mesh->Nelements,
       ins->fieldOffset,
       ins->NVfields,
       mesh->o_vgeo,
       ins->o_InvM,
       o_wrk);

  ins->curlKernel(mesh->Nelements,
                  mesh->o_vgeo,
                  mesh->o_Dmatrices,
                  mesh->o_MM, 
                  ins->fieldOffset,
                  o_wrk,
                  o_NC);

/*
  ogsGatherScatterMany(o_NC, ins->NVfields, ins->fieldOffset,
                       ogsDfloat, ogsAdd, mesh->ogs);

  ins->invMassMatrixKernel(
       mesh->Nelements,
       ins->fieldOffset,
       ins->NVfields,
       mesh->o_vgeo,
       ins->o_InvM,
       o_NC);
*/
}

#include "nrs.hpp"

void insGradient(ins_t *ins, dfloat time, occa::memory o_U, occa::memory o_GU);
void curlCurl(ins_t *ins, dfloat time, occa::memory o_wrk, occa::memory o_U, 
              occa::memory o_NC);

namespace tombo {

void pressureRhs(ins_t *ins, dfloat time, occa::memory o_rhsP)
{
  mesh_t *mesh = ins->mesh;
  occa::memory o_wrk = ins->o_scratch;

  ins->setEllipticCoeffPressureKernel(
       mesh->Np*mesh->Nelements,
       ins->fieldOffset,
       ins->o_rho,
       ins->o_ellipticCoeff);

  // o_NC = nu*( curl(curl(v)) - 4/3*grad(qtl) ) 
  insGradient(ins, time, ins->o_qtl, o_wrk);
  occa::memory o_irho = ins->o_ellipticCoeff.slice(0*ins->fieldOffset*sizeof(dfloat));
  occa::memory o_NC   = ins->o_scratch.slice(3*ins->fieldOffset*sizeof(dfloat));
  curlCurl(ins, time, o_wrk, ins->o_Ue, o_NC);
  ins->ncKernel(
       mesh->Np*mesh->Nelements,
       ins->fieldOffset,
       ins->o_mue,
       o_irho,
       o_wrk,
       o_NC);

  // o_wrk = f - o_NC
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
       ins->o_NU,
       o_NC,
       ins->o_FU,
       o_wrk);

  // TODO: dssum o_wrk

  // weak divergence term (no boundary contribution) 
  ins->divergenceVolumeKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_Dmatrices,
       ins->fieldOffset,
       o_wrk,
       o_rhsP);

  // we solve for delta p => o_rhsP -= (grad P, grad phi)
  ins->AxKernel(
       mesh->Nelements,
       mesh->o_ggeo,
       mesh->o_Dmatrices,
       mesh->o_Smatrices,
       mesh->o_MM,
       ins->o_P,
       ins->o_ellipticCoeff,
       o_rhsP);
 
  // o_rhsP += lambda*qtl
  const dfloat lambda = ins->g0*ins->idt; 
  ins->scaledAddKernel(
       mesh->Nelements*mesh->Np,
       lambda,
       ins->qtl,
       1.0,
       o_rhsP);

  // o_rho += surface term
  ins->divergenceSurfaceKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_sgeo,
       mesh->o_LIFTT,
       mesh->o_vmapM,
       mesh->o_vmapP,
       ins->o_EToB,
       time,
       lambda, 
       mesh->o_x,
       mesh->o_y,
       mesh->o_z,
       ins->fieldOffset,
       ins->o_usrwrk,
       o_wrk,
       o_rhsP);
}


void pressureSolve(ins_t *ins, dfloat time, occa::memory o_rhsP, occa::memory o_rkP)
{
  mesh_t *mesh = ins->mesh;
  elliptic_t *solver = ins->pSolver;

  ogsGatherScatter(o_rhsP, ogsDfloat, ogsAdd, mesh->ogs);
  if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, o_rhsP);
  if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, ins->o_PI);

  ins->NiterP = ellipticSolve(solver, 0.0, ins->presTOL, o_rhsP, ins->o_PI);

  ins->pressureAddBCKernel(mesh->Nelements,
                           time,
                           ins->dt,
                           ins->fieldOffset,
                           mesh->o_sgeo,
                           mesh->o_x,
                           mesh->o_y,
                           mesh->o_z,
                           mesh->o_vmapM,
                           ins->o_PmapB,
                           ins->o_EToB,
                           ins->o_usrwrk,
                           ins->o_U,
                           ins->o_PI);

  ins->pressureUpdateKernel(mesh->Nelements,
                            ins->fieldOffset,
                            ins->o_PmapB,
                            ins->o_PI,
                            ins->o_P,
                            o_rkP);
}


void velocityRhs(ins_t *ins, dfloat time, occa::memory o_rhsU) 
{
  mesh_t *mesh = ins->mesh;
  occa::memory o_wrk = ins->o_scratch;
  occa::memory o_PQ  = ins->o_scratch.slice(3*ins->fieldOffset*sizeof(dfloat));

  // grad(p - 1/3*mue*qtl)
  ins->pqKernel(
       mesh->Nelements*mesh->Np,
       ins->o_mue,
       ins->o_qtl,
       ins->o_P,
       o_PQ);  
  insGradient (ins, time, o_PQ, o_wrk);

  ins->velocityRhsKernel(mesh->Nelements,
                         mesh->o_vgeo,
                         mesh->o_MM,
                         ins->idt,
                         ins->o_extbdfA,
                         ins->o_extbdfB,
                         ins->fieldOffset,
                         ins->o_U,
                         ins->o_NU,
                         ins->o_FU,
                         o_wrk,
                         ins->o_rho,
                         o_rhsU);
}

void velocitySolve(ins_t *ins, dfloat time, occa::memory o_rhsU, occa::memory o_UH) 
{
  mesh_t *mesh = ins->mesh;
  dlong Ntotal = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;

  elliptic_t *usolver = ins->uSolver;
  elliptic_t *vsolver = ins->vSolver;
  elliptic_t *wsolver = ins->wSolver;

  occa::memory o_rhsu = o_rhsU.slice(0*ins->fieldOffset*sizeof(dfloat));
  occa::memory o_rhsv = o_rhsU.slice(1*ins->fieldOffset*sizeof(dfloat));
  occa::memory o_rhsw = o_rhsU.slice(2*ins->fieldOffset*sizeof(dfloat));

  occa::memory o_uh = o_UH.slice(0*ins->fieldOffset*sizeof(dfloat));
  occa::memory o_vh = o_UH.slice(1*ins->fieldOffset*sizeof(dfloat));
  occa::memory o_wh = o_UH.slice(2*ins->fieldOffset*sizeof(dfloat));

  ins->velocityRhsBCKernel(
       mesh->Nelements,
       ins->fieldOffset,
       mesh->o_ggeo,
       mesh->o_sgeo,
       mesh->o_Dmatrices,
       mesh->o_Smatrices,
       mesh->o_MM,
       mesh->o_vmapM,
       ins->o_EToB,
       mesh->o_sMT,
       time,
       mesh->o_x,
       mesh->o_y,
       mesh->o_z,
       ins->o_VmapB,
       ins->o_usrwrk,
       ins->o_U,
       ins->o_ellipticCoeff,
       o_rhsU);

  ogsGatherScatterMany(o_rhsU, ins->NVfields, ins->fieldOffset,
                       ogsDfloat, ogsAdd, mesh->ogs);

  // Use old velocity for velocity solver initial condition
  o_UH.copyFrom(ins->o_U,ins->NVfields*Ntotal*sizeof(dfloat));

  // TODO: fuse with rhs into single kernel 
  if (usolver->Nmasked) mesh->maskKernel(usolver->Nmasked, usolver->o_maskIds, o_uh);
  if (vsolver->Nmasked) mesh->maskKernel(vsolver->Nmasked, vsolver->o_maskIds, o_vh);
  if (ins->dim==3 && wsolver->Nmasked)
    mesh->maskKernel(wsolver->Nmasked, wsolver->o_maskIds, o_wh);

  if (usolver->Nmasked) mesh->maskKernel(usolver->Nmasked, usolver->o_maskIds, o_rhsu);
  if (vsolver->Nmasked) mesh->maskKernel(vsolver->Nmasked, vsolver->o_maskIds, o_rhsv);
  if (ins->dim==3 && wsolver->Nmasked)
    mesh->maskKernel(wsolver->Nmasked, wsolver->o_maskIds, o_rhsw);

  const dfloat lambda = 1; // dummy
  ins->NiterU = ellipticSolve(usolver, lambda, ins->velTOL, o_rhsu, o_uh);
  ins->NiterV = ellipticSolve(vsolver, lambda, ins->velTOL, o_rhsv, o_vh);
  if (ins->dim==3) 
    ins->NiterW = ellipticSolve(wsolver, lambda, ins->velTOL, o_rhsw, o_wh);
  
  ins->velocityAddBCKernel(mesh->Nelements,
                           ins->fieldOffset,
                           time,
                           mesh->o_sgeo,
                           mesh->o_x,
                           mesh->o_y,
                           mesh->o_z,
                           mesh->o_vmapM,
                           ins->o_VmapB,
                           ins->o_usrwrk,
                           ins->o_U,
                           o_UH);
}

} // namespace

void curlCurl(ins_t *ins, dfloat time, occa::memory o_wrk, occa::memory o_U, 
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
       ins->o_InvM, // mesh->o_MM, // should be invMM for tri/tet
       o_wrk);

  ins->curlKernel(mesh->Nelements,
                  mesh->o_vgeo,
                  mesh->o_Dmatrices,
                  mesh->o_MM, 
                  ins->fieldOffset,
                  o_wrk,
                  o_NC);

  ogsGatherScatterMany(o_NC, ins->NVfields, ins->fieldOffset,
                       ogsDfloat, ogsAdd, mesh->ogs);

  ins->invMassMatrixKernel(
       mesh->Nelements,
       ins->fieldOffset,
       ins->NVfields,
       mesh->o_vgeo,
       ins->o_InvM, // mesh->o_MM, // should be invMM for tri/tet
       o_NC);
}

void insGradient(ins_t *ins, dfloat time, occa::memory o_U, occa::memory o_GU)
{
  mesh_t *mesh = ins->mesh;

  ins->gradientVolumeKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_Dmatrices,
       ins->fieldOffset,
       o_U,
       o_GU);
}

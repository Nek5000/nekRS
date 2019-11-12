#include "nrs.hpp"
#include "udf.hpp"

void curlCurl(ins_t *ins, occa::memory o_wrk, occa::memory o_U, 
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

  // o_NC = nu*( JW*curl(curl(v)) - 4/3*JW*grad(qtl) ) 
  occa::memory o_NC = ins->o_scratch.slice(3*ins->fieldOffset*sizeof(dfloat));
  curlCurl(ins, o_wrk, ins->o_Ue, o_NC);
  ins->gradientVolumeKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_Dmatrices,
       ins->fieldOffset,
       ins->o_qtl,
       o_wrk);
  occa::memory o_irho = ins->o_ellipticCoeff.slice(0*ins->fieldOffset*sizeof(dfloat));
  ins->ncKernel(
       mesh->Np*mesh->Nelements,
       ins->fieldOffset,
       ins->o_mue,
       o_irho,
       o_wrk,
       o_NC);

  // o_wrk = JW*f - o_NC
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

  ogsGatherScatterMany(o_wrk, ins->NVfields, ins->fieldOffset,
                       ogsDfloat, ogsAdd, mesh->ogs);

  ins->invMassMatrixKernel(
       mesh->Nelements,
       ins->fieldOffset,
       ins->NVfields,
       mesh->o_vgeo,
       ins->o_InvM, // mesh->o_MM, // should be invMM for tri/tet
       o_wrk);

  // o_rhsP = div(o_wrk)
  ins->divergenceVolumeKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_Dmatrices,
       ins->fieldOffset,
       o_wrk,
       o_rhsP);

  // o_rhsP += surface term
  const dfloat lambda = ins->g0*ins->idt;
  ins->divergenceSurfaceKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_sgeo,
       mesh->o_LIFTT,
       mesh->o_vmapM,
       mesh->o_vmapP,
       ins->o_EToB,
       time,
       -lambda, 
       mesh->o_x,
       mesh->o_y,
       mesh->o_z,
       ins->fieldOffset,
       ins->o_usrwrk,
       o_wrk,
       o_rhsP);


  // we solve for delta p => o_rhsP += -(grad P, grad phi)
  ins->AxKernel(
       mesh->Nelements,
       ins->fieldOffset,
       mesh->o_ggeo,
       mesh->o_Dmatrices,
       mesh->o_Smatrices,
       mesh->o_MM,
       ins->o_P,
       ins->o_ellipticCoeff,
       o_rhsP);

  // o_rhsP += g0/dt*qtl
  ins->pressureAddQtlKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       lambda,
       ins->o_qtl,
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
  ins->gradientVolumeKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_Dmatrices,
       ins->fieldOffset,
       o_PQ,
       o_wrk);

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

  ins->setEllipticCoeffKernel(
       mesh->Np*mesh->Nelements,
       ins->g0*ins->idt,
       ins->fieldOffset,
       ins->o_mue,
       ins->o_rho,
       ins->o_ellipticCoeff);

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

// qtl = 1/(rho*cp*T) * (div[k*grad[T] ] + qvol)
void qthermal(ins_t *ins, dfloat time, occa::memory o_qtl)
{
  cds_t *cds = ins->cds;
  mesh_t *mesh = ins->mesh;

  occa::memory o_gradS = ins->o_scratch;
  
  ins->gradientVolumeKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_Dmatrices,
       ins->fieldOffset,
       cds->o_S,
       o_gradS);

  ogsGatherScatterMany(o_gradS, ins->NVfields, ins->fieldOffset,
                       ogsDfloat, ogsAdd, mesh->ogs);

   ins->invMassMatrixKernel(
        mesh->Nelements,
        ins->fieldOffset,
        ins->NVfields,
        mesh->o_vgeo,
        ins->o_InvM, 
        o_gradS);

   occa::memory o_src = cds->o_rkS;
   if(udf.sEqnSource)
     udf.sEqnSource(ins, time, cds->o_S, o_src);
   else
     ins->setScalarKernel(mesh->Nelements*mesh->Np, 0.0, o_src);

   ins->qtlKernel(
        mesh->Nelements,
        mesh->o_vgeo,
        mesh->o_Dmatrices,
        ins->fieldOffset,
        o_gradS,
        cds->o_S,
        cds->o_diff,
        cds->o_rho,
        o_src,
        o_qtl);

  ogsGatherScatterMany(o_qtl, 1, ins->fieldOffset,
                       ogsDfloat, ogsAdd, mesh->ogs);

   ins->invMassMatrixKernel(
        mesh->Nelements,
        ins->fieldOffset,
        1,
        mesh->o_vgeo,
        ins->o_InvM, 
        o_qtl);
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
       ins->o_InvM, // mesh->o_MM, // should be invMM for tri/tet
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
       ins->o_InvM, // mesh->o_MM, // should be invMM for tri/tet
       o_NC);
*/
}

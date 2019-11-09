#include "cds.h"

void cdsHelmholtzSolve(cds_t *cds, dfloat time, int stage,occa::memory o_rhsS,occa::memory o_Shat){
  
  mesh_t     *mesh   = cds->mesh; 
  elliptic_t *solver = cds->solver;

  cds->setEllipticCoeffKernel(
       mesh->Np*mesh->Nelements,
       cds->g0*cds->idt,
       cds->sOffset,
       cds->o_diff,
       cds->o_rho,
       cds->o_ellipticCoeff);

  cds->helmholtzRhsBCKernel(mesh->Nelements,
                            mesh->o_ggeo,
                            mesh->o_sgeo,
                            mesh->o_Dmatrices,
                            mesh->o_Smatrices,
                            mesh->o_MM,
                            mesh->o_vmapM,
                            cds->o_EToB,
                            mesh->o_sMT,
                            time,
                            cds->sOffset,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            cds->o_mapB,
                            cds->o_ellipticCoeff,
                            o_rhsS);

  ogsGatherScatter(o_rhsS, ogsDfloat, ogsAdd, mesh->ogs);
  if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, o_rhsS);

  //copy current solution fields as initial guess
  dlong Ntotal = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;
  o_Shat.copyFrom(cds->o_S,Ntotal*sizeof(dfloat),0,0*cds->sOffset*sizeof(dfloat)); 
 
  if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, o_Shat);

  const dfloat lambda = 1; // dummy value not used if coeff is variable
  cds->Niter = ellipticSolve(solver, lambda, cds->TOL, o_rhsS, o_Shat);

  cds->helmholtzAddBCKernel(mesh->Nelements,
                            time,
                            mesh->o_sgeo,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            mesh->o_vmapM,
                            cds->o_mapB,
                            o_Shat);
}

void cdsHelmholtzRhs(cds_t *cds, dfloat time, int stage, occa::memory o_rhsS){
  
  mesh_t *mesh = cds->mesh; 

  cds->helmholtzRhsKernel(mesh->Nelements,
                          mesh->o_vgeo,
                          mesh->o_MM,
                          cds->idt,
                          cds->o_extbdfA,
                          cds->o_extbdfB,
                          cds->o_extbdfC,
                          cds->sOffset,
                          cds->o_S,
                          cds->o_NS,
                          cds->o_FS,
                          cds->o_rho,
                          o_rhsS);
}

void cdsAdvection(cds_t *cds, dfloat time, occa::memory o_U, occa::memory o_S, occa::memory o_NS){

  mesh_t *mesh = cds->mesh;

  if(cds->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
    cds->advectionStrongCubatureVolumeKernel(
           cds->meshV->Nelements,
           mesh->o_vgeo,
           mesh->o_cubvgeo,
           mesh->o_cubDiffInterpT, // mesh->o_cubDWmatrices,
           mesh->o_cubInterpT,
           mesh->o_cubProjectT,
           cds->vOffset,
           cds->sOffset,
           o_U,
           o_S,
           cds->o_rho,
           o_NS);
  else
    cds->advectionStrongVolumeKernel(
           cds->meshV->Nelements,
           mesh->o_vgeo,
           mesh->o_Dmatrices,
           cds->vOffset,
           cds->sOffset,
           o_U,
           o_S,
           cds->o_rho,
           o_NS);

}

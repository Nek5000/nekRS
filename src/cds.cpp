#include "cds.h"

void cdsSolve(cds_t *cds, dfloat time, occa::memory o_wrk, occa::memory o_Shat){
  
  mesh_t     *mesh   = cds->mesh; 
  elliptic_t *solver = cds->solver;

  o_wrk.copyFrom(cds->o_BF, cds->Ntotal*sizeof(dfloat)); 

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
                            cds->fieldOffset,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            cds->o_mapB,
                            cds->o_ellipticCoeff,
                            o_wrk);

  ogsGatherScatter(o_wrk, ogsDfloat, ogsAdd, mesh->ogs);
  if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, o_wrk);

  //copy current solution fields as initial guess
  o_Shat.copyFrom(cds->o_S, cds->Ntotal*sizeof(dfloat)); 
 
  if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, o_Shat);

  const dfloat lambda = 1; // dummy value not used if coeff is variable
  cds->Niter = ellipticSolve(solver, lambda, cds->TOL, o_wrk, o_Shat);

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
           cds->vFieldOffset,
           cds->fieldOffset,
           o_U,
           o_S,
           cds->o_rho,
           o_NS);
  else
    cds->advectionStrongVolumeKernel(
           cds->meshV->Nelements,
           mesh->o_vgeo,
           mesh->o_Dmatrices,
           cds->vFieldOffset,
           cds->fieldOffset,
           o_U,
           o_S,
           cds->o_rho,
           o_NS);

}

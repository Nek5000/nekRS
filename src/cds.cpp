#include "nrs.hpp"

void cdsSolve(const int is, cds_t *cds, dfloat time, occa::memory o_wrk, occa::memory o_Shat){

  mesh_t *mesh; 
  (is) ? mesh = cds->meshV : mesh = cds->mesh;
  elliptic_t *solver = cds->solver[is];

  occa::memory o_BFi = cds->o_BF.slice(is*cds->fieldOffset*sizeof(dfloat));
  o_wrk.copyFrom(o_BFi, cds->Ntotal*sizeof(dfloat)); 

  cds->helmholtzRhsBCKernel(mesh->Nelements,
                            mesh->o_ggeo,
                            mesh->o_sgeo,
                            mesh->o_Dmatrices,
                            mesh->o_Smatrices,
                            mesh->o_MM,
                            mesh->o_vmapM,
                            mesh->o_EToB,
                            mesh->o_sMT,
                            is,
                            time,
                            cds->fieldOffset,
                            solver->Ntotal, // lambda offset required by elliptic
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            cds->o_mapB[is],
                            cds->o_ellipticCoeff,
                            cds->o_usrwrk,
                            o_wrk);

  ogsGatherScatter(o_wrk, ogsDfloat, ogsAdd, mesh->ogs);
  if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, o_wrk);

  //copy current solution fields as initial guess
  occa::memory o_Si = cds->o_S.slice(is*cds->fieldOffset*sizeof(dfloat));
  o_Shat.copyFrom(o_Si, cds->Ntotal*sizeof(dfloat)); 
  if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, o_Shat);

  const dfloat lambda = 1; // dummy value not used if coeff is variable
  cds->Niter[is] = ellipticSolve(solver, lambda, cds->TOL, o_wrk, o_Shat);

  cds->helmholtzAddBCKernel(mesh->Nelements,
                            cds->fieldOffset,
                            is,
                            time,
                            mesh->o_sgeo,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            mesh->o_vmapM,
                            mesh->o_EToB,
                            cds->o_EToB[is],
                            cds->o_usrwrk,
                            o_Shat);
}

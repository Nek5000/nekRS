#include "nrs.hpp"

occa::memory cdsSolve(const int is, cds_t* cds, dfloat time)
{
  mesh_t* mesh;
  (is) ? mesh = cds->meshV : mesh = cds->mesh;
  elliptic_t* solver = cds->solver[is];

  cds->o_wrk1.copyFrom(cds->o_BF, cds->Ntotal * sizeof(dfloat), 0,
                       is * cds->fieldOffset * sizeof(dfloat));

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
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            cds->o_mapB[is],
                            cds->o_ellipticCoeff,
                            *(cds->o_usrwrk),
                            cds->o_wrk1);

  ogsGatherScatter(cds->o_wrk1, ogsDfloat, ogsAdd, mesh->ogs);
  if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, cds->o_wrk1);

  //copy current solution fields as initial guess
  cds->o_wrk0.copyFrom(cds->o_S, cds->Ntotal * sizeof(dfloat), 0,
                       is * cds->fieldOffset * sizeof(dfloat));
  if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, cds->o_wrk0);

  cds->Niter[is] = ellipticSolve(solver, cds->TOL, cds->o_wrk1, cds->o_wrk0);

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
                            *(cds->o_usrwrk),
                            cds->o_wrk0);
  return cds->o_wrk0;
}

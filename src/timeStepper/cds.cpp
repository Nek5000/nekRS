#include <limits>
#include "nrs.hpp"

occa::memory cdsSolve(const int is, cds_t* cds, dfloat time)
{
  mesh_t* mesh;
  oogs_t* gsh;
  if(is) {
    mesh = cds->meshV;
    gsh = cds->gsh;
  } else {
    mesh = cds->mesh;
    gsh = cds->gshT;
  }
  elliptic_t* solver = cds->solver[is];

  cds->o_wrk0.copyFrom(cds->o_S, cds->Ntotal * sizeof(dfloat), 0, is * cds->fieldOffset * sizeof(dfloat));

  //enforce Dirichlet BCs
  cds->fillKernel(cds->fieldOffset, std::numeric_limits<dfloat>::min(), cds->o_wrk2); 
  for (int sweep = 0; sweep < 2; sweep++) {
    cds->dirichletBCKernel(mesh->Nelements,
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
                           cds->o_wrk2);

    //take care of Neumann-Dirichlet shared edges across elements
    if(sweep == 0) oogs::startFinish(cds->o_wrk2, 1, cds->fieldOffset, ogsDfloat, ogsMax, gsh);
    if(sweep == 1) oogs::startFinish(cds->o_wrk2, 1, cds->fieldOffset, ogsDfloat, ogsMin, gsh);
  }
  if (solver->Nmasked) cds->maskCopyKernel(solver->Nmasked, 0, solver->o_maskIds, cds->o_wrk2, cds->o_wrk0);

  //build RHS
  cds->o_wrk1.copyFrom(cds->o_BF, cds->Ntotal * sizeof(dfloat), 0, is * cds->fieldOffset * sizeof(dfloat));
  cds->helmholtzRhsBCKernel(mesh->Nelements,
                            mesh->o_sgeo,
                            mesh->o_vmapM,
                            mesh->o_EToB,
                            mesh->o_sMT,
                            is,
                            time,
                            cds->fieldOffset,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            cds->o_wrk0,
                            cds->o_mapB[is],
                            *(cds->o_usrwrk),
                            cds->o_wrk1);
  oogs::startFinish(cds->o_wrk1, 1, cds->fieldOffset, ogsDfloat, ogsAdd, gsh);
  if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, cds->o_wrk1);

  if(cds->options.compareArgs("SCALAR INITIAL GUESS DEFAULT", "EXTRAPOLATION")) {
    cds->o_wrk0.copyFrom(cds->o_Se, cds->Ntotal * sizeof(dfloat), 0, is * cds->fieldOffset * sizeof(dfloat));
    if (solver->Nmasked) cds->maskCopyKernel(solver->Nmasked, 0, solver->o_maskIds, cds->o_wrk2, cds->o_wrk0);
  }

  cds->Niter[is] = ellipticSolve(solver, cds->TOL, cds->o_wrk1, cds->o_wrk0);

  return cds->o_wrk0;
}

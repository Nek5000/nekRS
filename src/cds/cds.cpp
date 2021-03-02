#include <limits>
#include "nrs.hpp"
#include "linAlg.hpp"

occa::memory cdsSolve(const int is, cds_t* cds, dfloat time)
{
  
  mesh_t* mesh;
  oogs_t* gsh;
  if(is) {
    mesh = cds->meshV;
    gsh = cds->gsh;
  } else {
    mesh = cds->meshT[0];
    gsh = cds->gshT;
  }
  elliptic_t* solver = cds->solver[is];

  platform->o_slice0.copyFrom(cds->o_S, cds->Ntotal * sizeof(dfloat), 0, is * cds->fieldOffset * sizeof(dfloat));

  //enforce Dirichlet BCs
  platform->linAlg->fill(cds->fieldOffset, std::numeric_limits<dfloat>::min(), platform->o_slice2); 
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
                           platform->o_slice2);

    //take care of Neumann-Dirichlet shared edges across elements
    if(sweep == 0) oogs::startFinish(platform->o_slice2, 1, cds->fieldOffset, ogsDfloat, ogsMax, gsh);
    if(sweep == 1) oogs::startFinish(platform->o_slice2, 1, cds->fieldOffset, ogsDfloat, ogsMin, gsh);
  }
  if (solver->Nmasked) cds->maskCopyKernel(solver->Nmasked, 0, solver->o_maskIds, platform->o_slice2, platform->o_slice0);

  //build RHS
  platform->o_slice1.copyFrom(cds->o_BF, cds->Ntotal * sizeof(dfloat), 0, is * cds->fieldOffset * sizeof(dfloat));
  cds->helmholtzRhsBCKernel(mesh->Nelements,
                            mesh->o_sgeo,
                            mesh->o_vmapM,
                            mesh->o_EToB,
                            is,
                            time,
                            cds->fieldOffset,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            platform->o_slice0,
                            cds->o_EToB[is],
                            cds->o_mapB[is],
                            *(cds->o_usrwrk),
                            platform->o_slice1);

  std::stringstream ss;
  ss << std::setfill('0') << std::setw(2) << is;
  string sid = ss.str();
  if(cds->options[is].compareArgs("SCALAR" + sid + " INITIAL GUESS DEFAULT", "EXTRAPOLATION")) {
    platform->o_slice0.copyFrom(cds->o_Se, cds->Ntotal * sizeof(dfloat), 0, is * cds->fieldOffset * sizeof(dfloat));
    if (solver->Nmasked) cds->maskCopyKernel(solver->Nmasked, 0, solver->o_maskIds, platform->o_slice2, platform->o_slice0);
  }
  ellipticSolve(solver, platform->o_slice1, platform->o_slice0);

  return platform->o_slice0;
}



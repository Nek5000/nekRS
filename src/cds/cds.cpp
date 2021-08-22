#include <limits>
#include "nrs.hpp"
#include "linAlg.hpp"

occa::memory cdsSolve(const int is, cds_t* cds, dfloat time, int stage)
{
  
  mesh_t* mesh;
  oogs_t* gsh;
  if(is) {
    mesh = cds->meshV;
    gsh = cds->gsh;
  } else {
    mesh = cds->mesh[0];
    gsh = cds->gshT;
  }
  elliptic_t* solver = cds->solver[is];

  platform->o_mempool.slice0.copyFrom(cds->o_S, cds->fieldOffset[is] * sizeof(dfloat), 0, cds->fieldOffsetScan[is] * sizeof(dfloat));

  //enforce Dirichlet BCs
  platform->linAlg->fill(cds->fieldOffset[is], -1.0*std::numeric_limits<dfloat>::max(), platform->o_mempool.slice2); 
  for (int sweep = 0; sweep < 2; sweep++) {
    cds->dirichletBCKernel(mesh->Nelements,
                           cds->fieldOffset[is],
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
                           platform->o_mempool.slice2);

    //take care of Neumann-Dirichlet shared edges across elements
    if(sweep == 0) oogs::startFinish(platform->o_mempool.slice2, 1, cds->fieldOffset[is], ogsDfloat, ogsMax, gsh);
    if(sweep == 1) oogs::startFinish(platform->o_mempool.slice2, 1, cds->fieldOffset[is], ogsDfloat, ogsMin, gsh);
  }
  if (solver->Nmasked) cds->maskCopyKernel(solver->Nmasked, 0, solver->o_maskIds, platform->o_mempool.slice2, platform->o_mempool.slice0);

  //build RHS
  platform->o_mempool.slice1.copyFrom(cds->o_BF, cds->fieldOffset[is] * sizeof(dfloat), 0,  cds->fieldOffsetScan[is] * sizeof(dfloat));
  cds->helmholtzRhsBCKernel(mesh->Nelements,
                            mesh->o_sgeo,
                            mesh->o_vmapM,
                            mesh->o_EToB,
                            is,
                            time,
                            cds->fieldOffset[is],
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            platform->o_mempool.slice0,
                            cds->o_EToB[is],
                            *(cds->o_usrwrk),
                            platform->o_mempool.slice1);

  ellipticSolve(solver, platform->o_mempool.slice1, platform->o_mempool.slice0);

  return platform->o_mempool.slice0;
}



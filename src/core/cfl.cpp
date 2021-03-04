#include "nrs.hpp"
#include "platform.hpp"
#include "linAlg.hpp"

static int firstTime = 1;

void setup(nrs_t* nrs)
{
  mesh_t* mesh = nrs->meshV;
  

  dfloat* dH;
  if(nrs->elementType == QUADRILATERALS || nrs->elementType == HEXAHEDRA) {
    dH = (dfloat*) calloc((mesh->N + 1),sizeof(dfloat));

    for(int n = 0; n < (mesh->N + 1); n++) {
      if(n == 0)
        dH[n] = mesh->gllz[n + 1] - mesh->gllz[n];
      else if(n == mesh->N)
        dH[n] = mesh->gllz[n] - mesh->gllz[n - 1];
      else
        dH[n] = 0.5 * ( mesh->gllz[n + 1] - mesh->gllz[n - 1]);
    }
    for(int n = 0; n < (mesh->N + 1); n++)
      dH[n] = 1.0 / dH[n];

    nrs->o_idH = platform->device.malloc((mesh->N + 1) * sizeof(dfloat), dH);
    free(dH);
  }
  firstTime = 0;
}

dfloat computeCFL(nrs_t* nrs)
{
  mesh_t* mesh = nrs->meshV;
  
  if(firstTime) setup(nrs);

  // Compute cfl factors i.e. dt* U / h
  nrs->cflKernel(mesh->Nelements,
                 nrs->dt[0],
                 mesh->o_vgeo,
                 nrs->o_idH,
                 nrs->meshV->fieldOffset,
                 nrs->o_U,
                 platform->o_slice0);

  return platform->linAlg->max(mesh->Nlocal, platform->o_slice0, platform->comm.mpiComm);
}

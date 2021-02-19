#include "nrs.hpp"
#include "platform.hpp"
#include "linAlg.hpp"

static int firstTime = 1;
static dfloat* tmp;
static occa::memory o_tmp;

void setup(nrs_t* nrs)
{
  mesh_t* mesh = nrs->mesh;
  platform_t* platform = platform_t::getInstance();

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

  tmp = (dfloat*) calloc(nrs->Nblock, sizeof(dfloat));
  o_tmp = platform->device.malloc(nrs->Nblock * sizeof(dfloat), tmp);

  firstTime = 0;
}

dfloat computeCFL(nrs_t* nrs)
{
  mesh_t* mesh = nrs->mesh;
  linAlg_t* linAlg = linAlg_t::getInstance();
  if(firstTime) setup(nrs);

  // Compute cfl factors i.e. dt* U / h
  nrs->cflKernel(mesh->Nelements,
                 nrs->dt[0],
                 mesh->o_vgeo,
                 nrs->o_idH,
                 nrs->fieldOffset,
                 nrs->o_U,
                 nrs->o_wrk0);

  return linAlg->max(mesh->Nlocal, nrs->o_wrk0, mesh->comm);
}

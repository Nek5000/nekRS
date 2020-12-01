#include "nrs.hpp"

static int firstTime = 1;
static dfloat* tmp;
static occa::memory o_tmp;

void setup(nrs_t* nrs)
{
  mesh_t* mesh = nrs->mesh;

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

    nrs->o_idH = mesh->device.malloc((mesh->N + 1) * sizeof(dfloat), dH);
    free(dH);
  }

  tmp = (dfloat*) calloc(nrs->Nblock, sizeof(dfloat));
  o_tmp = mesh->device.malloc(nrs->Nblock * sizeof(dfloat), tmp);

  firstTime = 0;
}

dfloat computeCFL(nrs_t* nrs)
{
  mesh_t* mesh = nrs->mesh;
  if(firstTime) setup(nrs);

  // Compute cfl factors i.e. dt* U / h
  nrs->cflKernel(mesh->Nelements,
                 nrs->dt[0],
                 mesh->o_vgeo,
                 nrs->o_idH,
                 nrs->fieldOffset,
                 nrs->o_U,
                 nrs->o_wrk0);

  // find the local maximum of CFL number
  nrs->maxKernel(nrs->Nlocal, nrs->o_wrk0, o_tmp);
  o_tmp.copyTo(tmp);

  // finish reduction
  dfloat cfl = 0.f;
  for(dlong n = 0; n < nrs->Nblock; ++n)
    cfl  = mymax(cfl, tmp[n]);

  dfloat gcfl = 0.f;
  MPI_Allreduce(&cfl, &gcfl, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);

  return gcfl;
}

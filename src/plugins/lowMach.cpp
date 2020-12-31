#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"

#include "lowMach.hpp"

void lowMach::setup(nrs_t* nrs)
{
  mesh_t* mesh = nrs->mesh;
  int err = 1;
  if(nrs->options.compareArgs("SCALAR00 IS TEMPERATURE", "TRUE")) err = 0;
  if(err) {
    if(mesh->rank == 0) cout << "lowMach requires solving for temperature!\n";
    ABORT(1);
  } 
  udf.div = &lowMach::qtl;
  nrs->options.setArgs("LOWMACH", "TRUE"); 
}

// qtl = 1/(rho*cp*T) * (div[k*grad[T] ] + qvol)
void lowMach::qtl(nrs_t* nrs, dfloat time, occa::memory o_div)
{
  cds_t* cds = nrs->cds;
  mesh_t* mesh = nrs->mesh;

  nrs->gradientVolumeKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    mesh->o_Dmatrices,
    nrs->fieldOffset,
    cds->o_S,
    cds->o_wrk0);

  oogs::startFinish(cds->o_wrk0, nrs->NVfields, nrs->fieldOffset,ogsDfloat, ogsAdd, nrs->gsh);

  const dfloat one = 1.0;
  const dlong Nlocal = mesh->Nelements * mesh->Np;
  for(int field = 0 ; field < nrs->NVfields; ++field){
    occa::memory slice = nrs->o_wrk0 + field * nrs->fieldOffset * sizeof(dfloat);
    nrs->linAlg->axmy(Nlocal, one, nrs->mesh->o_invLMM, slice);
  }

  if(udf.sEqnSource) {
    timer::tic("udfSEqnSource", 1);
    udf.sEqnSource(nrs, time, cds->o_S, cds->o_wrk3);
    timer::toc("udfSEqnSource");
  } else {
    nrs->fillKernel(mesh->Nelements * mesh->Np, 0.0, cds->o_wrk3);
  }

  nrs->qtlKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    mesh->o_Dmatrices,
    nrs->fieldOffset,
    cds->o_wrk0,
    cds->o_S,
    cds->o_diff,
    cds->o_rho,
    cds->o_wrk3,
    o_div);

  oogs::startFinish(o_div, 1, nrs->fieldOffset, ogsDfloat, ogsAdd, nrs->gsh);

  nrs->linAlg->axmy(Nlocal, one, nrs->mesh->o_invLMM, o_div);
}

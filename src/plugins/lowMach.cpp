#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"

#include "lowMach.hpp"

void lowMach::setup(ins_t* ins)
{
  mesh_t* mesh = ins->mesh;
  int err = 0;
  if(ins->Nscalar) {
    if(!ins->cds->compute[0]) err = 1; 
  } else {
    err = 1;
  }
  if(err) {
    if(mesh->rank == 0) cout << "lowMach requires solving for temperature!\n";
    ABORT(1);
  } 
  udf.div = &lowMach::qtl;
  ins->options.setArgs("LOWMACH", "TRUE"); 
}

// qtl = 1/(rho*cp*T) * (div[k*grad[T] ] + qvol)
void lowMach::qtl(ins_t* ins, dfloat time, occa::memory o_div)
{
  cds_t* cds = ins->cds;
  mesh_t* mesh = ins->mesh;

  ins->gradientVolumeKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    mesh->o_Dmatrices,
    ins->fieldOffset,
    cds->o_S,
    cds->o_wrk0);

  oogs::startFinish(cds->o_wrk0, ins->NVfields, ins->fieldOffset,ogsDfloat, ogsAdd, ins->gsh);

  ins->invMassMatrixKernel(
    mesh->Nelements,
    ins->fieldOffset,
    ins->NVfields,
    mesh->o_vgeo,
    ins->mesh->o_invLMM,
    cds->o_wrk0);

  if(udf.sEqnSource) {
    timer::tic("udfSEqnSource", 1);
    udf.sEqnSource(ins, time, cds->o_S, cds->o_wrk3);
    timer::toc("udfSEqnSource");
  } else {
    ins->fillKernel(mesh->Nelements * mesh->Np, 0.0, cds->o_wrk3);
  }

  ins->qtlKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    mesh->o_Dmatrices,
    ins->fieldOffset,
    cds->o_wrk0,
    cds->o_S,
    cds->o_diff,
    cds->o_rho,
    cds->o_wrk3,
    o_div);

  oogs::startFinish(o_div, 1, ins->fieldOffset, ogsDfloat, ogsAdd, ins->gsh);

  ins->invMassMatrixKernel(
    mesh->Nelements,
    ins->fieldOffset,
    1,
    mesh->o_vgeo,
    ins->mesh->o_invLMM,
    o_div);
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"

#include "lowMach.hpp"

static nrs_t* nrs = nullptr;
static linAlg_t* linAlg = nullptr;
void lowMach::setup(nrs_t* _nrs)
{
  nrs = _nrs;
  linAlg = linAlg_t::getInstance();
  mesh_t* mesh = nrs->mesh;
  int err = 1;
  if(nrs->options.compareArgs("SCALAR00 IS TEMPERATURE", "TRUE")) err = 0;
  if(err) {
    if(mesh->rank == 0) cout << "lowMach requires solving for temperature!\n";
    ABORT(1);
  } 
  nrs->options.setArgs("LOWMACH", "TRUE"); 
}


// qtl = 1/(rho*cp*T) * (div[k*grad[T] ] + qvol)
void lowMach::qThermalPerfectGasSingleComponent(nrs_t* nrs, dfloat time, dfloat gamma, occa::memory o_div)
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

  nrs->invMassMatrixKernel(
    mesh->Nelements,
    nrs->fieldOffset,
    nrs->NVfields,
    mesh->o_vgeo,
    nrs->mesh->o_invLMM,
    cds->o_wrk0);

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

  nrs->invMassMatrixKernel(
    mesh->Nelements,
    nrs->fieldOffset,
    1,
    mesh->o_vgeo,
    nrs->mesh->o_invLMM,
    o_div);
  
  if(nrs->pSolver->allNeumann){
    const dfloat dd = (1.0 - gamma) / gamma;
    const dlong Nlocal = nrs->Nlocal;
    occa::memory& o_w1 = nrs->o_wrk0;
    occa::memory& o_w2 = nrs->o_wrk1;
    occa::memory& o_scratch = nrs->o_wrk2;

    // rho * cp = cds->o_rho
    // rho      = nrs->o_rho
    // cp       = cds->o_rho / nrs->o_rho
    nrs->p0thHelperKernel(Nlocal,
      dd,
      cds->o_rho,
      nrs->o_rho,
      nrs->mesh->o_LMM,
      o_w1,
      o_w2
    );

    linAlg_t * linAlg = linAlg_t::getInstance();

    const dfloat p0alpha1 = 1.0 / linAlg->sum(Nlocal, o_w1, mesh->comm);
    linAlg->axmyz(Nlocal, 1.0, mesh->o_LMM, o_div, o_w1);

    const dfloat termQ = linAlg->sum(Nlocal, o_w1, mesh->comm);
    dfloat* scratch = nrs->wrk;
    nrs->surfaceFluxKernel(
      mesh->Nelements,
      mesh->o_sgeo,
      mesh->o_vmapM,
      nrs->o_EToB,
      nrs->fieldOffset,
      nrs->o_Ue,
      o_scratch
    );
    o_scratch.copyTo(scratch, mesh->Nelements * sizeof(dfloat));
    dfloat termV = 0.0;
    for(int i = 0 ; i < mesh->Nelements; ++i) termV += scratch[i];
    MPI_Allreduce(MPI_IN_PLACE, &termV, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

    const dfloat prhs = p0alpha1 * (termQ - termV);
    const dfloat pcoef = nrs->g0 - nrs->dt[0] * prhs;

    dfloat Saqpq = 0.0;
    for(int i = 0 ; i < 3; ++i){
      Saqpq += nrs->extbdfB[i] * nrs->p0th[i];
    }
    const dfloat p0th = Saqpq / pcoef;
    nrs->p0th[2] = nrs->p0th[1];
    nrs->p0th[1] = nrs->p0th[0];
    nrs->p0th[0] = p0th;
    nrs->dp0thdt = prhs * p0th;

    const dfloat weight = -prhs;
    linAlg->axpby(Nlocal, weight, o_w2, 1.0, o_div);
  }
}
void lowMach::dpdt(dfloat gamma, occa::memory o_FU)
{
  linAlg->add(nrs->Nlocal, nrs->dp0thdt * (gamma - 1.0) / gamma, o_FU);
}

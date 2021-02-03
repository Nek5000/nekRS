#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"

#include "lowMach.hpp"

namespace{
static nrs_t* the_nrs = nullptr;
static linAlg_t* the_linAlg = nullptr;
static int qThermal = 0;
static occa::kernel qtlKernel;
static occa::kernel p0thHelperKernel;
void buildKernels(nrs_t* nrs)
{
  mesh_t* mesh = nrs->mesh;
  occa::properties kernelInfo = *(nrs->kernelInfo);
  string fileName;
  int rank = mesh->rank;
  fileName.assign(getenv("NEKRS_INSTALL_DIR"));
  fileName += "/okl/plugins/lowMach.okl";
  for (int r = 0; r < 2; r++) {
    if ((r == 0 && rank == 0) || (r == 1 && rank > 0)) {
      qtlKernel        = mesh->device.buildKernel(fileName.c_str(), "qtlHex3D"  , kernelInfo);
      p0thHelperKernel = mesh->device.buildKernel(fileName.c_str(), "p0thHelper", kernelInfo);
    }
    MPI_Barrier(mesh->comm);
  }
}
}
void lowMach::setup(nrs_t* nrs)
{
  the_nrs = nrs;
  the_linAlg = nrs->linAlg;
  mesh_t* mesh = nrs->mesh;
  int err = 1;
  if(nrs->options.compareArgs("SCALAR00 IS TEMPERATURE", "TRUE")) err = 0;
  if(err) {
    if(mesh->rank == 0) cout << "lowMach requires solving for temperature!\n";
    ABORT(1);
  } 
  buildKernels(nrs);
  nrs->options.setArgs("LOWMACH", "TRUE"); 
}

// qtl = 1/(rho*cp*T) * (div[k*grad[T] ] + qvol)
void lowMach::qThermalPerfectGasSingleComponent(nrs_t* nrs, dfloat time, dfloat gamma, occa::memory o_div)
{
  qThermal = 1;
  cds_t* cds = nrs->cds;
  mesh_t* mesh = nrs->mesh;
  linAlg_t * linAlg = nrs->linAlg;

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

  qtlKernel(
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

    linAlg->axmyz(Nlocal, 1.0, mesh->o_LMM, o_div, nrs->o_wrk0);
    const dfloat termQ = linAlg->sum(Nlocal, nrs->o_wrk0, mesh->comm);

    nrs->surfaceFluxKernel(
      mesh->Nelements,
      mesh->o_sgeo,
      mesh->o_vmapM,
      nrs->o_EToB,
      nrs->fieldOffset,
      nrs->o_Ue,
      nrs->o_wrk0
    );
    nrs->o_wrk0.copyTo(nrs->wrk, mesh->Nelements * sizeof(dfloat));
    dfloat termV = 0.0;
    for(int i = 0 ; i < mesh->Nelements; ++i) termV += nrs->wrk[i];
    MPI_Allreduce(MPI_IN_PLACE, &termV, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

    p0thHelperKernel(Nlocal,
      dd,
      cds->o_rho,
      nrs->o_rho,
      nrs->mesh->o_LMM,
      nrs->o_wrk0,
      nrs->o_wrk1 
    );
    const dfloat prhs = (termQ - termV)/linAlg->sum(Nlocal, nrs->o_wrk0, mesh->comm);;
    linAlg->axpby(Nlocal, -prhs, nrs->o_wrk1, 1.0, o_div);

    dfloat Saqpq = 0.0;
    for(int i = 0 ; i < 3; ++i){
      Saqpq += nrs->extbdfB[i] * nrs->p0th[i];
    }
    nrs->p0th[2] = nrs->p0th[1];
    nrs->p0th[1] = nrs->p0th[0];
    nrs->p0th[0] = Saqpq / (nrs->g0 - nrs->dt[0] * prhs);
    nrs->dp0thdt = prhs * nrs->p0th[0];
  }
  qThermal = 0;
}

void lowMach::dpdt(dfloat gamma, occa::memory o_FU)
{
  if(!qThermal)
    the_linAlg->add(the_nrs->Nlocal, the_nrs->dp0thdt * (gamma - 1.0) / gamma, o_FU);
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"

#include "lowMach.hpp"
#include "linAlg.hpp"

namespace{

static nrs_t* the_nrs = nullptr;
static linAlg_t* the_linAlg = nullptr;
static int qThermal = 0;
static dfloat gamma0 = 1;
static occa::kernel qtlKernel;
static occa::kernel p0thHelperKernel;
static occa::kernel surfaceFluxKernel;

}

void lowMach::buildKernel(nrs_t* nrs)
{
  mesh_t* mesh = nrs->meshV;
  occa::properties kernelInfo = *(nrs->kernelInfo);
  string fileName;
  int rank = platform->comm.mpiRank;
  fileName.assign(getenv("NEKRS_INSTALL_DIR"));
  fileName += "/okl/plugins/lowMach.okl";
  if( BLOCKSIZE < mesh->Nq * mesh->Nq ){
    if(rank == 0)
      printf("ERROR: nrsSurfaceFlux kernel requires BLOCKSIZE >= Nq * Nq."
        "BLOCKSIZE = %d, Nq*Nq = %d\n", BLOCKSIZE, mesh->Nq * mesh->Nq);
    ABORT(EXIT_FAILURE);
  }
  {
    qtlKernel        = platform->device.buildKernel(fileName, "qtlHex3D"  , kernelInfo);
    p0thHelperKernel = platform->device.buildKernel(fileName, "p0thHelper", kernelInfo);
    surfaceFluxKernel = platform->device.buildKernel(fileName, "surfaceFlux", kernelInfo);
  }
}

void lowMach::setup(nrs_t* nrs, dfloat gamma)
{
  the_nrs = nrs;
  gamma0 = gamma;
  the_linAlg = platform->linAlg;
  mesh_t* mesh = nrs->meshV;
  int err = 1;
  if(platform->options.compareArgs("SCALAR00 IS TEMPERATURE", "TRUE")) err = 0;
  if(err) {
    if(platform->comm.mpiRank == 0) cout << "lowMach requires solving for temperature!\n";
    ABORT(1);
  } 
  platform->options.setArgs("LOWMACH", "TRUE"); 
}

// qtl = 1/(rho*cp*T) * (div[k*grad[T] ] + qvol)
void lowMach::qThermalIdealGasSingleComponent(dfloat time, occa::memory o_div)
{
  qThermal = 1;
  nrs_t* nrs = the_nrs;
  cds_t* cds = nrs->cds;
  mesh_t* mesh = nrs->meshV;
  linAlg_t * linAlg = platform->linAlg;

  nrs->gradientVolumeKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    mesh->o_D,
    nrs->fieldOffset,
    cds->o_S,
    platform->o_mempool.slice0);

  oogs::startFinish(platform->o_mempool.slice0, nrs->NVfields, nrs->fieldOffset,ogsDfloat, ogsAdd, nrs->gsh);

  platform->linAlg->axmyVector(
    mesh->Nlocal,
    nrs->fieldOffset,
    0,
    1.0,
    nrs->meshV->o_invLMM,
    platform->o_mempool.slice0);

  if(udf.sEqnSource) {
    platform->timer.tic("udfSEqnSource", 1);
    udf.sEqnSource(nrs, time, cds->o_S, platform->o_mempool.slice3);
    platform->timer.toc("udfSEqnSource");
  } else {
    platform->linAlg->fill(mesh->Nelements * mesh->Np, 0.0, platform->o_mempool.slice3);
  }

  qtlKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    mesh->o_D,
    nrs->fieldOffset,
    platform->o_mempool.slice0,
    cds->o_S,
    cds->o_diff,
    cds->o_rho,
    platform->o_mempool.slice3,
    o_div);

  oogs::startFinish(o_div, 1, nrs->fieldOffset, ogsDfloat, ogsAdd, nrs->gsh);

  platform->linAlg->axmy(
    mesh->Nlocal,
    1.0,
    nrs->meshV->o_invLMM,
    o_div);
  
  if(nrs->pSolver->allNeumann){
    const dfloat dd = (1.0 - gamma0) / gamma0;
    const dlong Nlocal = mesh->Nlocal;

    linAlg->axmyz(Nlocal, 1.0, mesh->o_LMM, o_div, platform->o_mempool.slice0);
    const dfloat termQ = linAlg->sum(Nlocal, platform->o_mempool.slice0, platform->comm.mpiComm);

    surfaceFluxKernel(
      mesh->Nelements,
      mesh->o_sgeo,
      mesh->o_vmapM,
      nrs->o_EToB,
      nrs->fieldOffset,
      nrs->o_Ue,
      platform->o_mempool.slice0
    );
    platform->o_mempool.slice0.copyTo(platform->mempool.slice0, mesh->Nelements * sizeof(dfloat));
    dfloat termV = 0.0;
    for(int i = 0 ; i < mesh->Nelements; ++i) termV += platform->mempool.slice0[i];
    MPI_Allreduce(MPI_IN_PLACE, &termV, 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);

    p0thHelperKernel(Nlocal,
      dd,
      cds->o_rho,
      nrs->o_rho,
      nrs->meshV->o_LMM,
      platform->o_mempool.slice0,
      platform->o_mempool.slice1 
    );
    const dfloat prhs = (termQ - termV)/linAlg->sum(Nlocal, platform->o_mempool.slice0, platform->comm.mpiComm);
    linAlg->axpby(Nlocal, -prhs, platform->o_mempool.slice1, 1.0, o_div);

    dfloat Saqpq = 0.0;
    for(int i = 0 ; i < nrs->nBDF; ++i){
      Saqpq += nrs->coeffBDF[i] * nrs->p0th[i];
    }
    nrs->p0th[2] = nrs->p0th[1];
    nrs->p0th[1] = nrs->p0th[0];

    nrs->p0th[0] = Saqpq / (nrs->g0 - nrs->dt[0] * prhs);
    nrs->dp0thdt = prhs * nrs->p0th[0];
  }
  qThermal = 0;
}

void lowMach::dpdt(occa::memory o_FU)
{
  nrs_t* nrs = the_nrs;
  mesh_t* mesh = nrs->meshV;
  if(!qThermal)
    platform->linAlg->add(mesh->Nlocal, nrs->dp0thdt * (gamma0 - 1.0) / gamma0, o_FU);
}

#include "cds.hpp"
#include "avm.hpp"
#include "udf.hpp"
#include <limits>
#include <string>
#include <array>

/**
 * Persson's artificial viscosity method (http://persson.berkeley.edu/pub/persson06shock.pdf) with P1
 **/

namespace avm {


static occa::kernel relativeMassHighestModeKernel;
static occa::kernel computeMaxViscKernel;
static occa::kernel interpolateP1Kernel;

static occa::memory o_vertexIds;
static occa::memory o_r;
static occa::memory o_s;
static occa::memory o_t;
static occa::memory o_diffOld; // diffusion from initial state
static double cachedDt = -1.0;
static bool recomputeUrst = false;
namespace
{

void allocateMemory(cds_t* cds)
{
  mesh_t * mesh = cds->mesh[0];
  if(udf.properties == nullptr) o_diffOld = platform->device.malloc(cds->fieldOffsetSum, sizeof(dfloat), cds->o_diff);
  o_vertexIds = platform->device.malloc(mesh->Nverts * sizeof(int), mesh->vertexNodes);
  o_r = platform->device.malloc(mesh->Np * sizeof(dfloat), mesh->r);
  o_s = platform->device.malloc(mesh->Np * sizeof(dfloat), mesh->s);
  o_t = platform->device.malloc(mesh->Np * sizeof(dfloat), mesh->t);
}

void compileKernels(cds_t* cds)
{
  std::string kernelName;

  kernelName = "relativeMassHighestMode";
  relativeMassHighestModeKernel =
    platform->kernels.get(kernelName);

  kernelName = "computeMaxVisc";
  computeMaxViscKernel =
    platform->kernels.get(kernelName);

  kernelName = "interpolateP1";
  interpolateP1Kernel =
    platform->kernels.get(kernelName);
}

}
void setup(cds_t* cds)
{
  allocateMemory(cds);
  compileKernels(cds);
}

occa::memory computeEps(nrs_t* nrs, const dfloat time, const dlong scalarIndex, occa::memory o_S)
{
  cds_t* cds = nrs->cds;

  const dfloat TOL = 1e-10;
  if(std::abs(cachedDt - time) > TOL)
    recomputeUrst = true;
  
  cachedDt = time;
  
  mesh_t* mesh = cds->mesh[scalarIndex];
  int Nblock = (cds->mesh[scalarIndex]->Nlocal+BLOCKSIZE-1)/BLOCKSIZE;

  occa::memory& o_logRelativeMassHighestMode = platform->o_mempool.slice0;
  occa::memory& o_filteredField = platform->o_mempool.slice1;
  occa::memory& o_hpfResidual = platform->o_mempool.slice2;
  occa::memory& o_aliasedUrst = platform->o_mempool.slice6;

  // artificial viscosity magnitude
  occa::memory o_epsilon = platform->o_mempool.slice5;
  platform->linAlg->fill(cds->fieldOffset[scalarIndex], 0.0, o_epsilon);

  const dfloat p = mesh->N;
  
  relativeMassHighestModeKernel(
    mesh->Nelements,
    scalarIndex,
    cds->fieldOffsetScan[scalarIndex],
    cds->o_filterMT,
    mesh->o_LMM,
    o_S,
    o_filteredField,
    o_logRelativeMassHighestMode
  );

  const bool useHPFResidual = cds->options[scalarIndex].compareArgs("REGULARIZATION METHOD", "HPF_RESIDUAL");

  dfloat Uinf = 1.0;
  if(useHPFResidual){

    occa::memory o_rhoField = cds->o_rho + cds->fieldOffsetScan[scalarIndex] * sizeof(dfloat);
    
    if(recomputeUrst){
      nrs->UrstKernel(
        cds->meshV->Nelements,
        cds->meshV->o_vgeo,
        nrs->fieldOffset,
        nrs->o_U,
        cds->meshV->o_U,
        o_aliasedUrst
      );
      recomputeUrst = false;
    }

    cds->strongAdvectionVolumeKernel(cds->meshV->Nelements,
                                     cds->meshV->o_vgeo,
                                     mesh->o_D,
                                     cds->vFieldOffset,
                                     0,
                                     o_filteredField,
                                     o_aliasedUrst,
                                     o_rhoField,
                                     o_hpfResidual);

    occa::memory o_S_field = o_S + cds->fieldOffsetScan[scalarIndex] * sizeof(dfloat);
    
    const dfloat Uavg = platform->linAlg->weightedNorm2(
      mesh->Nlocal,
      mesh->o_LMM,
      o_S_field,
      platform->comm.mpiComm
    ) / sqrt(mesh->volume);

    platform->linAlg->fill(mesh->Nlocal,
      Uavg,
      o_filteredField);

    platform->linAlg->axpby(
      mesh->Nlocal,
      -1.0,
      o_S_field,
      1.0,
      o_filteredField
    );

    if(Uavg > 0.0){
      Uinf = platform->linAlg->max(mesh->Nlocal, o_filteredField, platform->comm.mpiComm);
    }
    Uinf = 1.0 / Uinf;
  }

  const dfloat logReferenceSensor = -4.0 * log10(p);

  dfloat coeff = 0.5;
  cds->options[scalarIndex].getArgs("REGULARIZATION VISMAX COEFF", coeff);

  dfloat rampParameter = 1.0;
  cds->options[scalarIndex].getArgs("REGULARIZATION RAMP CONSTANT", rampParameter);


  dfloat scalingCoeff = 1.0;
  cds->options[scalarIndex].getArgs("REGULARIZATION SCALING COEFF", scalingCoeff);

  computeMaxViscKernel(
    mesh->Nelements,
    cds->vFieldOffset,
    logReferenceSensor,
    rampParameter,
    coeff,
    scalingCoeff,
    Uinf,
    useHPFResidual,
    mesh->o_x,
    mesh->o_y,
    mesh->o_z,
    cds->o_U,
    o_hpfResidual,
    o_logRelativeMassHighestMode,
    o_epsilon // max(|df/du|) <- max visc
  );

  const bool makeCont = cds->options[scalarIndex].compareArgs("REGULARIZATION AVM C0", "TRUE");
  if(makeCont){
    oogs_t* gsh;
    if(scalarIndex){
      gsh = cds->gsh;
    } else {
      gsh = cds->gshT;
    }
    oogs::startFinish(o_epsilon, 1, cds->fieldOffset[scalarIndex], ogsDfloat, ogsMax, gsh);
    interpolateP1Kernel(
      mesh->Nelements,
      o_vertexIds,
      o_r,
      o_s,
      o_t,
      o_epsilon
    );
  }

  return o_epsilon;

}

void apply(nrs_t* nrs, const dfloat time, const dlong scalarIndex, occa::memory o_S)
{
  cds_t* cds = nrs->cds;
  const int verbose = platform->options.compareArgs("VERBOSE", "TRUE");
  mesh_t* mesh = cds->mesh[scalarIndex];
  const dlong scalarOffset = cds->fieldOffsetScan[scalarIndex];
  if(udf.properties == nullptr){ // restore viscosity from previous state
    if(verbose && platform->comm.mpiRank == 0) printf("Resetting properties...\n");
    cds->o_diff.copyFrom(o_diffOld,
      cds->fieldOffset[scalarIndex] * sizeof(dfloat),
      cds->fieldOffsetScan[scalarIndex] * sizeof(dfloat),
      cds->fieldOffsetScan[scalarIndex] * sizeof(dfloat)
    );
  }
  occa::memory o_eps = computeEps(nrs, time, scalarIndex, o_S);
  if(verbose && platform->comm.mpiRank == 0){
    const dfloat maxEps = platform->linAlg->max(
      mesh->Nlocal,
      o_eps,
      platform->comm.mpiComm
    );
    const dfloat minEps = platform->linAlg->min(
      mesh->Nlocal,
      o_eps,
      platform->comm.mpiComm
    );
    occa::memory o_S_slice = cds->o_diff + cds->fieldOffsetScan[scalarIndex] * sizeof(dfloat);
    const dfloat maxDiff = platform->linAlg->max(
      mesh->Nlocal,
      o_S_slice,
      platform->comm.mpiComm
    );
    const dfloat minDiff = platform->linAlg->min(
      mesh->Nlocal,
      o_S_slice,
      platform->comm.mpiComm
    );
    printf("Applying a min/max artificial viscosity of (%f,%f) to a field with min/max visc (%f,%f) (field = %d)\n",
      minEps, maxEps, minDiff, maxDiff, scalarIndex);
  }

  platform->linAlg->axpby(
    mesh->Nlocal,
    1.0,
    o_eps,
    1.0,
    cds->o_diff,
    0,
    scalarOffset
  );

  if(verbose && platform->comm.mpiRank == 0){
    occa::memory o_S_slice = cds->o_diff + cds->fieldOffsetScan[scalarIndex] * sizeof(dfloat);
    const dfloat maxDiff = platform->linAlg->max(
      mesh->Nlocal,
      o_S_slice,
      platform->comm.mpiComm
    );
    const dfloat minDiff = platform->linAlg->min(
      mesh->Nlocal,
      o_S_slice,
      platform->comm.mpiComm
    );
    printf("Field (%d) now has a min/max visc: (%f,%f)\n", scalarIndex, minDiff, maxDiff);
  }

}

} // namespace
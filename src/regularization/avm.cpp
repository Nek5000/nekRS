#include "cds.hpp"
#include "avm.hpp"
#include "udf.hpp"
#include <limits>
#include <string>

/**
 * Artificial viscosity method (https://arxiv.org/pdf/1810.02152.pdf), [0,1] viscosity weighting currently deemed as not useful
 **/

namespace avm {

static occa::kernel filterScalarNormKernel;
static occa::kernel applyAVMKernel;
static occa::kernel computeMaxViscKernel;

static occa::memory o_artVisc;
static occa::memory o_diffOld; // diffusion from initial state

void allocateMemory(cds_t* cds)
{
  if(udf.properties == nullptr) o_diffOld = platform->device.malloc(cds->fieldOffsetSum, sizeof(dfloat), cds->o_diff);

  o_artVisc = platform->device.malloc(cds->NSfields * cds->mesh[0]->Np, sizeof(dfloat));
  platform->linAlg->fill(cds->NSfields * cds->mesh[0]->Np, 1.0, o_artVisc);
}

void compileKernels(cds_t* cds)
{
  std::string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string oklpath = install_dir + "/okl/cds/regularization/";
  std::string filename = oklpath + "filterScalarNorm.okl";
  occa::properties info = platform->kernelInfo;
  info["defines/" "p_Nq"] = cds->mesh[0]->Nq;
  info["defines/" "p_Np"] = cds->mesh[0]->Np;
  filterScalarNormKernel =
    platform->device.buildKernel(filename,
                             "filterScalarNorm",
                             info);

  filename = oklpath + "applyAVM.okl";
  applyAVMKernel =
    platform->device.buildKernel(filename,
                             "applyAVM",
                             info);

  filename = oklpath + "computeMaxVisc.okl";
  computeMaxViscKernel =
    platform->device.buildKernel(filename,
                             "computeMaxVisc",
                             info);
}

void setup(cds_t* cds)
{
  allocateMemory(cds);
  compileKernels(cds);

  const dfloat alpha = -std::log(std::numeric_limits<dfloat>::epsilon());
  auto superGaussian = [alpha](const dfloat x, const dfloat lambda)
  {
    return std::exp(-alpha * std::pow(x, 2.0 * lambda));
  };

  dfloat* artVisc = (dfloat*) calloc(cds->mesh[0]->Np, sizeof(dfloat));
  for(dlong is = 0 ; is < cds->NSfields; ++is){
    if(cds->options[is].compareArgs("AVM C0", "TRUE")){
      dfloat lambda = 1.0;
      cds->options[is].getArgs("AVM LAMBDA", lambda);
      for(dlong point = 0; point < cds->mesh[is]->Np; ++point){
        const dfloat r = cds->mesh[is]->r[point];
        const dfloat s = cds->mesh[is]->s[point];
        const dfloat t = cds->mesh[is]->t[point];
        dfloat visc = superGaussian(r, lambda);
        visc *= superGaussian(s, lambda);
        visc *= superGaussian(t, lambda);
        visc = cbrt(visc);
        artVisc[point] = visc;
      }
      o_artVisc.copyFrom(artVisc, cds->mesh[is]->Np * sizeof(dfloat), is * cds->mesh[is]->Np * sizeof(dfloat));
    }
  }

  free(artVisc);

}

occa::memory computeEps(cds_t* cds, const dfloat time, const dlong scalarIndex, occa::memory o_S)
{
  dfloat sensorSensitivity = 1.0;
  cds->options[scalarIndex].getArgs("SENSITIVITY", sensorSensitivity);
  mesh_t* mesh = cds->mesh[scalarIndex];
  int Nblock = (cds->mesh[scalarIndex]->Nlocal+BLOCKSIZE-1)/BLOCKSIZE;

  occa::memory& o_logShockSensor = platform->o_mempool.slice0;

  // artificial viscosity magnitude
  occa::memory o_epsilon = platform->o_mempool.slice2;
  platform->linAlg->fill(cds->fieldOffset[scalarIndex], 0.0, o_epsilon);

  const dfloat p = mesh->N;

  dfloat sensorOrder = 4.0;
  cds->options[scalarIndex].getArgs("SENSOR ORDER", sensorOrder);
  
  filterScalarNormKernel(
    mesh->Nelements,
    scalarIndex,
    cds->fieldOffsetScan[scalarIndex],
    cds->o_filterMT,
    sensorSensitivity,
    sensorOrder,
    p,
    o_S,
    o_logShockSensor
  );

  // see footnote after eq 63
  const dfloat logReferenceSensor = -2.0;

  dfloat coeff = 0.5;
  cds->options[scalarIndex].getArgs("VISMAX COEFF", coeff);

  dfloat rampParameter = 1.0;
  cds->options[scalarIndex].getArgs("RAMP CONSTANT", rampParameter);

  computeMaxViscKernel(
    mesh->Nelements,
    cds->vFieldOffset,
    logReferenceSensor,
    rampParameter,
    coeff,
    mesh->o_x,
    mesh->o_y,
    mesh->o_z,
    cds->o_U,
    o_logShockSensor,
    o_epsilon // max(|df/du|) <- max visc
  );

  return o_epsilon;

}

void apply(cds_t* cds, const dfloat time, const dlong scalarIndex, occa::memory o_S)
{
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
  occa::memory o_eps = computeEps(cds, time, scalarIndex, o_S);
  if(verbose && platform->comm.mpiRank == 0){
    const dfloat maxEps = platform->linAlg->max(
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
    printf("Applying a max artificial viscosity of: %f to a field with max visc: %f\n", maxEps, maxDiff);
  }

  applyAVMKernel(
    mesh->Nelements,
    scalarOffset,
    scalarIndex,
    o_artVisc,
    o_eps,
    cds->o_diff
  );
  if(verbose && platform->comm.mpiRank == 0){
    occa::memory o_S_slice = cds->o_diff + cds->fieldOffsetScan[scalarIndex] * sizeof(dfloat);
    const dfloat maxDiff = platform->linAlg->max(
      mesh->Nlocal,
      o_S_slice,
      platform->comm.mpiComm
    );
    printf("Field now has a max visc: %f\n", maxDiff);
  }

}

} // namespace
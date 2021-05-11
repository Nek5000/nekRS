#include "cds.hpp"
#include "avm.hpp"
#include "udf.hpp"
#include <limits>
#include <string>

/**
 * Artificial viscosity method (https://arxiv.org/pdf/1810.02152.pdf)
 **/


namespace avm {


static occa::kernel filterScalarNormKernel;
static occa::kernel applyAVMKernel;
static occa::kernel computeMaxViscKernel;
static occa::kernel computeElementLengthScaleKernel;

static occa::memory o_artVisc;
static occa::memory o_diffOld; // diffusion from initial state
namespace
{

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
  filename = oklpath + "computeElementLengthScale.okl";
  computeElementLengthScaleKernel =
    platform->device.buildKernel(filename,
                             "computeElementLengthScale",
                             info);
}

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
  int useErrorIndicator = 1;
  useErrorIndicator = cds->options[scalarIndex].compareArgs("AVM ERROR INDICATOR", "TRUE");
  dfloat sensorSensitivity = 1.0;
  cds->options[scalarIndex].getArgs("SENSITIVITY", sensorSensitivity);
  mesh_t* mesh = cds->mesh[scalarIndex];
  int Nblock = (cds->mesh[scalarIndex]->Nlocal+BLOCKSIZE-1)/BLOCKSIZE;

  occa::memory& o_logShockSensor = platform->o_mempool.slice0;
  occa::memory& o_elemLengths = platform->o_mempool.slice1;

  occa::memory& o_filteredField = platform->o_mempool.slice3;
  occa::memory& o_advectedFilteredField = platform->o_mempool.slice4;
  occa::memory& o_ones = platform->o_mempool.slice5;
  occa::memory& o_wrk = platform->o_mempool.slice6;
  platform->linAlg->fill(cds->fieldOffset[scalarIndex], 1.0, o_ones);

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
    o_logShockSensor,
    o_filteredField
  );

  dfloat Sinf = 1.0;

  if(useErrorIndicator){
    // convect filtered field
    const dlong cubatureOffset = std::max(cds->vFieldOffset, cds->meshV->Nelements * cds->meshV->cubNp);
    if(cds->options[scalarIndex].compareArgs("ADVECTION TYPE", "CUBATURE"))
      cds->advectionStrongCubatureVolumeKernel(
        cds->meshV->Nelements,
        mesh->o_vgeo,
        mesh->o_cubDiffInterpT,
        mesh->o_cubInterpT,
        mesh->o_cubProjectT,
        cds->vFieldOffset,
        0,
        cubatureOffset,
        o_filteredField,
        cds->o_Urst,
        o_ones,
        o_advectedFilteredField);
    else
      cds->advectionStrongVolumeKernel(
        cds->meshV->Nelements,
        mesh->o_D,
        cds->vFieldOffset,
        0,
        o_filteredField,
        cds->o_Urst,
        o_ones,
        o_advectedFilteredField);
  
    occa::memory o_S_slice = o_S + cds->fieldOffsetScan[scalarIndex] * sizeof(dfloat);
    const dfloat Savg = 
      platform->linAlg->weightedNorm2(
        mesh->Nlocal,
        mesh->o_LMM,
        o_S_slice,
        platform->comm.mpiComm
      ) / sqrt(mesh->volume);

    if(Savg > 0.0){
      platform->linAlg->fill(cds->fieldOffset[scalarIndex], Savg, o_wrk);
      platform->linAlg->axpby(
        mesh->Nlocal,
        -1.0,
        o_S_slice,
        1.0,
        o_wrk
      );
      platform->linAlg->abs(mesh->Nlocal, o_wrk);
      Sinf = platform->linAlg->max(mesh->Nlocal, o_wrk, platform->comm.mpiComm);
    }
  }

  Sinf = 1.0 / Sinf;

  computeElementLengthScaleKernel(
    mesh->Nelements,
    mesh->o_x,
    mesh->o_y,
    mesh->o_z,
    o_elemLengths
  );

  // see footnote after eq 63
  const dfloat logReferenceSensor = -2.0;

  dfloat coeff = 0.5;
  cds->options[scalarIndex].getArgs("VISMAX COEFF", coeff);

  dfloat errCoeff = 1.0;
  cds->options[scalarIndex].getArgs("ERROR COEFF", errCoeff);

  dfloat rampParameter = 1.0;
  cds->options[scalarIndex].getArgs("RAMP CONSTANT", rampParameter);

  computeMaxViscKernel(
    mesh->Nelements,
    cds->vFieldOffset,
    logReferenceSensor,
    rampParameter,
    coeff,
    errCoeff,
    Sinf,
    useErrorIndicator,
    mesh->o_x,
    mesh->o_y,
    mesh->o_z,
    cds->o_U,
    o_logShockSensor,
    o_elemLengths,
    o_advectedFilteredField,
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
    const dfloat minDiff = platform->linAlg->min(
      mesh->Nlocal,
      o_S_slice,
      platform->comm.mpiComm
    );
    printf("Field (%d) now has a min/max visc: (%f,%f)\n", scalarIndex, minDiff, maxDiff);
  }

}

} // namespace
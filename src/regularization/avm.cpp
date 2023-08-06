#include <limits>
#include <string>
#include <array>

#include "nrssys.hpp"
#include "cds.hpp"
#include "avm.hpp"
#include "udf.hpp"
#include "hpf.hpp"

/**
 * Persson's artificial viscosity method (http://persson.berkeley.edu/pub/persson06shock.pdf) with P1
 **/

namespace avm
{

static occa::kernel relativeMassHighestModeKernel;
static occa::kernel computeMaxViscKernel;
static occa::kernel interpolateP1Kernel;

static occa::memory o_vertexIds;
static occa::memory o_r;
static occa::memory o_s;
static occa::memory o_t;
static std::vector<occa::memory> o_diff0;
static std::vector<occa::memory> o_filterMT;

static double cachedDt = -1.0;
static bool recomputeUrst = false;

namespace
{

}

void setup(cds_t *cds)
{
  mesh_t *mesh = cds->mesh[0];

  for (int is = 0; is < cds->NSfields; is++) {
    std::string sid = scalarDigitStr(is);

    if (platform->options.compareArgs("SCALAR" + sid + " REGULARIZATION METHOD", "AVM_RESIDUAL") ||
        platform->options.compareArgs("SCALAR" + sid + " REGULARIZATION METHOD", "AVM_HIGHEST_MODAL_DECAY")) {

      nrsCheck(cds->mesh[is]->N < 5,
               platform->comm.mpiComm,
               EXIT_FAILURE,
               "%s\n",
               "AVM requires polynomialOrder >= 5!");

      if (udf.properties == nullptr) {
        o_diff0.push_back(platform->device.malloc(cds->fieldOffset[is], sizeof(dfloat)));
        o_diff0[is].copyFrom(cds->o_diff,
                             cds->fieldOffset[is] * sizeof(dfloat),
                             0,
                             cds->fieldOffsetScan[is] * sizeof(dfloat));
      }

      int filterNc = -1;
      platform->options.getArgs("SCALAR" + sid + " REGULARIZATION HPF MODES", filterNc);
      o_filterMT.push_back(hpfSetup(cds->mesh[is], filterNc));
    } else {
      o_filterMT.push_back(o_NULL);
      if (udf.properties == nullptr) {
        o_diff0.push_back(o_NULL);
      }
    }
  }

  o_vertexIds = platform->device.malloc(mesh->Nverts * sizeof(int), mesh->vertexNodes);
  o_r = platform->device.malloc(mesh->Np * sizeof(dfloat), mesh->r);
  o_s = platform->device.malloc(mesh->Np * sizeof(dfloat), mesh->s);
  o_t = platform->device.malloc(mesh->Np * sizeof(dfloat), mesh->t);

  std::string kernelName;

  kernelName = "relativeMassHighestMode";
  relativeMassHighestModeKernel = platform->kernels.get(kernelName);

  kernelName = "computeMaxVisc";
  computeMaxViscKernel = platform->kernels.get(kernelName);

  kernelName = "interpolateP1";
  interpolateP1Kernel = platform->kernels.get(kernelName);
}

occa::memory computeEps(nrs_t *nrs, const dfloat time, const dlong scalarIndex, occa::memory o_S)
{
  cds_t *cds = nrs->cds;

  const dfloat TOL = 1e-10;
  if (std::abs(cachedDt - time) > TOL) {
    recomputeUrst = true;
  }

  cachedDt = time;

  mesh_t *mesh = cds->mesh[scalarIndex];

  occa::memory &o_logRelativeMassHighestMode = platform->o_mempool.slice0;
  occa::memory &o_filteredField = platform->o_mempool.slice1;
  occa::memory &o_hpfResidual = platform->o_mempool.slice2;
  occa::memory &o_epsilon = platform->o_mempool.slice5;
  occa::memory &o_aliasedUrst = platform->o_mempool.slice6;

  // artificial viscosity magnitude
  platform->linAlg->fill(cds->fieldOffset[scalarIndex], 0.0, o_epsilon);

  relativeMassHighestModeKernel(mesh->Nelements,
                                scalarIndex,
                                cds->fieldOffsetScan[scalarIndex],
                                o_filterMT[scalarIndex],
                                mesh->o_LMM,
                                o_S,
                                o_filteredField,
                                o_logRelativeMassHighestMode);

  std::string sid = scalarDigitStr(scalarIndex);

  const int useHPFResidual =
      platform->options.compareArgs("SCALAR" + sid + " REGULARIZATION METHOD", "AVM_RESIDUAL");

  dfloat Uinf = 1.0;
  if (useHPFResidual) {

    occa::memory o_rhoField = cds->o_rho + cds->fieldOffsetScan[scalarIndex] * sizeof(dfloat);

    if (recomputeUrst) {
      nrs->UrstKernel(cds->meshV->Nelements,
                      cds->meshV->o_vgeo,
                      nrs->fieldOffset,
                      nrs->o_U,
                      cds->meshV->o_U,
                      o_aliasedUrst);
      recomputeUrst = false;
    }

    cds->strongAdvectionVolumeKernel(cds->meshV->Nelements,
                                     1,
                                     cds->meshV->o_vgeo,
                                     mesh->o_D,
                                     cds->o_compute + scalarIndex * sizeof(dlong),
                                     cds->o_fieldOffsetScan + scalarIndex * sizeof(dlong),
                                     cds->vFieldOffset,
                                     o_filteredField,
                                     o_aliasedUrst,
                                     o_rhoField,
                                     o_hpfResidual);

    occa::memory o_S_field = o_S + cds->fieldOffsetScan[scalarIndex] * sizeof(dfloat);

    const dfloat Uavg =
        platform->linAlg->weightedNorm2(mesh->Nlocal, mesh->o_LMM, o_S_field, platform->comm.mpiComm) /
        sqrt(mesh->volume);

    platform->linAlg->fill(mesh->Nlocal, Uavg, o_filteredField);

    platform->linAlg->axpby(mesh->Nlocal, -1.0, o_S_field, 1.0, o_filteredField);

    if (Uavg > 0.0) {
      Uinf = platform->linAlg->max(mesh->Nlocal, o_filteredField, platform->comm.mpiComm);
    }
    Uinf = 1.0 / Uinf;
  }

  dfloat coeff = 0.5;
  platform->options.getArgs("SCALAR" + sid + " REGULARIZATION VISMAX COEFF", coeff);

  dfloat rampParameter = 1.0;
  platform->options.getArgs("SCALAR" + sid + " REGULARIZATION MDH ACTIVATION WIDTH", rampParameter);

  dfloat threshold = -4.0;
  platform->options.getArgs("SCALAR" + sid + " REGULARIZATION MDH THRESHOLD", threshold);

  dfloat scalingCoeff = 1.0;
  platform->options.getArgs("SCALAR" + sid + " REGULARIZATION SCALING COEFF", scalingCoeff);

  const dfloat logReferenceSensor = threshold * log10(mesh->N);

  computeMaxViscKernel(mesh->Nelements,
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

  const bool makeCont = platform->options.compareArgs("SCALAR" + sid + " REGULARIZATION AVM C0", "TRUE");
  if (makeCont) {
    oogs_t *gsh;
    if (scalarIndex) {
      gsh = cds->gsh;
    } else {
      gsh = cds->gshT;
    }
    oogs::startFinish(o_epsilon, 1, cds->fieldOffset[scalarIndex], ogsDfloat, ogsMax, gsh);
    interpolateP1Kernel(mesh->Nelements, o_vertexIds, o_r, o_s, o_t, o_epsilon);
  }

  return o_epsilon;
}

void apply(nrs_t *nrs, const dfloat time, const dlong scalarIndex, occa::memory o_S)
{
  cds_t *cds = nrs->cds;
  const int verbose = platform->options.compareArgs("VERBOSE", "TRUE");
  mesh_t *mesh = cds->mesh[scalarIndex];

  // restore inital viscosity
  if (udf.properties == nullptr) {
    cds->o_diff.copyFrom(o_diff0[scalarIndex],
                         cds->fieldOffset[scalarIndex] * sizeof(dfloat),
                         cds->fieldOffsetScan[scalarIndex] * sizeof(dfloat));
  }

  occa::memory o_eps = computeEps(nrs, time, scalarIndex, o_S);

  if (verbose) {
    const dfloat maxEps = platform->linAlg->max(mesh->Nlocal, o_eps, platform->comm.mpiComm);
    const dfloat minEps = platform->linAlg->min(mesh->Nlocal, o_eps, platform->comm.mpiComm);
    occa::memory o_S_slice = cds->o_diff + cds->fieldOffsetScan[scalarIndex] * sizeof(dfloat);
    const dfloat maxDiff = platform->linAlg->max(mesh->Nlocal, o_S_slice, platform->comm.mpiComm);
    const dfloat minDiff = platform->linAlg->min(mesh->Nlocal, o_S_slice, platform->comm.mpiComm);

    if (platform->comm.mpiRank == 0) {
      printf("Applying a min/max artificial viscosity of (%f,%f) to a field with min/max visc (%f,%f) (field "
             "= %d)\n",
             minEps,
             maxEps,
             minDiff,
             maxDiff,
             scalarIndex);
    }
  }

  platform->linAlg->axpby(mesh->Nlocal, 1.0, o_eps, 1.0, cds->o_diff, 0, cds->fieldOffsetScan[scalarIndex]);

  if (verbose) {
    occa::memory o_S_slice = cds->o_diff + cds->fieldOffsetScan[scalarIndex] * sizeof(dfloat);
    const dfloat maxDiff = platform->linAlg->max(mesh->Nlocal, o_S_slice, platform->comm.mpiComm);
    const dfloat minDiff = platform->linAlg->min(mesh->Nlocal, o_S_slice, platform->comm.mpiComm);

    if (platform->comm.mpiRank == 0) {
      printf("Field (%d) now has a min/max visc: (%f,%f)\n", scalarIndex, minDiff, maxDiff);
    }
  }
}

} // namespace avm

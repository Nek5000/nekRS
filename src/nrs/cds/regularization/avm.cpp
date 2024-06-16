#include <array>

#include "nekrsSys.hpp"
#include "cds.hpp"
#include "avm.hpp"
#include "udf.hpp"
#include "lowPassFilter.hpp"

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
static std::vector<occa::memory> o_filterRT;

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

      nekrsCheck(cds->mesh[is]->N < 5,
                 platform->comm.mpiComm,
                 EXIT_FAILURE,
                 "%s\n",
                 "AVM requires polynomialOrder >= 5!");

      if (cds->userProperties == nullptr) {
        o_diff0.push_back(platform->device.malloc<dfloat>(cds->fieldOffset[is]));
        o_diff0[is].copyFrom(cds->o_diff, cds->fieldOffset[is], 0, cds->fieldOffsetScan[is]);
      }

      int filterNc = -1;
      platform->options.getArgs("SCALAR" + sid + " REGULARIZATION HPF MODES", filterNc);
      o_filterRT.push_back(lowPassFilterSetup(cds->mesh[is], filterNc));
    } else {
      o_filterRT.push_back(o_NULL);
      if (cds->userProperties == nullptr) {
        o_diff0.push_back(o_NULL);
      }
    }
  }

  o_vertexIds = platform->device.malloc<int>(mesh->Nverts, mesh->vertexNodes);
  o_r = platform->device.malloc<dfloat>(mesh->Np, mesh->r);
  o_s = platform->device.malloc<dfloat>(mesh->Np, mesh->s);
  o_t = platform->device.malloc<dfloat>(mesh->Np, mesh->t);

  std::string kernelName;

  kernelName = "relativeMassHighestMode";
  relativeMassHighestModeKernel = platform->kernelRequests.load(kernelName);

  kernelName = "computeMaxVisc";
  computeMaxViscKernel = platform->kernelRequests.load(kernelName);

  kernelName = "interpolateP1";
  interpolateP1Kernel = platform->kernelRequests.load(kernelName);
}

occa::memory computeEps(cds_t *cds, const double time, const dlong scalarIndex, occa::memory o_S)
{
  if (std::abs(time - cachedDt) / time > 10 * std::numeric_limits<dfloat>::epsilon()) {
    recomputeUrst = true;
  }

  cachedDt = time;

  mesh_t *mesh = cds->mesh[scalarIndex];

  occa::memory o_logRelativeMassHighestMode =
      platform->o_memPool.reserve<dfloat>(cds->fieldOffset[scalarIndex]);
  occa::memory o_filteredField = platform->o_memPool.reserve<dfloat>(cds->fieldOffset[scalarIndex]);
  occa::memory o_hpfResidual = platform->o_memPool.reserve<dfloat>(cds->fieldOffset[scalarIndex]);
  occa::memory o_epsilon = platform->o_memPool.reserve<dfloat>(cds->fieldOffset[scalarIndex]);

  // artificial viscosity magnitude
  platform->linAlg->fill(cds->fieldOffset[scalarIndex], 0.0, o_epsilon);

  relativeMassHighestModeKernel(mesh->Nelements,
                                scalarIndex,
                                cds->fieldOffsetScan[scalarIndex],
                                o_filterRT[scalarIndex],
                                mesh->o_LMM,
                                o_S,
                                o_filteredField,
                                o_logRelativeMassHighestMode);

  std::string sid = scalarDigitStr(scalarIndex);

  dfloat Uinf = 1.0;

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
                       0 /* useHPFResidual */,
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

void apply(cds_t *cds, const double time, const dlong scalarIndex, occa::memory o_S)
{
  const int verbose = platform->options.compareArgs("VERBOSE", "TRUE");
  mesh_t *mesh = cds->mesh[scalarIndex];

  // restore inital viscosity
  if (cds->userProperties == nullptr) {
    cds->o_diff.copyFrom(o_diff0[scalarIndex],
                         cds->fieldOffset[scalarIndex],
                         cds->fieldOffsetScan[scalarIndex]);
  }

  occa::memory o_eps = computeEps(cds, time, scalarIndex, o_S);

  if (verbose) {
    const dfloat maxEps = platform->linAlg->max(mesh->Nlocal, o_eps, platform->comm.mpiComm);
    const dfloat minEps = platform->linAlg->min(mesh->Nlocal, o_eps, platform->comm.mpiComm);
    occa::memory o_S_slice = cds->o_diff + cds->fieldOffsetScan[scalarIndex];
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
    occa::memory o_S_slice = cds->o_diff + cds->fieldOffsetScan[scalarIndex];
    const dfloat maxDiff = platform->linAlg->max(mesh->Nlocal, o_S_slice, platform->comm.mpiComm);
    const dfloat minDiff = platform->linAlg->min(mesh->Nlocal, o_S_slice, platform->comm.mpiComm);

    if (platform->comm.mpiRank == 0) {
      printf("Field (%d) now has a min/max visc: (%f,%f)\n", scalarIndex, minDiff, maxDiff);
    }
  }
}

} // namespace avm

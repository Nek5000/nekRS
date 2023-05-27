#if !defined(nekrs_nekrs_hpp_)
#define nekrs_nekrs_hpp_

#include "nrssys.hpp"
#include "mesh3D.h"
#include "elliptic.h"
#include "cds.hpp"
#include "linAlg.hpp"
#include "timer.hpp"
#include "platform.hpp"
#include "neknek.hpp"
#include "cvode.hpp"
#include "fldFile.hpp"

std::vector<std::string> fieldsToSolve(setupAide& options);

struct nrs_t {

  static constexpr double targetTimeBenchmark {0.2};

  bool multiSession;

  int dim, elementType;

  mesh_t *_mesh = nullptr;
  mesh_t *meshV = nullptr;

  elliptic_t *uSolver = nullptr;
  elliptic_t *vSolver = nullptr;
  elliptic_t *wSolver = nullptr;
  elliptic_t *uvwSolver = nullptr;
  elliptic_t *pSolver = nullptr;
  elliptic_t *meshSolver = nullptr;

  cds_t *cds = nullptr;

  neknek_t *neknek = nullptr;
  cvode_t *cvode = nullptr;

  oogs_t *gsh = nullptr;
  oogs_t *gshMesh = nullptr;

  dlong ellipticWrkOffset;

  int flow;
  int cht;
  int Nscalar;
  int NVfields;

  dlong fieldOffset;
  dlong cubatureOffset;

  int timeStepConverged;

  dfloat dt[3], idt;
  dfloat g0, ig0;
  dfloat CFL, unitTimeCFL;

  dfloat timePrevious;

  dfloat p0th[3] = {0.0, 0.0, 0.0};
  dfloat p0the = 0.0;
  dfloat dp0thdt;

  dfloat alpha0Ref;

  int nEXT;
  int nBDF;

  int tstep;
  int lastStep;
  int isOutputStep;
  int outputForceStep;

  int nRK, Nsubsteps;
  dfloat *coeffsfRK, *weightsRK, *nodesRK;
  occa::memory o_coeffsfRK, o_weightsRK;

  dfloat *U, *P;
  occa::memory o_U, o_P;

  occa::memory o_Ue;

  occa::memory o_div;

  occa::memory o_rho, o_mue;
  occa::memory o_meshRho, o_meshMue;

  dfloat *usrwrk;
  occa::memory o_usrwrk;

  occa::memory o_idH;

  occa::memory o_BF;
  occa::memory o_FU;

  occa::memory o_prop, o_ellipticCoeff;

  dfloat *coeffEXT, *coeffBDF;
  occa::memory o_coeffEXT, o_coeffBDF;

  int *EToB;
  int *EToBMeshVelocity;
  occa::memory o_EToB;
  occa::memory o_EToBMeshVelocity;

  occa::memory o_EToBVVelocity;
  occa::memory o_EToBVMeshVelocity;

  occa::memory o_Uc, o_Pc;
  occa::memory o_prevProp;

  occa::memory o_relUrst;
  occa::memory o_Urst;

  occa::properties *kernelInfo;

  int filterNc;
  dfloat filterS;
  occa::memory o_filterMT;

  occa::kernel filterRTKernel;
  occa::kernel advectMeshVelocityKernel;
  occa::kernel pressureAddQtlKernel;
  occa::kernel pressureStressKernel;
  occa::kernel extrapolateKernel;
  occa::kernel subCycleRKKernel;
  occa::kernel subCycleInitU0Kernel;
  occa::kernel nStagesSum3Kernel;
  occa::kernel wgradientVolumeKernel;

  occa::kernel subCycleVolumeKernel, subCycleCubatureVolumeKernel;
  occa::kernel subCycleSurfaceKernel, subCycleCubatureSurfaceKernel;
  occa::kernel subCycleRKUpdateKernel;
  occa::kernel subCycleStrongCubatureVolumeKernel;
  occa::kernel subCycleStrongVolumeKernel;

  occa::kernel computeFaceCentroidKernel;
  occa::kernel computeFieldDotNormalKernel;

  occa::kernel UrstCubatureKernel;
  occa::kernel UrstKernel;

  occa::kernel advectionVolumeKernel;
  occa::kernel advectionCubatureVolumeKernel;

  occa::kernel strongAdvectionVolumeKernel;
  occa::kernel strongAdvectionCubatureVolumeKernel;

  occa::kernel gradientVolumeKernel;

  occa::kernel wDivergenceVolumeKernel;
  occa::kernel divergenceVolumeKernel;
  occa::kernel divergenceSurfaceKernel;

  occa::kernel divergenceStrongVolumeKernel;
  occa::kernel sumMakefKernel;
  occa::kernel pressureRhsKernel;
  occa::kernel pressureDirichletBCKernel;

  occa::kernel velocityRhsKernel;
  occa::kernel velocityNeumannBCKernel;
  occa::kernel velocityDirichletBCKernel;

  occa::kernel cflKernel;

  occa::kernel setEllipticCoeffKernel;
  occa::kernel setEllipticCoeffPressureKernel;

  occa::kernel curlKernel;

  occa::kernel SijOijKernel;

  occa::kernel maskCopyKernel;
  occa::kernel maskCopy2Kernel;
  occa::kernel maskKernel;

  occa::memory o_zeroNormalMaskVelocity;
  occa::memory o_zeroNormalMaskMeshVelocity;
  occa::kernel averageNormalBcTypeKernel;
  occa::kernel fixZeroNormalMaskKernel;
  occa::kernel initializeZeroNormalMaskKernel;

  occa::kernel applyZeroNormalMaskKernel;
};

int nrsFinalize(nrs_t *nrs);

void evaluateProperties(nrs_t *nrs, const double timeNew);

void compileKernels();

int numberActiveFields(nrs_t *nrs);

#endif

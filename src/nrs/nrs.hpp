#if !defined(nrs_nekrs_hpp_)
#define nrs_nekrs_hpp_

#include "platform.hpp"
#include "linAlg.hpp"
#include "elliptic.hpp"
#include "cds.hpp"
#include "neknek.hpp"
#include "randomVector.hpp"
#include "aeroForce.hpp"
#include "iofldFactory.hpp"

class nrs_t : public solver_t
{

  using userVelocitySource_t = std::function<void(double)>;
  using userScalarSource_t = std::function<void(double)>;
  using userProperties_t = std::function<void(double)>;
  using userDivergence_t = std::function<void(double)>;
  using preFluid_t = std::function<void(double, int)>;
  using postScalar_t = std::function<void(double, int)>;
  using userConvergenceCheck_t = std::function<bool(int)>;
  using userVelocityImplicitLinearTerm_t = std::function<occa::memory(double)>;
  using userScalarImplicitLinearTerm_t = std::function<occa::memory(double, int)>;

public:
  userVelocitySource_t userVelocitySource = nullptr;
  userScalarSource_t userScalarSource = nullptr;
  userProperties_t userProperties = nullptr;
  userDivergence_t userDivergence = nullptr;
  preFluid_t preFluid = nullptr;
  postScalar_t postScalar = nullptr;
  userConvergenceCheck_t userConvergenceCheck = nullptr;
  userVelocityImplicitLinearTerm_t userVelocityImplicitLinearTerm = nullptr;
  userScalarImplicitLinearTerm_t userScalarImplicitLinearTerm = nullptr;

  void addUserCheckpointField(const std::string& name, const std::vector<deviceMemory<dfloat>>& o_fld)
  {
    std::vector<occa::memory> o_fld_;
    for (const auto& entry : o_fld) o_fld_.push_back(entry);

    userCheckpointFields.push_back( {name, o_fld_} );
  };

  bool multiSession = false;

  int elementType = HEXAHEDRA;

  mesh_t *mesh = nullptr;
  mesh_t *meshV = nullptr; 

  elliptic *uSolver = nullptr;
  elliptic *vSolver = nullptr;
  elliptic *wSolver = nullptr;
  elliptic *uvwSolver = nullptr;
  elliptic *pSolver = nullptr;
  elliptic *meshSolver = nullptr;

  cds_t *cds = nullptr;

  neknek_t *neknek = nullptr;

  oogs_t *gsh = nullptr;

  QQt *qqt = nullptr;
  QQt *qqtT = nullptr;

  oogs_t *gshMesh = nullptr;

  int flow = 1;
  int cht = 0;
  int Nscalar = 0;
  int NVfields = 3;

  dlong fieldOffset;
  dlong cubatureOffset;

  std::unique_ptr<iofld> checkpointWriter;

  int timeStepConverged = 1;

  dfloat dt[3] = {0.0, 0.0, 0.0};
  dfloat g0 = 0;

  double timePrevious;

  dfloat p0th[3] = {0.0, 0.0, 0.0};
  dfloat dp0thdt = 0;

  dfloat alpha0Ref = 1;

  int nEXT = 3;
  int nBDF = 3;

  int tstep = 0;
  int lastStep = 0;
  int outerCorrector = 1;
  int checkpointStep = 0;
  int outputForceStep = 0;

  int Nsubsteps = 0;

  occa::memory o_U;
  occa::memory o_Ue;

  occa::memory o_P;
  occa::memory o_div;

  occa::memory o_prop;
  occa::memory o_rho, o_mue;
  occa::memory o_meshRho, o_meshMue;

  deviceMemory<dfloat> o_usrwrk;

  occa::memory o_idH;

  occa::memory o_JwF;
  occa::memory o_NLT;

  dfloat *coeffEXT, *coeffBDF;
  occa::memory o_coeffEXT, o_coeffBDF;

  occa::memory o_EToB;
  occa::memory o_EToBMeshVelocity;

  occa::memory o_EToBVVelocity;
  occa::memory o_EToBVMeshVelocity;

  occa::memory o_Uc, o_Pc;
  occa::memory o_prevProp;

  occa::memory o_relUrst;
  occa::memory o_Urst;

  dfloat filterS = 0;
  occa::memory o_filterRT;

  occa::kernel filterRTKernel;
  occa::kernel advectMeshVelocityKernel;
  occa::kernel pressureAddQtlKernel;
  occa::kernel pressureStressKernel;
  occa::kernel extrapolateKernel;

  occa::kernel subCycleRKKernel;
  occa::kernel subCycleInitU0Kernel;
  occa::kernel nStagesSum3Kernel;
  occa::kernel wgradientVolumeKernel;

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

  nrs_t();

  void init();

  std::string id() const
  {
    return "nrs";
  };

  void initStep(double time, dfloat dt, int tstep);
  dfloat adjustDt(int tstep);
  bool runStep(std::function<bool(int)> convergenceCheck, int stage);
  void finishStep();

  void saveSolutionState();
  void restoreSolutionState();

  void printMinMax();
  void printRunStat(int step);
  void printStepInfo(double time, int tstep, bool printStepInfo, bool printVerboseInfo);

  void makeNLT(double time, int tstep, occa::memory &o_Usubcycling);
  dfloat computeCFL();
  dfloat computeCFL(dfloat dt);

  void flowRatePrintInfo(bool verboseInfo);
  bool adjustFlowRate(int tstep, double time);
  dfloat flowRatescaleFactor();

  void evaluateProperties(const double time);
  void evaluateDivergence(const double time);

  occa::memory advectionSubcycling(int nEXT, double time);

  int numberActiveFields();

  AeroForce *aeroForces(int nbID, const occa::memory &o_bID, const occa::memory &o_Sij_ = o_NULL);

  // output in row-major order 
  occa::memory strainRotationRate(const occa::memory &o_U, bool smooth = true);
  occa::memory strainRotationRate(bool smooth = true);
  occa::memory strainRate(const occa::memory &o_U, bool smooth = true);
  occa::memory strainRate(bool smooth = true);

  void Qcriterion(occa::memory &o_Q);
  void Qcriterion(const occa::memory &o_U, occa::memory &o_Q);
  occa::memory Qcriterion(const occa::memory &o_U);
  occa::memory Qcriterion();

  void restartFromFile(const std::string& restartStr);
  void writeCheckpoint(double t, int step, bool enforceOutXYZ = false, bool enforceFP64 = false, int Nout = 0, bool uniform = false);

  void finalize();
  int setLastStep(double timeNew, int tstep, double elapsedTime);
  int lastStepLocalSession(double timeNew, int tstep, double elapsedTime);

  void copyToNek(double time, int tstep, bool updateMesh = false);
  void copyToNek(double time, bool updateMesh = false);

  void copyFromNek(double &time);
  void copyFromNek();
  void getICFromNek();

private:
  void initInnerStep(double time, dfloat dt, int tstep);
  bool runInnerStep(std::function<bool(int)> convergenceCheck, int stage, bool outerConverged);
  void finishInnerStep();

  void initOuterStep(double time, dfloat dt, int tstep);
  void runOuterStep(std::function<bool(int)> convergenceCheck, int stage);
  void finishOuterStep();

  int tStepOuterStart;
  double timeOuterStart;

  occa::memory o_Usave;
  occa::memory o_Psave;
  occa::memory o_NLTsave;
  occa::memory o_Urstsave;
  occa::memory o_Upropsave;

  occa::memory o_LMMsave;
  occa::memory o_Umeshsave;
  occa::memory o_xsave;
  occa::memory o_ysave;
  occa::memory o_zsave;

  void setIC();

  std::vector<std::pair<std::string, std::vector<occa::memory>>> userCheckpointFields;

};

void nrsSetDefaultSettings(setupAide *options);

static std::vector<std::string> nrsFieldsToSolve(setupAide &options)
{
  int Nscalar = 0;
  options.getArgs("NUMBER OF SCALARS", Nscalar);

  std::vector<std::string> fields;

  if (!options.compareArgs("MESH SOLVER", "NONE")) {
    fields.push_back("mesh");
  }

  if (!options.compareArgs("VELOCITY SOLVER", "NONE")) {
    fields.push_back("velocity");
  }

  for (int i = 0; i < Nscalar; i++) {
    const auto sid = scalarDigitStr(i);
    if (!options.compareArgs("SCALAR" + sid + " SOLVER", "NONE")) {
      fields.push_back("scalar" + sid);
    }
  }
  return fields;
}

#endif

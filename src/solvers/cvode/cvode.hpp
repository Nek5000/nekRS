#ifndef CVODE_SOLVER_HPP_
#define CVODE_SOLVER_HPP_

#include <limits>
#include "nrssys.hpp"
#include "occa.hpp"
#include <functional>
#include <map>
#include <array>
#include <vector>
#include <tuple>
#include <memory>

#ifdef ENABLE_CVODE
#include <cvode/cvode.h>
#endif

class nrs_t;

class cvode_t {
public:
#ifdef ENABLE_CVODE
  static constexpr bool enabled = true;
#else
  static constexpr bool enabled = false;
#endif
  using userRHS_t = std::function<
      void(nrs_t *nrs, dfloat time, dfloat t0, occa::memory o_y, occa::memory o_ydot)>;
  using userJacobian_t = std::function<
      void(nrs_t *nrs, dfloat time, dfloat t0, occa::memory o_y, occa::memory o_ydot)>;
  using userLocalPointSource_t =
      std::function<void(nrs_t *nrs, dlong LFieldOffset, occa::memory o_y, occa::memory o_ydot)>;
  using userPostNrsToCv_t = std::function<void(nrs_t *nrs, occa::memory o_LField)>;
  using userPostCvToNrs_t = std::function<void(nrs_t *nrs, occa::memory o_EField)>;

  cvode_t(nrs_t *nrs);
  ~cvode_t();

  void initialize();
  void solve(dfloat t0, dfloat t1, int tstep);

  void setRHS(userRHS_t _userRHS) { userRHS = _userRHS; }
  void setJacobian(userJacobian_t _userJacobian) { userJacobian = _userJacobian; }
  void setLocalPointSource(userLocalPointSource_t _userLocalPointSource);

  void setUserPostCvToNrs(userPostCvToNrs_t _userPostCvToNrs) { userPostCvToNrs = _userPostCvToNrs; }
  void setUserPostNrsToCv(userPostCvToNrs_t _userPostNrsToCv) { userPostNrsToCv = _userPostNrsToCv; }

  void printInfo(bool printVerboseInfo) const;

  bool isRhsEvaluation() const { return rhsEval; }
  void setIsRhsEvaluation(bool _rhsEval) { rhsEval = _rhsEval; }

  bool areDetailedTimersEnabled() const { return detailedTimersEnabled; }
  void enableDetailedTimers() { detailedTimersEnabled = true; }

  void setIsJacobianEvaluation(bool _jacEval) { jacEval = _jacEval; }
  bool isJacobianEvaluation() const { return jacEval; }

  int timeStep() const { return externalTStep; }
  void setTimeStep(int tstep) { externalTStep = tstep; }
  double time() const { return tnekRS; }

  void rhs(dfloat time, occa::memory o_y, occa::memory o_ydot);
  void jtvRHS(dfloat time, occa::memory o_y, occa::memory o_ydot);
  dlong numEquations() const { return nEq; }

  // getters needed for CVLsJacTimesVecFn
  void *getCvodeMem() { return cvodeMem; }
  double sigmaScale() const { return sigScale; }

  bool mixedPrecisionJtv() const { return mixedPrecisionJtvEnabled; }

  // returns array in E-vector layout that maps E-vector points
  // to the corresponding L-vector point (if unique),
  // else contains -1
  occa::memory &Lpoints() { return o_EToLUnique; }

  dfloat &g0() { return _g0; }
  dfloat &dt() { return dtCvode[0]; }
  dfloat *coeffBDF() { return _coeffBDF.data(); }
  dfloat *coeffEXT() { return _coeffEXT.data(); }

  // compute error weights for CVODE
  void computeErrorWeight(occa::memory o_y, occa::memory o_ewt);

  // CVODE solver statistics
  long numSteps() const;
  long numRHSEvals() const;
  long numNonlinSolveIters() const;
  long numLinIters() const;

  // hand nekRS-style E-vector to CVODE, which uses an L-vector
  void nrsToCv(occa::memory o_EFeild, occa::memory o_LField);

  // unpack CVODE L-vector into nekRS-style E-vector
  void cvToNrs(occa::memory o_LField, occa::memory o_EField);

  void defaultRHS(dfloat time, dfloat t0, occa::memory o_y, occa::memory o_ydot);

  void printTimers();
  void resetTimers();

  std::string scope() const { return timerScope; }
  void setTimerScope(std::string scope) { timerScope = scope; }

  occa::memory o_pointSource; // scratch field for point source
  occa::memory o_vgeoPfloat;

private:

  nrs_t* nrs;

  std::string timerName = "cvode_t::";
  std::string timerScope;
  std::string rhsTagName() const;

  // package data to pass in as user data to cvode
  struct userData_t {

    userData_t(platform_t *_platform, nrs_t *_nrs, cvode_t *_cvode)
        : platform(_platform), nrs(_nrs), cvode(_cvode)
    {
    }

    platform_t *platform;
    nrs_t *nrs;
    cvode_t *cvode;
  };
  std::shared_ptr<userData_t> userdata;

  dlong LFieldOffset;

  // most recent time from nekRS -- used to compute dt in CVODE integration call
  mutable double tnekRS;
  mutable int externalTStep;
  mutable bool localPointSourceEval = false;
  mutable bool jacEval = false;
  mutable bool rhsEval = false;

  bool detailedTimersEnabled = false;

  bool mixedPrecisionJtvEnabled = false;

  bool verboseCVODE = false;

  bool sharedRho = false;

  int minCvodeScalarId;
  int maxCvodeScalarId;

  mutable long int prevNsteps = 0;
  mutable long int prevNrhs = 0;
  mutable long int prevNli = 0;
  mutable long int prevNni = 0;

  bool recycleProperties;

  bool isInitialized = false;

  dfloat tprev = std::numeric_limits<dfloat>::max();
  occa::memory o_U;     // CVODE is responsible for correctly handling the fluid velocity state
  occa::memory o_meshU; // CVODE is responsible for correctly handling the mesh velocity state
  occa::memory o_xyz0;

  void setupEToLMapping();
  void setupDirichletMask();
  void applyDirichlet(dfloat time);

  userRHS_t userRHS;
  userJacobian_t userJacobian;

  userLocalPointSource_t userLocalPointSource;
  userPostCvToNrs_t userPostCvToNrs;
  userPostNrsToCv_t userPostNrsToCv;

  void makeq(dfloat time);

  static constexpr int maxTimestepperOrder = 3;
  std::array<dfloat, maxTimestepperOrder> _coeffBDF;
  std::array<dfloat, maxTimestepperOrder> _coeffEXT;
  std::array<dfloat, maxTimestepperOrder> dtCvode;

  dfloat relTol;
  occa::memory o_absTol; // absolute tolerances for CVODE, per scalar

  dfloat _g0;

  dlong Nscalar;

  occa::memory o_coeffExt;
  occa::memory o_EToLUnique;
  occa::memory o_EToL;

  occa::memory o_cvodeScalarIds;
  occa::memory o_scalarIds;
  // host shadows
  std::vector<dlong> scalarIds;
  std::vector<dlong> cvodeScalarIds;

  // time-lagged Dirichlet values
  occa::memory o_maskValues;

  // Dirichlet values, extrapolated to the current time
  occa::memory o_maskIds;
  dlong Nmasked;
  dlong maskOffset;           // page-aligned offset for indexing into o_maskValues

  occa::kernel weakLaplacianKernel;
  occa::kernel nrsToCvKernel;
  occa::kernel cvToNrsKernel;
  occa::kernel mapToMaskedPointKernel;
  occa::kernel extrapolateDirichletKernel;
  occa::kernel errorWeightKernel;
  occa::kernel fusedAddRhoDivKernel;

  dlong nEq;

  long long int nEqTotal;

  dfloat sigScale = 1.0;

  // cvode internals
  void *cvodeMem;
#ifdef ENABLE_CVODE
  N_Vector y;
  N_Vector cvodeY;
#endif
  occa::memory o_cvodeY;
  occa::memory o_invDegree; // in L-vector format
};

#endif

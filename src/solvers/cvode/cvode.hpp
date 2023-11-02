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
#include "LVector.hpp"

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
      void(nrs_t *nrs, double time, double t0, const  LVector_t<dfloat> & o_y,  LVector_t<dfloat> & o_ydot)>;
  using userJacobian_t = std::function<
      void(nrs_t *nrs, double time, double t0, const  LVector_t<dfloat> & o_y,  LVector_t<dfloat> & o_ydot)>;
  using userLocalPointSourceL_t =
      std::function<void(nrs_t *nrs, double time, const  LVector_t<dfloat> & o_y,  LVector_t<dfloat> & o_ydot)>;
  using userLocalPointSourceE_t =
      std::function<void(nrs_t *nrs, double time, const occa::memory& o_y, occa::memory &_ydot)>;

  using userPostNrsToCv_t = std::function<void(nrs_t *nrs,  LVector_t<dfloat> & o_LField, bool isYdot)>;
  using userPostCvToNrs_t = std::function<void(nrs_t *nrs, occa::memory o_EField, bool isYdot)>;
  using userMakeq_t = std::function<void(nrs_t *nrs, double time, occa::memory& o_FS)>;
  using userPreSolve_t = std::function<void(nrs_t *nrs)>;
  using userPostSolve_t = std::function<void(nrs_t *nrs)>;

  cvode_t(nrs_t *nrs);
  ~cvode_t();

  void initialize();
  void solve(double t0, double t1, int tstep);

  void setRHS(userRHS_t _userRHS) { userRHS = _userRHS; }
  void setJacobian(userJacobian_t _userJacobian) { userJacobian = _userJacobian; }
  void setLocalPointSource(userLocalPointSourceL_t _userLocalPointSource);
  void setLocalPointSource(userLocalPointSourceE_t _userLocalPointSource);

  void setUserPostCvToNrs(userPostCvToNrs_t _userPostCvToNrs) { userPostCvToNrs = _userPostCvToNrs; }
  void setUserPostNrsToCv(userPostNrsToCv_t _userPostNrsToCv) { userPostNrsToCv = _userPostNrsToCv; }
  void setUserMakeq(userMakeq_t _userMakeq) { userMakeq = _userMakeq; }
  void setUserPreSolve(userPreSolve_t _userPreSolve) { userPreSolve = _userPreSolve; }

  void setUserPostSolve(userPostSolve_t _userPostSolve) { userPostSolve = _userPostSolve; }

  void printInfo(bool printVerboseInfo);

  bool isRhsEvaluation() const { return rhsEval; }
  void setIsRhsEvaluation(bool _rhsEval) { rhsEval = _rhsEval; }

  bool areDetailedTimersEnabled() const { return detailedTimersEnabled; }
  void enableDetailedTimers() { detailedTimersEnabled = true; }

  void setIsJacobianEvaluation(bool _jacEval) { jacEval = _jacEval; }
  bool isJacobianEvaluation() const { return jacEval; }

  int timeStep() const { return externalTStep; }
  void setTimeStep(int tstep) { externalTStep = tstep; }
  double time() const { return tExternal; }

  void rhs(double time, const  LVector_t<dfloat> & o_y,  LVector_t<dfloat> & o_ydot);
  void jtvRhs(double time, const  LVector_t<dfloat> & o_y,  LVector_t<dfloat> & o_ydot);
  dlong numEquations() const { return nEq; }
  int numScalars() const { return Nscalar; } 

  bool mixedPrecisionJtv() const { return mixedPrecisionJtvEnabled; }

  dfloat &g0() { return _g0; }
  dfloat &dt() { return dtCvode[0]; }
  dfloat *coeffBDF() { return _coeffBDF.data(); }
  dfloat *coeffEXT() { return _coeffEXT.data(); }

  void updateCounters();
  void resetCounters();



  long numSteps() const;
  long numRHSEvals() const;
  long numNonlinSolveIters() const;
  long numLinIters() const;

  void printTimers();
  void resetTimers();

  std::string scope() const { return timerScope; }
  void setTimerScope(std::string scope) { timerScope = scope; }

  occa::memory o_pointSource; // scratch field for point source
  occa::memory o_vgeoPfloat;

private:
  long int nsteps;
  long int nrhs;
  long int nni;
  long int nli;

  std::shared_ptr<LVector_t<dfloat>> YLVec;
  std::shared_ptr<LVector_t<dfloat>> YdotLVec;

  oogs_t *gsh;

#ifdef ENABLE_CVODE
  // CVODE function pointers (required to access private data members)
  int cvodeRHS(double time, N_Vector Y, N_Vector Ydot);
  int cvodeJtv(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, N_Vector work);
  int cvodeErrorWt(N_Vector y, N_Vector ewt);
#endif

  int jtvRhs(double time, const occa::memory& o_y, occa::memory& o_ydot);
  int jtv(double t, const occa::memory& o_v, const occa::memory& o_y, const occa::memory& o_fy, 
          occa::memory& o_work, occa::memory& o_Jv);
 
  void defaultRHS(double time, double t0, const  LVector_t<dfloat> & o_y,  LVector_t<dfloat> & o_ydot);
  
  // hand nekRS-style E-vector to CVODE, which uses an L-vector
  void nrsToCv(occa::memory o_EFeild,  LVector_t<dfloat> & o_LField, bool isYdot);

  // unpack CVODE L-vector into nekRS-style E-vector
  void cvToNrs(const  LVector_t<dfloat> & o_LField, occa::memory o_EField, bool isYdot);

  nrs_t* nrs;

  std::string timerName = "cvode_t::";
  std::string timerScope;
  std::string rhsTagName() const;

  std::string linearSolverType;

  // most recent time from nekRS -- used to compute dt in CVODE integration call
  mutable double tExternal;
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

  double tprev = std::numeric_limits<double>::max();
  occa::memory o_U0;
  occa::memory o_meshU0;
  occa::memory o_xyz0;

  void setupDirichletMask();
  void applyDirichlet(double time);

  userRHS_t userRHS;
  userJacobian_t userJacobian;

  userLocalPointSourceL_t userLocalPointSourceL;
  userLocalPointSourceE_t userLocalPointSourceE;

  userPostCvToNrs_t userPostCvToNrs;
  userPostNrsToCv_t userPostNrsToCv;
  userMakeq_t userMakeq;
  userPreSolve_t userPreSolve;
  userPostSolve_t userPostSolve;

  void makeq(double time, occa::memory& o_FS);

  static constexpr int maxTimestepperOrder = 3;
  std::array<dfloat, maxTimestepperOrder> _coeffBDF;
  std::array<dfloat, maxTimestepperOrder> _coeffEXT;
  std::array<dfloat, maxTimestepperOrder> dtCvode;

  dfloat relTol;
  occa::memory o_absTol; // absolute tolerances for CVODE, per scalar

  dfloat _g0;

  int Nscalar;

  occa::memory o_rhoCpAvg;

  occa::memory o_coeffExt;

  occa::memory o_cvodeScalarIds;
  occa::memory o_scalarIds;
  // host shadows
  std::vector<dlong> scalarIds;
  std::vector<dlong> cvodeScalarIds;

  occa::memory o_ewt;

  // time-lagged Dirichlet values
  occa::memory o_maskValues;

  // Dirichlet values, extrapolated to the current time
  occa::memory o_maskIds;
  dlong Nmasked;
  dlong maskOffset;           // page-aligned offset for indexing into o_maskValues

  occa::kernel weakLaplacianKernel;
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

};

#endif

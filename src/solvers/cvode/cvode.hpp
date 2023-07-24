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
  
  using LVec = LVector_t<dfloat>;

#ifdef ENABLE_CVODE
  static constexpr bool enabled = true;
#else
  static constexpr bool enabled = false;
#endif
  using userRHS_t = std::function<
      void(nrs_t *nrs, double time, double t0, const  LVector_t<dfloat> & o_y,  LVector_t<dfloat> & o_ydot)>;
  using userJacobian_t = std::function<
      void(nrs_t *nrs, double time, double t0, const  LVector_t<dfloat> & o_y,  LVector_t<dfloat> & o_ydot)>;
  using userLocalPointSource_t =
      std::function<void(nrs_t *nrs, const  LVector_t<dfloat> & o_y,  LVector_t<dfloat> & o_ydot)>;
  using userPostNrsToCv_t = std::function<void(nrs_t *nrs,  LVector_t<dfloat> & o_LField, bool isYdot)>;
  using userPostCvToNrs_t = std::function<void(nrs_t *nrs, occa::memory o_EField, bool isYdot)>;

  cvode_t(nrs_t *nrs);
  ~cvode_t();

  void initialize();
  void solve(double t0, double t1, int tstep);

  void setRHS(userRHS_t _userRHS) { userRHS = _userRHS; }
  void setJacobian(userJacobian_t _userJacobian) { userJacobian = _userJacobian; }
  void setLocalPointSource(userLocalPointSource_t _userLocalPointSource);

  void setUserPostCvToNrs(userPostCvToNrs_t _userPostCvToNrs) { userPostCvToNrs = _userPostCvToNrs; }
  void setUserPostNrsToCv(userPostNrsToCv_t _userPostNrsToCv) { userPostNrsToCv = _userPostNrsToCv; }

  void printInfo(bool printVerboseInfo) const;

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
  void jtvRHS(double time, const  LVector_t<dfloat> & o_y,  LVector_t<dfloat> & o_ydot);
  dlong numEquations() const { return nEq; }

  bool mixedPrecisionJtv() const { return mixedPrecisionJtvEnabled; }

  dfloat &g0() { return _g0; }
  dfloat &dt() { return dtCvode[0]; }
  dfloat *coeffBDF() { return _coeffBDF.data(); }
  dfloat *coeffEXT() { return _coeffEXT.data(); }

  // CVODE solver statistics
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

  // L-vector specific accessors
  const auto & meshes() const { return YLVec->meshes(); }
  const auto & offsets() const { return YLVec->offsets(); }
  auto getLocalPointSource() { return userLocalPointSource; }

private:
  
  std::shared_ptr<LVec> YLVec;
  std::shared_ptr<LVec> YdotLVec;

#ifdef ENABLE_CVODE
  // CVODE function pointers (required to access private data members)
  int cvodeRHS(double time, N_Vector Y, N_Vector Ydot);
  int cvodeJtvRHS(double time, N_Vector Y, N_Vector Ydot);
  int cvodeJtv(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, N_Vector work);
  int cvodeErrorWt(N_Vector y, N_Vector ewt);
#endif
  
  void defaultRHS(double time, double t0, const  LVector_t<dfloat> & o_y,  LVector_t<dfloat> & o_ydot);
  
  // hand nekRS-style E-vector to CVODE, which uses an L-vector
  void nrsToCv(occa::memory o_EFeild,  LVector_t<dfloat> & o_LField, bool isYdot);

  // unpack CVODE L-vector into nekRS-style E-vector
  void cvToNrs(const  LVector_t<dfloat> & o_LField, occa::memory o_EField, bool isYdot);

  nrs_t* nrs;

  std::string timerName = "cvode_t::";
  std::string timerScope;
  std::string rhsTagName() const;

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
  occa::memory o_U;     // CVODE is responsible for correctly handling the fluid velocity state
  occa::memory o_meshU; // CVODE is responsible for correctly handling the mesh velocity state
  occa::memory o_xyz0;

  void setupDirichletMask();
  void applyDirichlet(double time);

  userRHS_t userRHS;
  userJacobian_t userJacobian;

  userLocalPointSource_t userLocalPointSource;
  userPostCvToNrs_t userPostCvToNrs;
  userPostNrsToCv_t userPostNrsToCv;

  void makeq(double time);

  static constexpr int maxTimestepperOrder = 3;
  std::array<dfloat, maxTimestepperOrder> _coeffBDF;
  std::array<dfloat, maxTimestepperOrder> _coeffEXT;
  std::array<dfloat, maxTimestepperOrder> dtCvode;

  dfloat relTol;
  occa::memory o_absTol; // absolute tolerances for CVODE, per scalar

  dfloat _g0;

  dlong Nscalar;

  occa::memory o_rhoCpAvg;

  occa::memory o_coeffExt;

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

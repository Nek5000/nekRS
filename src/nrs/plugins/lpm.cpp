#include <inttypes.h>

#include "platform.hpp"
#include "nekInterfaceAdapter.hpp" // for nek::coeffAB
#include "lpm.hpp"
#include "neknek.hpp"
#include "nrs.hpp"
#include "pointInterpolation.hpp"
#include "tuple_for_each.hpp"

#include "gslib.h" // needed for sarray_transfer

lpm_t::lpm_t(dfloat bb_tol_, dfloat newton_tol_)
    : nrs(dynamic_cast<nrs_t *>(platform->solver)), 
      solverOrder(nrs->nEXT), 
      bb_tol(bb_tol_),
      newton_tol(newton_tol_), 
      interp(std::make_unique<pointInterpolation_t>(nrs->mesh, platform->comm.mpiComm, 
                                                    true, std::vector<int>{}, bb_tol, newton_tol))
{
  nekrsCheck(!kernelsRegistered_,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "lpm_t::registerKernels has not been called prior to constructing lpm_t!");

  nekrsCheck(neknekCoupled(),
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "lpm_t + neknek is not supported!");

  nEXT = nrs->nEXT;
  nBDF = nrs->nBDF;

  dtEXT.resize(nEXT + 1);
  coeffEXT.resize(nEXT);
  o_coeffEXT = platform->device.malloc<dfloat>(nEXT);

  coeffRK.resize(std::max(solverOrder, bootstrapRKOrder));
  o_coeffRK = platform->device.malloc<dfloat>(coeffRK.size());

  coeffAB.resize(solverOrder);
  dt.resize(solverOrder + 1);
  o_coeffAB = platform->device.malloc<dfloat>(solverOrder);

  // coordinates are registered by default
  registerDOF("x");
  registerDOF("y");
  registerDOF("z");

  nStagesSumManyKernel = platform->kernelRequests.load("nStagesSumMany");
  remapParticlesKernel = platform->kernelRequests.load("lpm", "remapParticles");

  setTimerLevel(timerLevel);
  setTimerName(timerName);
}

void lpm_t::abOrder(int order)
{
  nekrsCheck(order <= 0 && order <= 3,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "Invalid integration order (%d)!\n",
             order);
  nekrsCheck(initialized_, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "lpm_t already initialized!");
  solverOrder = order;

  dt.resize(solverOrder + 1);
  coeffAB.resize(solverOrder);

  if (o_coeffAB.byte_size()) {
    o_coeffAB.free();
  }

  o_coeffAB = platform->device.malloc<dfloat>(solverOrder);
}

void lpm_t::rkOrder(int order)
{
  nekrsCheck(order <= 0,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "Integration order (%d) must be positive!\n",
             order);

  bool supported = false;
  constexpr int maxRKOrder = 4;
  for (int ord = 1; ord <= maxRKOrder; ++ord) {
    supported |= order == ord;
  }

  nekrsCheck(!supported, platform->comm.mpiComm, EXIT_FAILURE, "RK order (%d) is not supported!\n", order);
  nekrsCheck(initialized_, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "lpm_t already initialized!");
  solverOrder = order;

  dt.resize(solverOrder + 1);
  coeffRK.resize(solverOrder);

  if (o_coeffRK.byte_size()) {
    o_coeffRK.free();
  }

  o_coeffRK = platform->device.malloc<dfloat>(solverOrder);
}

lpm_t::SolverType lpm_t::stringToSolverType(const std::string &_solverType)
{
  std::string solverType = _solverType;
  lowerCase(solverType);
  if (solverType == "ab") {
    return SolverType::AB;
  }

  if (solverType == "rk") {
    return SolverType::RK;
  }

  return SolverType::INVALID;
}

void lpm_t::setSolver(const std::string &_solver)
{
  std::string solver = _solver;
  lowerCase(solver);
  this->solverType = stringToSolverType(solver);
  nekrsCheck(this->solverType == SolverType::INVALID,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "Solver (%s) is not supported.!",
             solver.c_str());
}

void lpm_t::registerDOF(const std::string &dofName, bool output)
{
  registerDOF(1, dofName, output);
}

void lpm_t::registerDOF(dlong Nfields, const std::string &_dofName, bool output)
{
  std::string dofName = _dofName;
  nekrsCheck(this->initialized(),
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "cannot register DOF %s after calling initialize!\n",
             dofName.c_str());

  lowerCase(dofName);
  const auto nDOFs = dofIds.size();
  if (dofIds.count(dofName) == 0) {
    dofNames.push_back(dofName);
    dofIds[dofName] = nDOFs;
    outputDofs[dofName] = output;
    dofCounts[dofName] = Nfields;
    nDOFs_ += Nfields;
    fieldType[dofName] = FieldType::DOF;
  }
}

int lpm_t::dofId(const std::string &_dofName) const
{
  std::string dofName = _dofName;
  lowerCase(dofName);
  nekrsCheck(dofIds.count(dofName) == 0,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "DOF %s not registered!\n",
             dofName.c_str());
  return dofIds.at(dofName);
}

int lpm_t::numDOFs(const std::string &_dofName) const
{
  std::string dofName = _dofName;
  lowerCase(dofName);
  nekrsCheck(dofIds.count(dofName) == 0,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "DOF %s not registered!\n",
             dofName.c_str());
  return dofCounts.at(dofName);
}

void lpm_t::registerProp(const std::string &propName, bool output)
{
  registerProp(1, propName, output);
}

void lpm_t::registerProp(dlong Nfields, const std::string &_propName, bool output)
{
  std::string propName = _propName;
  lowerCase(propName);
  nekrsCheck(this->initialized(),
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "cannot register prop %s after calling initialize!\n",
             propName.c_str());

  const auto nprops = propIds.size();
  if (propIds.count(propName) == 0) {
    propIds[propName] = nprops;
    outputProps[propName] = output;
    propCounts[propName] = Nfields;
    nProps_ += Nfields;
    fieldType[propName] = FieldType::PROP;
  }
}

int lpm_t::propId(const std::string &_propName) const
{
  std::string propName = _propName;
  lowerCase(propName);
  nekrsCheck(propIds.count(propName) == 0,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "prop %s not registered!\n",
             propName.c_str());
  return propIds.at(propName);
}

int lpm_t::numProps(const std::string &_propName) const
{
  std::string propName = _propName;
  lowerCase(propName);
  nekrsCheck(propIds.count(propName) == 0,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "prop %s not registered!\n",
             propName.c_str());
  return propCounts.at(propName);
}

void lpm_t::registerInterpField(const std::string &_interpFieldName,
                                int Nfields,
                                const occa::memory &o_fld,
                                bool output)
{
  std::string interpFieldName = _interpFieldName;
  lowerCase(interpFieldName);
  nekrsCheck(this->initialized(),
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "cannot register interpField %s after calling initialize!\n",
             interpFieldName.c_str());

  if (interpFieldIds.count(interpFieldName) == 0) {
    interpFieldIds[interpFieldName] = nInterpFields_;
    nInterpFields_ += Nfields;
    interpFieldCounts[interpFieldName] = Nfields;
    outputInterpFields[interpFieldName] = output;
    interpFieldInputs[interpFieldName] = o_fld;
    fieldType[interpFieldName] = FieldType::INTERP_FIELD;
  }
}

void lpm_t::registerInterpField(const std::string &interpFieldName, const occa::memory &o_fld, bool output)
{
  registerInterpField(interpFieldName, 1, o_fld, output);
}

int lpm_t::interpFieldId(const std::string &_interpFieldName) const
{
  std::string interpFieldName = _interpFieldName;
  lowerCase(interpFieldName);
  nekrsCheck(interpFieldIds.count(interpFieldName) == 0,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "interpField %s not registered!\n",
             interpFieldName.c_str());
  return interpFieldIds.at(interpFieldName);
}

int lpm_t::numFieldsInterp(const std::string &_interpFieldName) const
{
  std::string interpFieldName = _interpFieldName;
  lowerCase(interpFieldName);
  nekrsCheck(interpFieldIds.count(interpFieldName) == 0,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "interpField %s not registered!\n",
             interpFieldName.c_str());
  return interpFieldCounts.at(interpFieldName);
}

void lpm_t::setUserRHS(lpm_t::rhsFunc_t userRHS)
{
  userRHS_ = userRHS;
}

void lpm_t::addUserData(void *userdata)
{
  userdata_ = userdata;
}

occa::memory lpm_t::getDOF(const std::string &dofName) const
{
  const auto id = dofId(dofName);
  auto Nfields = numDOFs(dofName);

  return o_y.slice(id * fieldOffset_, Nfields * fieldOffset_);
}

const std::vector<dfloat> lpm_t::getDOFHost(const std::string &dofName)
{
  auto o_dof = getDOF(dofName);
  auto Nfields = numDOFs(dofName);

  std::vector<dfloat> h_dof(Nfields * fieldOffset_);
  o_dof.copyTo(h_dof.data(), Nfields * fieldOffset_);
  return h_dof;
}

occa::memory lpm_t::getProp(const std::string &propName) const
{
  const auto id = propId(propName);
  auto Nfields = numProps(propName);
  return o_prop.slice(id * fieldOffset_, Nfields * fieldOffset_);
}

std::vector<dfloat> lpm_t::getPropHost(const std::string &propName) const
{
  auto o_propEntry = getProp(propName);
  auto Nfields = numProps(propName);

  std::vector<dfloat> h_prop(Nfields * fieldOffset_);
  o_propEntry.copyTo(h_prop.data(), Nfields * fieldOffset_);
  return h_prop;
}

void lpm_t::setProp(const std::string &propName, const occa::memory &o_fld, dlong fldOffset)
{
  auto o_propEntry = getProp(propName);
  auto Nfields = numProps(propName);
  const auto offset = (fldOffset > 0) ? fldOffset : numParticles();

  if (o_fld.byte_size()) {
    o_propEntry.copyFrom(o_fld, Nfields * offset);
  }
}

const occa::memory lpm_t::getInterpField(const std::string &interpFieldName)
{
  const auto id = interpFieldId(interpFieldName);
  auto Nfields = numFieldsInterp(interpFieldName);
  return o_interpFld.slice(id * fieldOffset_, Nfields * fieldOffset_);
}

const std::vector<dfloat> lpm_t::getInterpFieldHost(const std::string &interpFieldName)
{
  auto o_interpFldEntry = getInterpField(interpFieldName);
  auto Nfields = numFieldsInterp(interpFieldName);

  std::vector<dfloat> h_interpField(Nfields * fieldOffset_);
  o_interpFldEntry.copyTo(h_interpField.data(), Nfields * fieldOffset_);
  return h_interpField;
}

void lpm_t::handleAllocation(size_t offset)
{
  o_y = platform->device.malloc<dfloat>(nDOFs_ * offset);
  o_ytmp = platform->device.malloc<dfloat>(nDOFs_ * offset);
  o_ydot = platform->device.malloc<dfloat>((solverOrder * nDOFs_) * offset);
  o_k = platform->device.malloc<dfloat>((std::max(solverOrder, bootstrapRKOrder) * nDOFs_) * offset);

  if (nProps_) {
    o_prop = platform->device.malloc<dfloat>(nProps_ * offset);
  }
  if (nInterpFields_) {
    o_interpFld = platform->device.malloc<dfloat>(nInterpFields_ * offset);
  }
}

void lpm_t::initialize(int nParticles, double t0, const std::vector<dfloat> &y0)
{
  nekrsCheck(initialized_, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "lpm_t already initialized!");
  nekrsCheck(y0.size() != nParticles * nDOFs_,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "y0.size() = %ld, while expecting %d entries!\n",
             y0.size(),
             nParticles * nDOFs_);

  auto o_y0 = platform->device.malloc<dfloat>(y0.size());
  o_y0.copyFrom(y0.data());
  this->initialize(nParticles, t0, o_y0);
}

void lpm_t::initialize(int nParticles, double t0, const occa::memory &o_y0)
{
  nekrsCheck(initialized_, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "lpm_t already initialized!");
  nekrsCheck(o_y0.length() != nParticles * nDOFs_,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "o_y0.length() = %llu , while expecting %d words!\n",
             o_y0.length(),
             nParticles * nDOFs_);

  time = t0;

  nParticles_ = nParticles;
  fieldOffset_ = alignStride<dfloat>(nParticles);

  handleAllocation(fieldOffset_);

  for (auto [fieldName, nFields] : interpFieldCounts) {
    laggedInterpFields[fieldName] = platform->device.malloc<dfloat>(nEXT * nFields * nrs->fieldOffset);
    extrapolatedInterpFields[fieldName] = platform->device.malloc<dfloat>(nFields * nrs->fieldOffset);
  }

  // set initial condition based on user-input
  if (nParticles_ > 0) {
    for (auto &dof : dofNames) {
      const auto id = dofId(dof);
      auto o_y_dof = getDOF(dof);
      auto o_y0_dof = o_y0 + id * nParticles_;
      o_y_dof.copyFrom(o_y0_dof, nParticles_);
    }
  }

  // do first findpts evaluation
  this->find(this->o_y);

  initialized_ = true;
}

void lpm_t::abCoeff(dfloat *dt, int tstep)
{
  const int order = std::min(tstep, this->solverOrder);
  nek::coeffAB(coeffAB.data(), dt, order);
  for (int i = 0; i < order; ++i) {
    coeffAB[i] *= dt[0];
  }
  for (int i = order; i > order; i--) {
    coeffAB[i - 1] = 0.0;
  }
  o_coeffAB.copyFrom(coeffAB.data(), solverOrder);
}

void lpm_t::interpolate()
{
  for (auto [interpFieldName, interpFieldId] : interpFieldInputs) {
    interpolate(interpFieldName);
  }
}

void lpm_t::interpolate(const std::string &interpFieldName)
{
  if (timerLevel != TimerLevel::None) {
    platform->timer.tic(timerName + "integrate::userRHS::interpolate", 1);
  }
  auto o_fld = extrapolatedInterpFields.at(interpFieldName);
  auto o_interpFld = getInterpField(interpFieldName);
  const auto Nfields = numFieldsInterp(interpFieldName);

  interp->setTimerName(timerName + "integrate::userRHS::interpolate::");
  deviceMemory<dfloat> d_fld(o_fld);
  deviceMemory<dfloat> d_interpFld(o_interpFld);
  interp->eval(Nfields, nrs->fieldOffset, d_fld, fieldOffset_, d_interpFld);
  if (timerLevel != TimerLevel::None) {
    platform->timer.toc(timerName + "integrate::userRHS::interpolate");
  }
}

void lpm_t::integrate(double tf)
{
  nekrsCheck(!initialized_,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "cannot call integrate before calling initialize!");
  nekrsCheck(!userRHS_,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "cannot call integrate without setting a userRHS!");

  if (timerLevel != TimerLevel::None) {
    platform->timer.tic(timerName + "integrate", 1);
  }

  // past integration time, exit without integrating
  // NOTE: need to also cache the initial condition of the user-added
  // interpolated fields
  if (time >= tf) {
    for (auto [fieldName, o_field] : laggedInterpFields) {
      const auto Nfields = numFieldsInterp(fieldName);
      auto o_currentField = interpFieldInputs.at(fieldName);
      o_field.copyFrom(o_currentField, Nfields * nrs->fieldOffset);
    }

    if (timerLevel != TimerLevel::None) {
      platform->timer.toc(timerName + "integrate");
    }
    return;
  }

  // set time step
  dfloat dtStep = tf - time;

  tstep++;

  if (platform->options.compareArgs("MOVING MESH", "TRUE")) {
    interp.reset();

    interp = std::make_unique<pointInterpolation_t>(nrs->mesh, platform->comm.mpiComm,
     true, std::vector<int>{}, bb_tol, newton_tol);
  }

  // set extrapolated state to t^n (copy from laggedInterpFields)
  for (auto [fieldName, o_field] : laggedInterpFields) {
    const auto Nfields = numFieldsInterp(fieldName);
    auto o_extField = extrapolatedInterpFields.at(fieldName);
    o_extField.copyFrom(o_field, Nfields * nrs->fieldOffset);
  }

  // set EXT dt's
  // lag previous time states
  for (int i = nEXT; i > 1; --i) {
    dtEXT[i] = dtEXT[i - 1];
  }
  dtEXT[0] = 0.0; // dtEXT[0] just used for extrapolation computation

  // lag previous time states
  for (int i = solverOrder; i > 0; --i) {
    dt[i] = dt[i - 1];
  }
  dt[0] = dtStep;

  if (userODESolver_) {
    deviceMemory<dfloat> o_y_(o_y);
    deviceMemory<dfloat> o_ydot_(o_ydot);
    userODESolver_(nrs, this, time, tf, tstep, o_y_, userdata_, o_ydot_);
  } else {
    if (solverType == SolverType::AB) {
      integrateAB();
    } else if (solverType == SolverType::RK) {
      integrateRK();
    }
  }

  dtEXT[1] = tf - time;
  time = tf;

  // lag previous time states in laggedInterpFields
  for (auto [fieldName, o_field] : laggedInterpFields) {
    const auto Nfields = numFieldsInterp(fieldName);
    const auto N = Nfields * nrs->fieldOffset;
    for (int s = nEXT; s > 1; s--) {
      o_field.copyFrom(o_field, N, (s - 1) * N, (s - 2) * N);
    }

    auto o_currentField = interpFieldInputs.at(fieldName);

    // update most recent time state
    o_field.copyFrom(o_currentField, N);
  }

  // always provide (t^n,y^n) for next step
  this->find(this->o_y);

  if (timerLevel != TimerLevel::None) {
    platform->timer.toc(timerName + "integrate");
  }
}

// setup new particle coordinates
void lpm_t::find(const occa::memory &o_yNew)
{
  occa::memory o_xCoord, o_yCoord, o_zCoord;
  if (fieldOffset_) {
    o_xCoord = o_yNew.slice(0 * fieldOffset_, nParticles_);
    o_yCoord = o_yNew.slice(1 * fieldOffset_, nParticles_);
    o_zCoord = o_yNew.slice(2 * fieldOffset_, nParticles_);
  }
  interp->setPoints(o_xCoord, o_yCoord, o_zCoord);
  platform->timer.tic(timerName + "integrate::find", 1);
  interp->setTimerName(timerName + "integrate::find::");
  interp->find(verbosityLevel);
  platform->timer.toc(timerName + "integrate::find");
}

// extrapolate fluid state to specified time state
void lpm_t::extrapolateFluidState(dfloat tEXT)
{
  const auto extOrder = std::min(tstep, nEXT);
  const auto bdfOrder = std::min(tstep, nBDF);

  auto mesh = nrs->mesh;

  std::copy(nrs->dt, nrs->dt + 3, dtEXT.begin());
  dtEXT[0] = tEXT - time;

  nek::extCoeff(coeffEXT.data(), dtEXT.data(), extOrder, bdfOrder);
  for (int i = nEXT; i > extOrder; i--) {
    coeffEXT[i - 1] = 0.0;
  }
  o_coeffEXT.copyFrom(coeffEXT.data(), nEXT);

  for (auto [fieldName, o_field] : laggedInterpFields) {
    const auto Nfields = numFieldsInterp(fieldName);
    auto o_extField = extrapolatedInterpFields.at(fieldName);
    nrs->extrapolateKernel(mesh->Nlocal, Nfields, nEXT, nrs->fieldOffset, o_coeffEXT, o_field, o_extField);
  }
}

void lpm_t::integrateAB()
{
  abCoeff(dt.data(), tstep);

  // lag derivatives
  if (nParticles_ > 0) {
    for (int s = solverOrder; s > 1; s--) {
      const auto N = nDOFs_ * fieldOffset_;
      o_ydot.copyFrom(o_ydot, N, (s - 1) * N, (s - 2) * N);
    }
  }

  // boostrap using RK method
  if (tstep <= solverOrder) {
    const auto saveSolverOrder = solverOrder;
    this->solverOrder = bootstrapRKOrder;
    integrateRK4();
    this->solverOrder = saveSolverOrder;
    if (fieldOffset_ > 0) {
      o_ydot.copyFrom(o_k, nDOFs_ * fieldOffset_); // for later lagging
    }
  } else {
    platform->timer.tic(timerName + "integrate::userRHS", 1);
    deviceMemory<dfloat> o_ydot_(o_ydot);
    deviceMemory<dfloat> o_y_(o_y);
    userRHS_(nrs, this, time, o_y_, userdata_, o_ydot_);
    platform->timer.toc(timerName + "integrate::userRHS");

    if (nParticles_ > 0) {
      nStagesSumManyKernel(nParticles_, fieldOffset_, solverOrder, nDOFs_, o_coeffAB, o_ydot, o_y);
    }
  }
}

void lpm_t::integrateRK()
{
  if (solverOrder == 1) {
    integrateRK1();
    return;
  }

  if (solverOrder == 2) {
    integrateRK2();
    return;
  }

  if (solverOrder == 3) {
    integrateRK3();
    return;
  }

  if (solverOrder == 4) {
    integrateRK4();
    return;
  }

  nekrsCheck(false, platform->comm.mpiComm, EXIT_FAILURE, "RK solver order %d not supported!\n", solverOrder);
}

void lpm_t::integrateRK1()
{
  platform->timer.tic(timerName + "integrate::userRHS", 1);
  deviceMemory<dfloat> o_k_(o_k);
  deviceMemory<dfloat> o_y_(o_y);
  userRHS_(nrs, this, time, o_y_, userdata_, o_k_);
  platform->timer.toc(timerName + "integrate::userRHS");

  dfloat rkCoeff = dt[0];

  // o_y += dt[0] * o_ydot
  if (nParticles_) {
    platform->linAlg->axpbyMany(nParticles_, nDOFs_, fieldOffset_, rkCoeff, o_k, 1.0, o_y);
  }
}

void lpm_t::integrateRK2()
{
  occa::memory o_k1, o_k2;
  if (fieldOffset_ > 0) {
    o_k1 = o_k + 0 * nDOFs_ * fieldOffset_;
    o_k2 = o_k + 1 * nDOFs_ * fieldOffset_;
  }
  deviceMemory<dfloat> o_k1_(o_k1);
  deviceMemory<dfloat> o_k2_(o_k2);
  deviceMemory<dfloat> o_y_(o_y);

  platform->timer.tic(timerName + "integrate::userRHS", 1);

  userRHS_(nrs, this, time, o_y_, userdata_, o_k1_);
  platform->timer.toc(timerName + "integrate::userRHS");

  // o_ytmp = o_y + dt[0] * o_k1
  if (nParticles_) {
    platform->linAlg->axpbyzMany(nParticles_, nDOFs_, fieldOffset_, 1.0, o_y, dt[0], o_k1, o_ytmp);
  }

  this->find(o_ytmp);

  // extrapolate to t^{n+1} using t^{n}, t^{n-1}, ...
  extrapolateFluidState(time + dt[0]);

  platform->timer.tic(timerName + "integrate::userRHS", 1);
  deviceMemory<dfloat> o_ytmp_(o_ytmp);
  userRHS_(nrs, this, time + dt[0], o_ytmp_, userdata_, o_k2_);
  platform->timer.toc(timerName + "integrate::userRHS");

  // o_y = o_y + 0.5 * dt[0] * (o_k1 + o_k2)
  coeffRK[0] = 0.5 * dt[0], coeffRK[1] = 0.5 * dt[0];
  o_coeffRK.copyFrom(coeffRK.data(), coeffRK.size());

  if (nParticles_) {
    this->nStagesSumManyKernel(nParticles_, fieldOffset_, solverOrder, nDOFs_, o_coeffRK, o_k, o_y);
  }
}

void lpm_t::integrateRK3()
{
  occa::memory o_k1, o_k2, o_k3;
  if (fieldOffset_ > 0) {
    o_k1 = o_k + 0 * nDOFs_ * fieldOffset_;
    o_k2 = o_k + 1 * nDOFs_ * fieldOffset_;
    o_k3 = o_k + 2 * nDOFs_ * fieldOffset_;
  }
  deviceMemory<dfloat> o_k1_(o_k1);
  deviceMemory<dfloat> o_k2_(o_k2);
  deviceMemory<dfloat> o_k3_(o_k3);
  deviceMemory<dfloat> o_y_(o_y);

  platform->timer.tic(timerName + "integrate::userRHS", 1);
  userRHS_(nrs, this, time, o_y_, userdata_, o_k1_);
  platform->timer.toc(timerName + "integrate::userRHS");

  // o_ytmp = o_y + dt[0]/2 * o_k1
  if (nParticles_) {
    platform->linAlg->axpbyzMany(nParticles_, nDOFs_, fieldOffset_, 1.0, o_y, 0.5 * dt[0], o_k1, o_ytmp);
  }

  this->find(o_ytmp);

  // extrapolate to t^{n+1/2} using t^{n}, t^{n-1}, ...
  extrapolateFluidState(time + 0.5 * dt[0]);

  platform->timer.tic(timerName + "integrate::userRHS", 1);
  userRHS_(nrs, this, time + 0.5 * dt[0], o_y_, userdata_, o_k2_);
  platform->timer.toc(timerName + "integrate::userRHS");

  // o_ytmp = o_y - dt[0] * o_k1 + 2 * dt[0] * o_k2
  if (nParticles_) {
    platform->linAlg->axpbyzMany(nParticles_, nDOFs_, fieldOffset_, 1.0, o_y, -dt[0], o_k1, o_ytmp);
    platform->linAlg->axpbyMany(nParticles_, nDOFs_, fieldOffset_, 2.0 * dt[0], o_k2, 1.0, o_ytmp);
  }

  this->find(o_ytmp);

  // extrapolate to t^{n+1} using t^{n}, t^{n-1}, ...
  extrapolateFluidState(time + dt[0]);

  platform->timer.tic(timerName + "integrate::userRHS", 1);
  deviceMemory<dfloat> o_ytmp_(o_ytmp);
  userRHS_(nrs, this, time + dt[0], o_ytmp_, userdata_, o_k3_);
  platform->timer.toc(timerName + "integrate::userRHS");

  // o_y = o_y + dt[0] * (1/6 * o_k1 + 2/3 * o_k2 + 1/6 * o_k3)
  coeffRK[0] = 1.0 / 6.0 * dt[0], coeffRK[1] = 2.0 / 3.0 * dt[0], coeffRK[2] = 1.0 / 6.0 * dt[0];
  o_coeffRK.copyFrom(coeffRK.data(), coeffRK.size());

  if (nParticles_) {
    this->nStagesSumManyKernel(nParticles_, fieldOffset_, solverOrder, nDOFs_, o_coeffRK, o_k, o_y);
  }
}

void lpm_t::integrateRK4()
{
  occa::memory o_k1, o_k2, o_k3, o_k4;
  if (fieldOffset_ > 0) {
    o_k1 = o_k + 0 * nDOFs_ * fieldOffset_;
    o_k2 = o_k + 1 * nDOFs_ * fieldOffset_;
    o_k3 = o_k + 2 * nDOFs_ * fieldOffset_;
    o_k4 = o_k + 3 * nDOFs_ * fieldOffset_;
  }
  deviceMemory<dfloat> o_k1_(o_k1);
  deviceMemory<dfloat> o_k2_(o_k2);
  deviceMemory<dfloat> o_k3_(o_k3);
  deviceMemory<dfloat> o_k4_(o_k4);
  deviceMemory<dfloat> o_y_(o_y);

  platform->timer.tic(timerName + "integrate::userRHS", 1);
  userRHS_(nrs, this, time, o_y_, userdata_, o_k1_);
  platform->timer.toc(timerName + "integrate::userRHS");

  // o_ytmp = o_y + dt[0]/2 * o_k1
  if (nParticles_) {
    platform->linAlg->axpbyzMany(nParticles_, nDOFs_, fieldOffset_, 1.0, o_y, 0.5 * dt[0], o_k1, o_ytmp);
  }

  this->find(o_ytmp);

  // extrapolate to t^{n+1/2} using t^{n}, t^{n-1}, ...
  extrapolateFluidState(time + 0.5 * dt[0]);

  platform->timer.tic(timerName + "integrate::userRHS", 1);
  userRHS_(nrs, this, time + 0.5 * dt[0], o_y_, userdata_, o_k2_);
  platform->timer.toc(timerName + "integrate::userRHS");

  // o_ytmp = o_y + dt[0]/2 * o_k2
  if (nParticles_) {
    platform->linAlg->axpbyzMany(nParticles_, nDOFs_, fieldOffset_, 1.0, o_y, 0.5 * dt[0], o_k2, o_ytmp);
  }

  this->find(o_ytmp);

  platform->timer.tic(timerName + "integrate::userRHS", 1);
  userRHS_(nrs, this, time + 0.5 * dt[0], o_y_, userdata_, o_k3_);
  platform->timer.toc(timerName + "integrate::userRHS");

  // o_ytmp = o_y + dt[0] * o_k3
  if (nParticles_) {
    platform->linAlg->axpbyzMany(nParticles_, nDOFs_, fieldOffset_, 1.0, o_y, dt[0], o_k3, o_ytmp);
  }

  this->find(o_ytmp);

  // extrapolate to t^{n+1} using t^{n}, t^{n-1}, ...
  extrapolateFluidState(time + dt[0]);

  platform->timer.tic(timerName + "integrate::userRHS", 1);
  deviceMemory<dfloat> o_ytmp_(o_ytmp);
  userRHS_(nrs, this, time + dt[0], o_ytmp_, userdata_, o_k4_);
  platform->timer.toc(timerName + "integrate::userRHS");

  // o_y = o_y + dt[0] * (1/6 * o_k1 + 1/3 * o_k2 + 1/3 * o_k3 + 1/6 * o_k4)
  coeffRK[0] = 1.0 / 6.0 * dt[0];
  coeffRK[1] = 1.0 / 3.0 * dt[0];
  coeffRK[2] = 1.0 / 3.0 * dt[0];
  coeffRK[3] = 1.0 / 6.0 * dt[0];
  o_coeffRK.copyFrom(coeffRK.data(), coeffRK.size());

  if (nParticles_) {
    this->nStagesSumManyKernel(nParticles_, fieldOffset_, solverOrder, nDOFs_, o_coeffRK, o_k, o_y);
  }
}

std::set<std::string> lpm_t::nonCoordinateOutputDOFs() const
{
  std::set<std::string> outputDofFields;
  for (auto [dofName, outputDof] : outputDofs) {
    if (outputDof) {
      outputDofFields.insert(dofName);
    }
  }

  outputDofFields.erase("x");
  outputDofFields.erase("y");
  outputDofFields.erase("z");

  return outputDofFields;
}

long long int lpm_t::numGlobalParticles() const
{
  long long int numGlobalParticles = this->numParticles();
  MPI_Allreduce(MPI_IN_PLACE, &numGlobalParticles, 1, MPI_LONG_LONG_INT, MPI_SUM, platform->comm.mpiComm);
  return numGlobalParticles;
}

int lpm_t::numNonLocalParticles() const
{
  auto &data = interp->data();

  int numNonLocal = 0;
  for (int i = 0; i < this->numParticles(); ++i) {
    numNonLocal += (data.proc[i] != platform->comm.mpiRank);
  }

  return numNonLocal;
}

int lpm_t::numUnfoundParticles() const
{
  auto &data = interp->data();

  int numUnfound = 0;
  for (int i = 0; i < this->numParticles(); ++i) {
    numUnfound += (data.code[i] == findpts::CODE_NOT_FOUND);
  }

  return numUnfound;
}

namespace
{
template <int N> struct particle_t_N {
  dfloat r[3];
  dlong code, proc, el;

  dfloat data[N];
};
} // namespace

template <int N>
void lpm_t::sendReceiveDataImpl(const std::vector<dfloat> &sendData,
                                const std::vector<dfloat> &r,
                                const std::vector<dlong> &proc,
                                const std::vector<dlong> &code,
                                const std::vector<dlong> &elem,
                                std::vector<dfloat> &recvData,
                                std::vector<dfloat> &recvR,
                                std::vector<dlong> &recvCode,
                                std::vector<dlong> &recvElem,
                                int entriesPerParticle)
{
  if (N != entriesPerParticle) {
    return; // nothing to do in this instantiation
  }

  const auto nSend = elem.size();

  using particle_t = particle_t_N<N>;

  struct array transfer;
  array_init(particle_t, &transfer, nSend);
  transfer.n = nSend;

  constexpr int dim = 3;

  particle_t *particles = (particle_t *)transfer.ptr;

  for (int pid = 0; pid < nSend; ++pid) {
    auto &p = particles[pid];
    for (int d = 0; d < dim; ++d) {
      p.r[d] = r[dim * pid + d];
    }
    p.code = code[pid];
    p.proc = proc[pid];
    p.el = elem[pid];
    for (int j = 0; j < N; ++j) {
      p.data[j] = sendData[N * pid + j];
    }
  }

  sarray_transfer(particle_t, &transfer, proc, true, this->interp->ptr()->crystalRouter());

  particles = (particle_t *)transfer.ptr;

  const auto nRecv = transfer.n;

  // allocate host buffers
  recvData.resize(nRecv * N);
  recvR.resize(nRecv * dim);
  recvCode.resize(nRecv);
  recvElem.resize(nRecv);

  // unpack from particles into host buffers
  for (int pid = 0; pid < nRecv; ++pid) {
    auto &p = particles[pid];
    for (int d = 0; d < dim; ++d) {
      recvR[dim * pid + d] = p.r[d];
    }
    recvCode[pid] = p.code;
    recvElem[pid] = p.el;
    for (int j = 0; j < N; ++j) {
      recvData[N * pid + j] = p.data[j];
    }
  }
}

std::tuple<std::vector<dfloat>, std::vector<dfloat>, std::vector<dlong>, std::vector<dlong>>
lpm_t::sendReceiveData(const std::vector<dfloat> &sendData,
                       const std::vector<dfloat> &r,
                       const std::vector<dlong> &proc,
                       const std::vector<dlong> &code,
                       const std::vector<dlong> &elem,
                       int entriesPerParticle)
{

  if (timerLevel != TimerLevel::None) {
    platform->timer.tic(timerName + "migrate::sendReceiveData", 1);
  }

  std::vector<dfloat> recvData;
  std::vector<dfloat> recvR;
  std::vector<dlong> recvCode;
  std::vector<dlong> recvElem;

  auto entries = n_tuple<int, lpm_t::maxEntriesPerParticleMigration>{};
  tuple_for_each(entries, [&](auto T) {
    sendReceiveDataImpl<decltype(T)::value>(sendData,
                                            r,
                                            proc,
                                            code,
                                            elem,
                                            recvData,
                                            recvR,
                                            recvCode,
                                            recvElem,
                                            entriesPerParticle);
  });

  if (timerLevel != TimerLevel::None) {
    platform->timer.toc(timerName + "migrate::sendReceiveData");
  }

  return std::make_tuple(recvData, recvR, recvCode, recvElem);
}

void lpm_t::migrate()
{

  const int entriesPerParticle = nDOFs_ + solverOrder * nDOFs_ + nProps_ + nInterpFields_;
  nekrsCheck(entriesPerParticle > lpm_t::maxEntriesPerParticleMigration,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "entriesPerParticle (%d) > lpm_t::maxEntriesPerParticleMigration (%d)!\n",
             entriesPerParticle,
             lpm_t::maxEntriesPerParticleMigration);

  if (timerLevel != TimerLevel::None) {
    platform->timer.tic(timerName + "migrate", 1);
  }

  if (timerLevel != TimerLevel::None) {
    platform->timer.tic(timerName + "migrate::find", 1);
  }

  {
    auto o_xcoord = getDOF("x");
    auto o_ycoord = getDOF("y");
    auto o_zcoord = getDOF("z");
    interp->setPoints(o_xcoord, o_ycoord, o_zcoord);
  }

  // disable findpts kernel timer for this call
  auto saveLevel = getTimerLevel();
  setTimerLevel(TimerLevel::None);
  interp->find(VerbosityLevel::None);
  setTimerLevel(saveLevel);
  if (timerLevel != TimerLevel::None) {
    platform->timer.toc(timerName + "migrate::find");
  }

  auto &data = interp->data();

  const auto nUnfound = numUnfoundParticles();
  const auto nNonLocal = numNonLocalParticles();

  // host allocations
  auto ySend = std::vector<dfloat>(nNonLocal * nDOFs_);
  auto ydotSend = std::vector<dfloat>(solverOrder * nNonLocal * nDOFs_);
  auto propSend = std::vector<dfloat>(nNonLocal * nProps_);
  auto interpFldSend = std::vector<dfloat>(nNonLocal * nInterpFields_);

  std::vector<dfloat> r(3 * nNonLocal);
  std::vector<dlong> proc(nNonLocal);
  std::vector<dlong> code(nNonLocal);
  std::vector<dlong> elem(nNonLocal);
  std::vector<dfloat> sendingData(ySend.size() + ydotSend.size() + propSend.size() + interpFldSend.size());

  // Pack data to send to other ranks
  if (timerLevel != TimerLevel::None) {
    platform->timer.tic(timerName + "migrate::packSending", 1);
  }
  if (nNonLocal) {
    // pack buffers for sending to other procs
    auto o_ySend = platform->device.malloc<dfloat>(nNonLocal * nDOFs_);
    auto o_ydotSend = platform->device.malloc<dfloat>(solverOrder * nNonLocal * nDOFs_);

    occa::memory o_propSend, o_interpFldSend;
    if (nProps_) {
      o_propSend = platform->device.malloc<dfloat>(nNonLocal * nProps_);
    }
    if (nInterpFields_) {
      o_interpFldSend = platform->device.malloc<dfloat>(nNonLocal * nInterpFields_);
    }

    std::vector<dlong> sendRankMap(this->numParticles(), -1);
    unsigned ctr = 0;
    for (int pid = 0; pid < this->numParticles(); ++pid) {
      if (data.proc[pid] != platform->comm.mpiRank && data.code[pid] != findpts::CODE_NOT_FOUND) {
        sendRankMap[pid] = ctr;
        r[3 * ctr + 0] = data.r[3 * pid + 0];
        r[3 * ctr + 1] = data.r[3 * pid + 1];
        r[3 * ctr + 2] = data.r[3 * pid + 2];
        proc[ctr] = data.proc[pid];
        code[ctr] = data.code[pid];
        elem[ctr] = data.el[pid];
        ctr++;
      }
    }

    if (o_sendRankMap.length() < sendRankMap.size()) {
      if (o_sendRankMap.length()) {
        o_sendRankMap.free();
      }
      o_sendRankMap = platform->device.malloc<dlong>(sendRankMap.size());
    }

    o_sendRankMap.copyFrom(sendRankMap.data(), sendRankMap.size());

    remapParticlesKernel(this->numParticles(),
                         fieldOffset_,
                         nNonLocal,
                         nProps_,
                         nInterpFields_,
                         nDOFs_,
                         solverOrder,
                         o_sendRankMap,
                         o_y,
                         o_ydot,
                         o_prop,
                         o_interpFld,
                         o_ySend,
                         o_ydotSend,
                         o_propSend,
                         o_interpFldSend);

    // copy to host
    o_ySend.copyTo(ySend.data(), nNonLocal * nDOFs_);
    o_ydotSend.copyTo(ydotSend.data(), solverOrder * nNonLocal * nDOFs_);
    if (nProps_) {
      o_propSend.copyTo(propSend.data(), nNonLocal * nProps_);
    }
    if (nInterpFields_) {
      o_interpFldSend.copyTo(interpFldSend.data(), nNonLocal * nInterpFields_);
    }
  }

  // copy host data into single long array, converting from packing as SoA to AoS structure
  unsigned ctr = 0;
  for (int pid = 0; pid < nNonLocal; ++pid) {
    for (int dof = 0; dof < nDOFs_; ++dof) {
      sendingData[ctr++] = ySend[pid + nNonLocal * dof];
    }
    for (int dof = 0; dof < nDOFs_; ++dof) {
      for (int s = 0; s < solverOrder; ++s) {
        sendingData[ctr++] = ydotSend[pid + nNonLocal * dof + s * nNonLocal * nDOFs_];
      }
    }
    for (int prop = 0; prop < nProps_; ++prop) {
      sendingData[ctr++] = propSend[pid + nNonLocal * prop];
    }
    for (int fld = 0; fld < nInterpFields_; ++fld) {
      sendingData[ctr++] = interpFldSend[pid + nNonLocal * fld];
    }
  }

  if (timerLevel != TimerLevel::None) {
    platform->timer.toc(timerName + "migrate::packSending");
  }

  auto [receivingData, recvR, recvCode, recvElem] =
      sendReceiveData(sendingData, r, proc, code, elem, entriesPerParticle);

  if (timerLevel != TimerLevel::None) {
    platform->timer.tic(timerName + "migrate::unpackReceiving", 1);
  }

  // copy data into host arrays, converting from packing as AoS to SoA structure
  const int nReceived = recvElem.size();

  auto yRecv = std::vector<dfloat>(nReceived * nDOFs_);
  auto ydotRecv = std::vector<dfloat>(solverOrder * nReceived * nDOFs_);
  auto propRecv = std::vector<dfloat>(nReceived * nProps_);
  auto interpFldRecv = std::vector<dfloat>(nReceived * nInterpFields_);

  ctr = 0;
  for (int pid = 0; pid < nReceived; ++pid) {
    for (int dof = 0; dof < nDOFs_; ++dof) {
      yRecv[pid + nReceived * dof] = receivingData[ctr++];
    }
    for (int dof = 0; dof < nDOFs_; ++dof) {
      for (int s = 0; s < solverOrder; ++s) {
        ydotRecv[pid + nReceived * dof + s * nReceived * nDOFs_] = receivingData[ctr++];
      }
    }
    for (int prop = 0; prop < nProps_; ++prop) {
      propRecv[pid + nReceived * prop] = receivingData[ctr++];
    }
    for (int fld = 0; fld < nInterpFields_; ++fld) {
      interpFldRecv[pid + nReceived * fld] = receivingData[ctr++];
    }
  }

  // start by deleting unfound particles -- local particles are stacked first
  const dlong newNParticles = this->numParticles() - nUnfound - nNonLocal + nReceived;
  const dlong newFieldOffset = alignStride<dfloat>(newNParticles);

  // allocate new fields
  auto o_propOld = this->o_prop;
  auto o_interpFldOld = this->o_interpFld;
  auto o_yOld = this->o_y;
  auto o_ydotOld = this->o_ydot;

  ctr = 0;
  if (newFieldOffset) {

    handleAllocation(newFieldOffset);

    if (this->numParticles()) {
      std::vector<dlong> migrateMap(this->numParticles(), -1);
      // map local particles to new location first
      for (int pid = 0; pid < this->numParticles(); ++pid) {
        if (data.proc[pid] == platform->comm.mpiRank && data.code[pid] != findpts::CODE_NOT_FOUND) {
          migrateMap[pid] = ctr;
          ctr++;
        }
      }

      if (o_migrateMap.length() < migrateMap.size()) {
        if (o_migrateMap.length()) {
          o_migrateMap.free();
        }
        o_migrateMap = platform->device.malloc<dlong>(migrateMap.size());
      }
      if (this->numParticles() > 0) {
        o_migrateMap.copyFrom(migrateMap.data(), migrateMap.size());
      }

      remapParticlesKernel(this->numParticles(),
                           fieldOffset_,
                           newFieldOffset,
                           nProps_,
                           nInterpFields_,
                           nDOFs_,
                           solverOrder,
                           o_migrateMap,
                           o_yOld,
                           o_ydotOld,
                           o_propOld,
                           o_interpFldOld,
                           this->o_y,
                           this->o_ydot,
                           this->o_prop,
                           this->o_interpFld);
    }
  }

  // Unpack received data into device arrays
  if (nReceived) {
    // pack buffers for sending to other procs
    auto o_yRecv = platform->device.malloc<dfloat>(nReceived * nDOFs_, yRecv.data());
    auto o_ydotRecv = platform->device.malloc<dfloat>(solverOrder * nReceived * nDOFs_, ydotRecv.data());

    occa::memory o_propRecv, o_interpFldRecv;
    if (nProps_) {
      o_propRecv = platform->device.malloc<dfloat>(nReceived * nProps_, propRecv.data());
    }
    if (nInterpFields_) {
      o_interpFldRecv = platform->device.malloc<dfloat>(nReceived * nInterpFields_, interpFldRecv.data());
    }

    std::vector<dlong> recvRankMap(nReceived, -1);
    for (int pid = 0; pid < nReceived; ++pid) {
      recvRankMap[pid] = ctr;
      ctr++;
    }

    if (o_recvRankMap.length() < recvRankMap.size()) {
      if (o_recvRankMap.length()) {
        o_recvRankMap.free();
      }
      o_recvRankMap = platform->device.malloc<dlong>(recvRankMap.size());
    }

    o_recvRankMap.copyFrom(recvRankMap.data(), recvRankMap.size());

    remapParticlesKernel(nReceived,
                         nReceived,
                         newFieldOffset,
                         nProps_,
                         nInterpFields_,
                         nDOFs_,
                         solverOrder,
                         o_recvRankMap,
                         o_yRecv,
                         o_ydotRecv,
                         o_propRecv,
                         o_interpFldRecv,
                         this->o_y,
                         this->o_ydot,
                         this->o_prop,
                         this->o_interpFld);
  }

  if (timerLevel != TimerLevel::None) {
    platform->timer.toc(timerName + "migrate::unpackReceiving");
  }

  nParticles_ = newNParticles;
  fieldOffset_ = newFieldOffset;

  if (o_yOld.byte_size()) {
    o_yOld.free();
  }
  if (o_ydotOld.byte_size()) {
    o_ydotOld.free();
  }
  if (o_propOld.byte_size()) {
    o_propOld.free();
  }
  if (o_interpFldOld.byte_size()) {
    o_interpFldOld.free();
  }

  // do an additional findpts call
  if (timerLevel != TimerLevel::None) {
    platform->timer.tic(timerName + "migrate::find", 1);
  }

  {
    auto o_xcoord = getDOF("x");
    auto o_ycoord = getDOF("y");
    auto o_zcoord = getDOF("z");
    interp->setPoints(o_xcoord, o_ycoord, o_zcoord);
  }

  // disable findpts kernel timer for this call
  saveLevel = getTimerLevel();
  setTimerLevel(TimerLevel::None);
  interp->find(VerbosityLevel::None);
  setTimerLevel(saveLevel);
  if (timerLevel != TimerLevel::None) {
    platform->timer.toc(timerName + "migrate::find");
  }

  if (timerLevel != TimerLevel::None) {
    platform->timer.toc(timerName + "migrate");
  }
}

void lpm_t::addParticles(int newNParticles,
                         const std::vector<dfloat> &yNewPart,
                         const std::vector<dfloat> &propNewPart)
{
  const auto expectedYdotSize = solverOrder * newNParticles * nDOFs_;
  std::vector<dfloat> ydotNewPart(expectedYdotSize, 0.0);
  addParticles(newNParticles, yNewPart, propNewPart, ydotNewPart);
}

void lpm_t::addParticles(int newNParticles,
                         const std::vector<dfloat> &yNewPart,
                         const std::vector<dfloat> &propNewPart,
                         const std::vector<dfloat> &ydotNewPart)
{
  const auto expectedYSize = newNParticles * nDOFs_;
  const auto expectedPropSize = newNParticles * nProps_;
  const auto expectedYdotSize = solverOrder * newNParticles * nDOFs_;
  nekrsCheck(yNewPart.size() < expectedYSize,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "yNewPart size is %ld but expected %d words!\n",
             yNewPart.size(),
             expectedYSize);
  nekrsCheck(propNewPart.size() < expectedPropSize,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "propNewPart size is %ld but expected %d words!\n",
             propNewPart.size(),
             expectedPropSize);
  nekrsCheck(ydotNewPart.size() < expectedYdotSize,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "ydotNewPart size is %ld but expected %d words!\n",
             ydotNewPart.size(),
             expectedYdotSize);

  auto o_yNewPart = platform->device.malloc<dfloat>(expectedYSize);
  o_yNewPart.copyFrom(yNewPart.data());
  auto o_propNewPart = platform->device.malloc<dfloat>(expectedPropSize);
  o_propNewPart.copyFrom(propNewPart.data());
  auto o_ydotNewPart = platform->device.malloc<dfloat>(expectedYdotSize);

  addParticles(newNParticles, o_yNewPart, o_propNewPart, o_ydotNewPart);

  o_yNewPart.free();
  o_propNewPart.free();
  o_ydotNewPart.free();
}

void lpm_t::addParticles(int newNParticles, const occa::memory &o_yNewPart, const occa::memory &o_propNewPart)
{
  const auto expectedYdotSize = solverOrder * newNParticles * nDOFs_;
  auto o_ydotNewPart = platform->device.malloc<dfloat>(expectedYdotSize);
  addParticles(newNParticles, o_yNewPart, o_propNewPart, o_ydotNewPart);
}

void lpm_t::addParticles(int newNParticles,
                         const occa::memory &o_yNewPart,
                         const occa::memory &o_propNewPart,
                         const occa::memory &o_ydotNewPart)
{
  if (timerLevel != TimerLevel::None) {
    platform->timer.tic(timerName + "addParticles", 1);
  }

  std::array<dlong, 2> counts = {newNParticles, this->numParticles()};
  MPI_Allreduce(MPI_IN_PLACE, counts.data(), 2, MPI_DLONG, MPI_SUM, platform->comm.mpiComm);

  if (platform->comm.mpiRank == 0) {
    std::cout << "Adding " << counts[0] << " to " << counts[1] << " particles!\n";
  }

  int incomingOffset = newNParticles;
  int newOffset = alignStride<dfloat>(this->nParticles_ + newNParticles);

  // check that the sizes of o_yNewPart, o_propNewPart are correct
  auto expectedYSize = incomingOffset * nDOFs_;
  auto expectedPropSize = incomingOffset * nProps_;
  auto expectedYdotSize = solverOrder * incomingOffset * nDOFs_;
  nekrsCheck(o_yNewPart.length() < expectedYSize,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "o_yNewPart length is %llu but expected %d words!\n",
             o_yNewPart.length(),
             expectedYSize);
  nekrsCheck(o_propNewPart.length() < expectedPropSize,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "o_propNewPart length is %llu but expected %d words!\n",
             o_propNewPart.length(),
             expectedPropSize);
  nekrsCheck(o_ydotNewPart.length() < expectedYdotSize,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "o_ydotNewPart length is %llu but expected %d words!\n",
             o_ydotNewPart.length(),
             expectedYdotSize);

  std::vector<dlong> remainingMap(this->nParticles_, 0);
  std::vector<dlong> insertMap(newNParticles, 0);

  // remainingMap[id] = id for existing particles
  if (this->nParticles_) {
    std::iota(remainingMap.begin(), remainingMap.end(), 0);
    if (o_remainingMap.length() < this->nParticles_) {
      if (o_remainingMap.length()) {
        o_remainingMap.free();
      }
      o_remainingMap = platform->device.malloc<dlong>(this->nParticles_);
    }
    o_remainingMap.copyFrom(remainingMap.data(), this->nParticles_);
  }

  // insertMap[id] = id + nParticles_ for incoming particles
  if (newNParticles) {
    std::iota(insertMap.begin(), insertMap.end(), this->nParticles_);
    if (o_insertMap.length() < newNParticles) {
      if (o_insertMap.length()) {
        o_insertMap.free();
      }
      o_insertMap = platform->device.malloc<dlong>(newNParticles);
    }
    o_insertMap.copyFrom(insertMap.data(), newNParticles);
  }

  auto o_propOld = this->o_prop;
  auto o_interpFldOld = this->o_interpFld;
  auto o_yOld = this->o_y;
  auto o_ydotOld = this->o_ydot;

  if (newOffset) {
    handleAllocation(newOffset);

    // dummy arrays for remapParticlesKernel for new particles
    occa::memory o_interpFldDummy;
    if (nInterpFields_) {
      o_interpFldDummy = platform->device.malloc<dfloat>(incomingOffset * nInterpFields_);
    }

    if (this->numParticles()) {
      // map existing particles to new data
      remapParticlesKernel(this->numParticles(),
                           fieldOffset_,
                           newOffset,
                           nProps_,
                           nInterpFields_,
                           nDOFs_,
                           solverOrder,
                           o_remainingMap,
                           o_yOld,
                           o_ydotOld,
                           o_propOld,
                           o_interpFldOld,
                           o_y,
                           o_ydot,
                           o_prop,
                           o_interpFld);
    }

    if (newNParticles) {
      // map new particles to new data
      remapParticlesKernel(newNParticles,
                           incomingOffset,
                           newOffset,
                           nProps_,
                           nInterpFields_,
                           nDOFs_,
                           solverOrder,
                           o_insertMap,
                           o_yNewPart,
                           o_ydotNewPart,
                           o_propNewPart,
                           o_interpFldDummy,
                           o_y,
                           o_ydot,
                           o_prop,
                           o_interpFld);
    }
    o_interpFldDummy.free();
  }

  if (o_propOld.byte_size()) {
    o_propOld.free();
  }
  if (o_interpFldOld.byte_size()) {
    o_interpFldOld.free();
  }
  if (o_yOld.byte_size()) {
    o_yOld.free();
  }
  if (o_ydotOld.byte_size()) {
    o_ydotOld.free();
  }

  nParticles_ += newNParticles;
  fieldOffset_ = newOffset;

  // update findpts results for newly added particles
  {
    auto o_xcoord = getDOF("x");
    auto o_ycoord = getDOF("y");
    auto o_zcoord = getDOF("z");
    interp->setPoints(o_xcoord, o_ycoord, o_zcoord);

    if (timerLevel != TimerLevel::None) {
      platform->timer.tic(timerName + "addParticles::find", 1);
    }

    // disable findpts kernel timer for this call
    auto saveLevel = getTimerLevel();
    setTimerLevel(TimerLevel::None);
    interp->find(VerbosityLevel::None);
    setTimerLevel(saveLevel);
    if (timerLevel != TimerLevel::None) {
      platform->timer.toc(timerName + "addParticles::find");
    }
  }

  if (timerLevel != TimerLevel::None) {
    platform->timer.toc(timerName + "addParticles");
  }
}

void lpm_t::deleteParticles()
{
  // apply tic here to get correct number of deletion events in timer output
  if (timerLevel != TimerLevel::None) {
    platform->timer.tic(timerName + "deleteParticles", 1);
  }

  int nDelete = numUnfoundParticles();
  long long int nDeleteGlobal = nDelete;
  MPI_Allreduce(MPI_IN_PLACE, &nDeleteGlobal, 2, MPI_LONG_LONG_INT, MPI_SUM, platform->comm.mpiComm);

  const auto nParticlesGlobal = numGlobalParticles();

  if (platform->comm.mpiRank == 0 && nDeleteGlobal > 0) {
    std::cout << "Deleting " << nDeleteGlobal << " of " << nParticlesGlobal << " particles!\n";
  }

  // Nothing to do on the current rank, exit early
  if (nDelete == 0) {
    if (timerLevel != TimerLevel::None) {
      platform->timer.toc(timerName + "deleteParticles");
    }
    return;
  }

  auto &code = interp->data().code;

  std::vector<dlong> remainingMap(this->numParticles(), -1);

  dlong newId = 0;
  for (int pid = 0; pid < this->numParticles(); ++pid) {
    if (code[pid] != findpts::CODE_NOT_FOUND) {
      remainingMap[pid] = newId;
      newId++;
    }
  }

  const auto Nwords = this->numParticles();
  if (Nwords > 0) {
    if (o_remainingMap.length() < Nwords) {
      if (o_remainingMap.length()) {
        o_remainingMap.free();
      }
      o_remainingMap = platform->device.malloc<dlong>(Nwords);
    }
    o_remainingMap.copyFrom(remainingMap.data(), Nwords);
  }

  auto o_propOld = this->o_prop;
  auto o_interpFldOld = this->o_interpFld;
  auto o_yOld = this->o_y;
  auto o_ydotOld = this->o_ydot;

  const auto newNParticles = this->numParticles() - nDelete;
  const auto newOffset = alignStride<dfloat>(newNParticles);

  if (newOffset) {
    handleAllocation(newOffset);

    remapParticlesKernel(this->numParticles(),
                         fieldOffset_,
                         newOffset,
                         nProps_,
                         nInterpFields_,
                         nDOFs_,
                         solverOrder,
                         o_remainingMap,
                         o_yOld,
                         o_ydotOld,
                         o_propOld,
                         o_interpFldOld,
                         o_y,
                         o_ydot,
                         o_prop,
                         o_interpFld);
  }

  if (o_propOld.byte_size()) {
    o_propOld.free();
  }
  if (o_interpFldOld.byte_size()) {
    o_interpFldOld.free();
  }
  if (o_yOld.byte_size()) {
    o_yOld.free();
  }
  if (o_ydotOld.byte_size()) {
    o_ydotOld.free();
  }

  nParticles_ = newNParticles;
  fieldOffset_ = newOffset;

  // update findpts results
  {
    auto o_xcoord = getDOF("x");
    auto o_ycoord = getDOF("y");
    auto o_zcoord = getDOF("z");
    interp->setPoints(o_xcoord, o_ycoord, o_zcoord);

    if (timerLevel != TimerLevel::None) {
      platform->timer.tic(timerName + "deleteParticles::find", 1);
    }

    // disable findpts kernel timer for this call
    auto saveLevel = getTimerLevel();
    setTimerLevel(TimerLevel::None);
    interp->find(VerbosityLevel::None);
    setTimerLevel(saveLevel);
    if (timerLevel != TimerLevel::None) {
      platform->timer.toc(timerName + "deleteParticles::find");
    }
  }

  if (timerLevel != TimerLevel::None) {
    platform->timer.toc(timerName + "deleteParticles");
  }
}

namespace
{
std::string lpm_vtu_data(const std::string &fieldName, long long int nComponent, long long int distance)
{
  return "<DataArray type=\"Float32\" Name=\"" + fieldName + "\" NumberOfComponents=\"" +
         std::to_string(nComponent) + "\" format=\"append\" offset=\"" + std::to_string(distance) + "\"/>\n";
}
} // namespace

void lpm_t::writeFld()
{
  if (timerLevel != TimerLevel::None) {
    platform->timer.tic(timerName + "write", 1);
  }

  // Required to determine if points are outside of the domain
  // Do not output points outside of the domain

  long long int nPartOutput = 0;
  auto &code = interp->data().code;
  {
    auto o_xcoord = getDOF("x");
    auto o_ycoord = getDOF("y");
    auto o_zcoord = getDOF("z");
    interp->setPoints(o_xcoord, o_ycoord, o_zcoord);

    // disable findpts kernel timer for this call
    auto saveLevel = getTimerLevel();
    setTimerLevel(TimerLevel::None);
    interp->find(VerbosityLevel::None);
    setTimerLevel(saveLevel);

    for (int pid = 0; pid < this->numParticles(); ++pid) {
      if (code[pid] != findpts::CODE_NOT_FOUND) {
        ++nPartOutput;
      }
    }
  }

  static int out_step = 0;
  ++out_step;

  MPI_Comm mpi_comm = platform->comm.mpiComm;

  long long int globalNPartOutput = nPartOutput;

  MPI_Allreduce(MPI_IN_PLACE, &globalNPartOutput, 1, MPI_LONG_LONG_INT, MPI_SUM, platform->comm.mpiComm);

  if (globalNPartOutput == 0) {
    if (platform->comm.mpiRank == 0) {
      std::cout << "No particles to output, skipping output step " << out_step << std::endl;
    }
    return;
  }

  std::ostringstream output;
  output << "par" << std::setw(5) << std::setfill('0') << out_step << ".vtu";
  std::string fname = output.str();

  long long int pOffset = 0;
  MPI_Exscan(&nPartOutput, &pOffset, 1, MPI_LONG_LONG_INT, MPI_SUM, mpi_comm);

  if (platform->comm.mpiRank == 0) {
    std::ofstream file(fname, std::ios::trunc);
    file.close();
  }

  MPI_File file_out;
  MPI_File_open(mpi_comm, fname.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_out);

  long long offset = 0;
  constexpr long long int dim = 3;

  // particles DOFs, sans coordinates
  auto particleOutputDOFs = nonCoordinateOutputDOFs();

  std::string message = "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" "
                        "header_type=\"UInt64\">\n";
  message += "\t<UnstructuredGrid>\n";
  message += "\t\t<FieldData>\n";
  message += "\t\t\t<DataArray type=\"Float32\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\"> " +
             std::to_string(time) + " </DataArray>\n";
  message += "\t\t\t<DataArray type=\"Int32\" Name=\"CYCLE\" NumberOfTuples=\"1\" format=\"ascii\"> " +
             std::to_string(out_step) + " </DataArray>\n";
  message += "\t\t</FieldData>\n";
  message += "\t\t<Piece NumberOfPoints=\"" + std::to_string(globalNPartOutput) + "\" NumberOfCells=\"0\">\n";
  message += "\t\t\t<Points>\n";
  message += "\t\t\t\t" + lpm_vtu_data("Position", dim, offset);
  offset += (dim * globalNPartOutput) * sizeof(float) + 1 * sizeof(long long int);
  message += "\t\t\t</Points>\n";

  message += "\t\t\t<PointData>\n";

  // output particle DOFs
  for (auto &&dofName : particleOutputDOFs) {
    const auto Nfields = dofCounts.at(dofName);
    message += "\t\t\t\t" + lpm_vtu_data(dofName, Nfields, offset);
    offset += (Nfields * globalNPartOutput) * sizeof(float) + 1 * sizeof(long long int);
  }

  // output particle properties
  for (auto [propName, isOutput] : outputProps) {
    if (!isOutput) {
      continue;
    }
    const auto Nfields = propCounts.at(propName);
    message += "\t\t\t\t" + lpm_vtu_data(propName, Nfields, offset);
    offset += (Nfields * globalNPartOutput) * sizeof(float) + 1 * sizeof(long long int);
  }

  // output interpolated fields
  for (auto [interpFieldName, isOutput] : outputInterpFields) {
    if (!isOutput) {
      continue;
    }
    const auto Nfields = interpFieldCounts.at(interpFieldName);
    message += "\t\t\t\t" + lpm_vtu_data(interpFieldName, Nfields, offset);
    offset += (Nfields * globalNPartOutput) * sizeof(float) + 1 * sizeof(long long int);
  }

  message += "\t\t\t</PointData>\n";

  message += "\t\t\t<Cells>\n";
  message += "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"/>\n";
  message += "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"/>\n";
  message += "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\"/>\n";
  message += "\t\t\t</Cells>\n";
  message += "\t\t</Piece>\n";
  message += "\t</UnstructuredGrid>\n";
  message += "\t<AppendedData encoding=\"raw\">\n";
  message += "_";

  if (platform->comm.mpiRank == 0) {
    MPI_File_write(file_out, message.c_str(), message.length(), MPI_CHAR, MPI_STATUS_IGNORE);
  }

  auto writeField = [&](int nFields, std::vector<float> &field) {
    MPI_Barrier(platform->comm.mpiComm);
    MPI_Offset position;
    MPI_File_get_size(file_out, &position);
    MPI_File_set_view(file_out, position, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
    if (platform->comm.mpiRank == 0) {
      unsigned long long int nbyte = nFields * globalNPartOutput * sizeof(float);
      MPI_File_write(file_out, &nbyte, 1, MPI_UNSIGNED_LONG_LONG, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(platform->comm.mpiComm);
    MPI_File_get_size(file_out, &position);

    position += sizeof(float) * (nFields * pOffset);
    MPI_File_set_view(file_out, position, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
    MPI_File_write_all(file_out, field.data(), field.size(), MPI_FLOAT, MPI_STATUS_IGNORE);
  };

  // output coordinates (required)
  {
    auto xCoord = getDOFHost("x");
    auto yCoord = getDOFHost("y");
    auto zCoord = getDOFHost("z");

    std::vector<float> positions(dim * nPartOutput, 0.0);
    dlong pid = 0;
    for (int particle = 0; particle < this->numParticles(); ++particle) {
      if (code[particle] != findpts::CODE_NOT_FOUND) {
        positions[dim * pid + 0] = static_cast<float>(xCoord[particle]);
        positions[dim * pid + 1] = static_cast<float>(yCoord[particle]);
        positions[dim * pid + 2] = static_cast<float>(zCoord[particle]);
        pid++;
      }
    }

    writeField(dim, positions);
  }

  // other particle DOFs
  for (auto &&dofName : particleOutputDOFs) {
    auto dofHost = getDOFHost(dofName);
    auto Nfields = numDOFs(dofName);

    std::vector<float> dofFloat(Nfields * nPartOutput, 0.0);
    dlong pid = 0;
    for (int particle = 0; particle < this->numParticles(); ++particle) {
      if (code[particle] != findpts::CODE_NOT_FOUND) {
        for (int fld = 0; fld < Nfields; ++fld) {
          dofFloat[Nfields * pid + fld] = static_cast<float>(dofHost[particle + fld * fieldOffset_]);
        }
        pid++;
      }
    }

    writeField(Nfields, dofFloat);
  }

  // particle properties
  for (auto [propName, isOutput] : outputProps) {
    if (!isOutput) {
      continue;
    }
    auto propHost = getPropHost(propName);
    auto Nfields = numProps(propName);

    std::vector<float> propFloat(Nfields * nPartOutput, 0.0);
    dlong pid = 0;
    for (int particle = 0; particle < this->numParticles(); ++particle) {
      if (code[particle] != findpts::CODE_NOT_FOUND) {
        for (int fld = 0; fld < Nfields; ++fld) {
          propFloat[Nfields * pid + fld] = static_cast<float>(propHost[particle + fld * fieldOffset_]);
        }
        pid++;
      }
    }

    writeField(Nfields, propFloat);
  }

  // interpolated fields
  for (auto [interpFieldName, isOutput] : outputInterpFields) {
    if (!isOutput) {
      continue;
    }
    auto interpFieldHost = getInterpFieldHost(interpFieldName);
    auto Nfields = numFieldsInterp(interpFieldName);

    std::vector<float> interpFieldFloat(Nfields * nPartOutput, 0.0);
    dlong pid = 0;
    for (int particle = 0; particle < this->numParticles(); ++particle) {
      if (code[particle] != findpts::CODE_NOT_FOUND) {
        for (int fld = 0; fld < Nfields; ++fld) {
          interpFieldFloat[Nfields * pid + fld] =
              static_cast<float>(interpFieldHost[particle + fld * fieldOffset_]);
        }
        pid++;
      }
    }

    writeField(Nfields, interpFieldFloat);
  }

  MPI_Barrier(platform->comm.mpiComm);
  MPI_Offset position;
  MPI_File_get_size(file_out, &position);
  MPI_File_set_view(file_out, position, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
  if (platform->comm.mpiRank == 0) {
    message = "";
    message += "</AppendedData>\n";
    message += "</VTKFile>";
    MPI_File_write(file_out, message.c_str(), message.length(), MPI_CHAR, MPI_STATUS_IGNORE);
  }

  MPI_File_close(&file_out);

  if (timerLevel != TimerLevel::None) {
    platform->timer.toc(timerName + "write");
  }
}

void lpm_t::registerKernels(occa::properties &kernelInfo)
{
  kernelsRegistered_ = true;

  std::string installDir(getenv("NEKRS_HOME"));
  // build kernels
  std::string fileName, kernelName;
  const std::string suffix = "Hex3D";
  const std::string oklpath(getenv("NEKRS_KERNEL_DIR"));

  fileName = oklpath + "/nrs/plugins/" + "lpm.okl";
  platform->kernelRequests.add("lpm", fileName, kernelInfo);
}

void lpm_t::setTimerLevel(TimerLevel level)
{
  timerLevel = level;
  interp->setTimerLevel(level);
}

TimerLevel lpm_t::getTimerLevel() const
{
  return timerLevel;
}

void lpm_t::setTimerName(const std::string &name)
{
  timerName = name;
  interp->setTimerName(name);
}

namespace
{
long long int parseNumParticles(const std::string &restartfile, const std::string &header)
{
  std::smatch npartmatch;
  bool found = std::regex_search(header, npartmatch, std::regex(R"(<Piece NumberOfPoints=\"(\d+)\")"));

  std::ostringstream errorLogger;
  long long int nparticles = -1;

  if (!found) {
    errorLogger << "Could not read number of particles while reading " << restartfile << "!\n";
  }

  try {
    nparticles = std::stoll(npartmatch[1].str());
  } catch (std::invalid_argument e) {
    errorLogger << "Could not read number of particles while reading " << restartfile << "!\n";
    errorLogger << "Exception said:\n" << e.what() << std::endl;
  }

  auto errorString = errorLogger.str();
  int errorLength = errorString.length();
  MPI_Allreduce(MPI_IN_PLACE, &errorLength, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);

  nekrsCheck(errorLength > 0, platform->comm.mpiComm, EXIT_FAILURE, "%s", errorString.c_str());

  return nparticles;
}

auto parsePointData(const std::string &restartfile, const std::string &pointData)
{
  std::ostringstream errorLogger;
  std::string fieldName = "";
  long long int numComponents = -1;
  long long int offset = -1;

  std::smatch match;
  bool found = std::regex_search(pointData, match, std::regex(R"(\s*Name=\"(.+?)\")"));
  if (!found) {
    errorLogger << "Could not parse pointData while reading " << restartfile << "!\n";
  }
  fieldName = match[1].str();

  found = std::regex_search(pointData, match, std::regex(R"(\s*NumberOfComponents=\"(\d+)\")"));
  if (!found) {
    errorLogger << "Could not parse " << fieldName << " number of components while reading " << restartfile
                << "\n";
  }

  try {
    numComponents = std::stoll(match[1].str());
  } catch (std::invalid_argument &e) {
    errorLogger << "Could not parse " << fieldName << " number of components while reading " << restartfile
                << "!\n";
    errorLogger << "Exception said:\n" << e.what() << std::endl;
  }

  found = std::regex_search(pointData, match, std::regex(R"(\s*offset=\"(\d+)\")"));
  if (!found) {
    errorLogger << "Could not parse " << fieldName << " offset while reading " << restartfile << "!\n";
  }

  try {
    offset = std::stoll(match[1].str());
  } catch (std::invalid_argument &e) {
    errorLogger << "Could not parse " << fieldName << " offset while reading " << restartfile << "!\n";
    errorLogger << "Exception said:\n" << e.what() << std::endl;
  }

  auto errorString = errorLogger.str();
  int errorLength = errorString.length();
  MPI_Allreduce(MPI_IN_PLACE, &errorLength, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);

  nekrsCheck(errorLength > 0, platform->comm.mpiComm, EXIT_FAILURE, "%s", errorString.c_str());

  return std::make_tuple(fieldName, numComponents, offset);
}

auto readHeader(const std::string &restartFile)
{
  // read header of VTK UnstructuredGrid file format until reading after the <AppendedData encoding=\"raw\">
  // line
  std::string header;

  // associated with DataArray's inside <PointData> tag
  std::vector<std::string> pointData;

  // read metadata from restart file
  std::ifstream file(restartFile);
  std::string line;
  bool insidePointData = false;

  while (std::getline(file, line)) {
    header += line + "\n";
    if (line.find("<AppendedData encoding=\"raw\">") != std::string::npos) {
      break;
    }

    if (line.find("<PointData>") != std::string::npos) {
      insidePointData = true;
    }

    if (line.find("</PointData>") != std::string::npos) {
      insidePointData = false;
    }

    // gather PointData attributes to read later
    if (insidePointData && line.find("<DataArray") != std::string::npos) {
      pointData.push_back(line);
    }
  }
  file.close();

  header += "_";

  return std::make_tuple(header, pointData);
}

} // namespace

void lpm_t::restart(const std::string &restartFile)
{
  bool fileExists = fs::exists(restartFile);
  nekrsCheck(!fileExists,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "Restart file %s does not exist!\n",
             restartFile.c_str());

  constexpr long long int dim = 3;
  auto [header, pointData] = readHeader(restartFile);

  // from header, extract number of particles stored in NumberOfPoints
  const auto nPartGlobal = parseNumParticles(restartFile, header);
  long long int nPartLocal = nPartGlobal / platform->comm.mpiCommSize;
  // distribute remaining particles
  const auto remainder = nPartGlobal % platform->comm.mpiCommSize;
  if (platform->comm.mpiRank < remainder) {
    nPartLocal++;
  }

  // initial conditions are read from VTK file below
  // pass in zeros at the moment -- these will be overwritten
  {
    std::vector<dfloat> dummy_y0(nPartLocal * this->nDOFs(), 0.0);
    double t0;
    platform->options.getArgs("START TIME", t0);
    this->initialize(nPartLocal, t0, dummy_y0);
  }

  long long int pOffset = 0;
  MPI_Exscan(&nPartLocal, &pOffset, 1, MPI_LONG_LONG_INT, MPI_SUM, platform->comm.mpiComm);

  std::map<std::string, std::tuple<long long int, long long int>> fieldToInfo;
  for (auto &&field : pointData) {
    auto [fieldName, numComponents, offset] = parsePointData(restartFile, field);
    fieldToInfo[fieldName] = std::make_tuple(numComponents, offset);
  }

  // start by reading coordinates, starting at the position left by the header
  MPI_File file_in;
  MPI_File_open(platform->comm.mpiComm, restartFile.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_in);

  MPI_Offset position = header.length();
  MPI_File_set_view(file_in, position, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

  // first field is number of bytes in coordinate data
  long long int nPointData = 0;
  MPI_File_read(file_in, &nPointData, 1, MPI_LONG_LONG_INT, MPI_STATUS_IGNORE);
  nPointData /= dim;
  nPointData /= sizeof(float);

  nekrsCheck(nPointData != nPartGlobal,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "Number of particles in header (%lld) does not match number of particles in file (%lld)",
             nPartGlobal,
             nPointData);

  position = header.length();
  position += sizeof(float) * (dim * pOffset) + 1 * sizeof(long long int);
  MPI_File_set_view(file_in, position, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

  std::vector<float> coords(nPartLocal * dim);

  std::vector<double> xCoord(nPartLocal);
  std::vector<double> yCoord(nPartLocal);
  std::vector<double> zCoord(nPartLocal);

  MPI_File_read(file_in, coords.data(), nPartLocal * dim, MPI_FLOAT, MPI_STATUS_IGNORE);
  for (int pid = 0; pid < nPartLocal; ++pid) {
    xCoord[pid] = static_cast<dfloat>(coords[dim * pid + 0]);
    yCoord[pid] = static_cast<dfloat>(coords[dim * pid + 1]);
    zCoord[pid] = static_cast<dfloat>(coords[dim * pid + 2]);
  }

  auto o_xCoord = getDOF("x");
  auto o_yCoord = getDOF("y");
  auto o_zCoord = getDOF("z");

  o_xCoord.copyFrom(xCoord.data(), nPartLocal);
  o_yCoord.copyFrom(yCoord.data(), nPartLocal);
  o_zCoord.copyFrom(zCoord.data(), nPartLocal);

  auto readField = [&, &header = header](std::string fieldName,
                                         long long int expectedNumComponents,
                                         long long int offset) {
    nekrsCheck(fieldType.count(fieldName) == 0,
               platform->comm.mpiComm,
               EXIT_FAILURE,
               "Encountered unregisterd field %s while reading restart %s!\n",
               fieldName.c_str(),
               restartFile.c_str());

    auto type = fieldType[fieldName];

    long long int nComponents = -1;
    occa::memory o_fld;
    if (type == FieldType::DOF) {
      nComponents = dofCounts.at(fieldName);
      o_fld = getDOF(fieldName);
    } else if (type == FieldType::PROP) {
      nComponents = propCounts.at(fieldName);
      o_fld = getProp(fieldName);
    } else if (type == FieldType::INTERP_FIELD) {
      nComponents = interpFieldCounts.at(fieldName);
      o_fld = getInterpField(fieldName);
    }

    nekrsCheck(
        nComponents != expectedNumComponents,
        platform->comm.mpiComm,
        EXIT_FAILURE,
        "Expected number of components for field %s (%lld) does not match number of components (%lld) in "
        "restart file %s!\n",
        fieldName.c_str(),
        expectedNumComponents,
        nComponents,
        restartFile.c_str());

    position = header.length();
    position += offset;

    MPI_File_set_view(file_in, position, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

    // first field is number of bytes in coordinate data
    long long int nPointData = 0;
    MPI_File_read(file_in, &nPointData, 1, MPI_LONG_LONG_INT, MPI_STATUS_IGNORE);
    nPointData /= nComponents;
    nPointData /= sizeof(float);

    nekrsCheck(nPointData != nPartGlobal,
               platform->comm.mpiComm,
               EXIT_FAILURE,
               "Number of particles in header (%lld) does not match number of particles in file (%lld) when "
               "reading field %s!\n",
               nPartGlobal,
               nPointData,
               fieldName.c_str());

    std::vector<float> fld(nPartLocal * nComponents);
    std::vector<dfloat> fldHost(this->fieldOffset() * nComponents);

    position = header.length();
    position += offset;
    position += sizeof(float) * (nComponents * pOffset) + 1 * sizeof(long long int);

    MPI_File_set_view(file_in, position, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

    MPI_File_read(file_in, fld.data(), nPartLocal * nComponents, MPI_FLOAT, MPI_STATUS_IGNORE);
    for (int pid = 0; pid < nPartLocal; ++pid) {
      for (int component = 0; component < nComponents; ++component) {
        fldHost[pid + component * this->fieldOffset()] =
            static_cast<dfloat>(fld[nComponents * pid + component]);
      }
    }

    o_fld.copyFrom(fldHost.data(), this->fieldOffset() * nComponents);
  };

  for (auto &&[fieldName, info] : fieldToInfo) {
    auto [numComponents, offset] = info;
    readField(fieldName, numComponents, offset);
  }

  MPI_File_close(&file_in);
}

void lpm_t::printTimers()
{
  const auto timerTags = platform->timer.tags();

  // filter out tags that do not start with timerName
  std::vector<std::string> filteredTags;
  std::copy_if(timerTags.begin(),
               timerTags.end(),
               std::back_inserter(filteredTags),
               [&](const std::string &tag) { return tag.find(timerName) == 0; });

  // construct tree where parent entries are the portion of the tag left of the last '::'
  std::map<std::string, std::vector<std::string>> tree;

  for (auto &&tag : filteredTags) {
    auto pos = tag.rfind("::");
    if (pos == std::string::npos) {
      tree[""].push_back(tag);
    } else {
      auto parent = tag.substr(0, pos);
      tree[parent].push_back(tag);
    }
  }

  // print tree
  std::function<void(std::string, int)> printTree;
  printTree = [&](std::string tag, int level) {
    if (level > 0) {
      const auto tTag = platform->timer.query(tag, "DEVICE:MAX");
      if (platform->comm.mpiRank == 0) {
        for (int i = 0; i < level; ++i) {
          std::cout << "> ";
        }
        std::cout << tag << " " << tTag << "s\n";
      }
    }

    for (auto &&child : tree[tag]) {
      printTree(child, level + 1);
    }
  };

  auto pos = timerName.rfind("::");
  const auto start = timerName.substr(0, pos);

  if (platform->comm.mpiRank == 0) {
    std::cout << "\n";
    std::cout << "Detailed timers for particles " << start << ":\n";
  }

  printTree(start, 0);

  if (platform->comm.mpiRank == 0) {
    std::cout << "\n";
  }
}

void lpm_t::resetTimers()
{
  const auto timerTags = platform->timer.tags();

  // filter out tags that do not start with timerName
  std::vector<std::string> filteredTags;
  std::copy_if(timerTags.begin(),
               timerTags.end(),
               std::back_inserter(filteredTags),
               [&](const std::string &tag) { return tag.find(timerName) == 0; });

  for (auto &&tag : filteredTags) {
    platform->timer.reset(tag);
  }
}

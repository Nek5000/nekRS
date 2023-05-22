#include "nekInterfaceAdapter.hpp" // for nek::coeffAB
#include "lpm.hpp"
#include "neknek.hpp"
#include "nrs.hpp"
#include "pointInterpolation.hpp"
#include <algorithm>
#include <cstdlib>
#include <numeric>
#include <regex>
#include <tuple>
#include <filesystem>
#include "tuple_for_each.hpp"
#include "gslib.h" // needed for sarray_transfer

#include <inttypes.h>

namespace {

int computeFieldOffset(int n)
{
  auto offset = n;
  const int pageW = ALIGN_SIZE / sizeof(dfloat);
  if (offset % pageW)
    offset = (offset / pageW + 1) * pageW;
  return offset;
}
} // namespace

lpm_t::lpm_t(nrs_t *nrs_, dfloat bb_tol_, dfloat newton_tol_)
    : nrs(nrs_), solverOrder(nrs->nEXT), bb_tol(bb_tol_), newton_tol(newton_tol_),
      interp(std::make_unique<pointInterpolation_t>(nrs, bb_tol, newton_tol))
{
  nrsCheck(!kernelsRegistered_,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "%s\n",
           "lpm_t::registerKernels has not been called prior to constructing lpm_t!");
  
  nrsCheck(neknekCoupled(),
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "%s\n",
           "lpm_t + neknek is not supported!");

  nEXT = nrs->nEXT;
  nBDF = nrs->nBDF;

  dtEXT.resize(nEXT + 1);
  coeffEXT.resize(nEXT);
  o_coeffEXT = platform->device.malloc(nEXT * sizeof(dfloat));

  coeffRK.resize(std::max(solverOrder, bootstrapRKOrder));
  o_coeffRK = platform->device.malloc(coeffRK.size() * sizeof(dfloat));

  coeffAB.resize(solverOrder);
  dt.resize(solverOrder + 1);
  o_coeffAB = platform->device.malloc(solverOrder * sizeof(dfloat));

  // coordinates are registered by default
  registerDOF("x");
  registerDOF("y");
  registerDOF("z");

  nStagesSumManyKernel = platform->kernels.get("nStagesSumMany");
  remapParticlesKernel = platform->kernels.get("remapParticles");

  setTimerLevel(timerLevel);
  setTimerName(timerName);
}

void lpm_t::abOrder(int order)
{
  nrsCheck(order <= 0,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "Integration order (%d) must be positive!\n",
           order);
  nrsCheck(initialized_, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "lpm_t already initialized!");
  solverOrder = order;

  dt.resize(solverOrder + 1);
  coeffAB.resize(solverOrder);

  if (o_coeffAB.size()) {
    o_coeffAB.free();
  }

  o_coeffAB = platform->device.malloc(solverOrder * sizeof(dfloat));
}

void lpm_t::rkOrder(int order)
{
  nrsCheck(order <= 0,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "Integration order (%d) must be positive!\n",
           order);

  bool supported = false;
  constexpr int maxRKOrder = 4;
  for (int ord = 1; ord <= maxRKOrder; ++ord) {
    supported |= order == ord;
  }

  nrsCheck(!supported,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "RK order (%d) is not supported!\n",
           order);
  nrsCheck(initialized_, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "lpm_t already initialized!");
  solverOrder = order;

  dt.resize(solverOrder + 1);
  coeffRK.resize(solverOrder);

  if (o_coeffRK.size()) {
    o_coeffRK.free();
  }

  o_coeffRK = platform->device.malloc(solverOrder * sizeof(dfloat));
}

lpm_t::SolverType lpm_t::stringToSolverType(std::string solverType)
{
  lowerCase(solverType);
  if (solverType == "ab")
    return SolverType::AB;

  if (solverType == "rk")
    return SolverType::RK;

  return SolverType::INVALID;
}

void lpm_t::setSolver(std::string solver)
{
  lowerCase(solver);
  this->solverType = stringToSolverType(solver);
  nrsCheck(this->solverType == SolverType::INVALID,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "Solver (%s) is not supported.!",
           solver.c_str());
}

void lpm_t::registerDOF(std::string dofName, bool output) { registerDOF(1, dofName, output); }

void lpm_t::registerDOF(dlong Nfields, std::string dofName, bool output)
{
  nrsCheck(this->initialized(),
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "cannot register DOF %s after calling initialize!\n",
           dofName.c_str());

  lowerCase(dofName);
  const auto nDOFs = dofIds.size();
  if (dofIds.count(dofName) == 0) {
    dofIds[dofName] = nDOFs;
    outputDofs[dofName] = output;
    dofCounts[dofName] = Nfields;
    nDOFs_ += Nfields;
    fieldType[dofName] = FieldType::DOF;
  }
}

int lpm_t::dofId(std::string dofName) const
{
  lowerCase(dofName);
  nrsCheck(dofIds.count(dofName) == 0,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "DOF %s not registered!\n",
           dofName.c_str());
  return dofIds.at(dofName);
}

int lpm_t::numDOFs(std::string dofName) const
{
  lowerCase(dofName);
  nrsCheck(dofIds.count(dofName) == 0,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "DOF %s not registered!\n",
           dofName.c_str());
  return dofCounts.at(dofName);
}

void lpm_t::registerProp(std::string propName, bool output) { registerProp(1, propName, output); }

void lpm_t::registerProp(dlong Nfields, std::string propName, bool output)
{
  lowerCase(propName);
  nrsCheck(this->initialized(),
           platform->comm.mpiComm,
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

int lpm_t::propId(std::string propName) const
{
  lowerCase(propName);
  nrsCheck(propIds.count(propName) == 0,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "prop %s not registered!\n",
           propName.c_str());
  return propIds.at(propName);
}

int lpm_t::numProps(std::string propName) const
{
  lowerCase(propName);
  nrsCheck(propIds.count(propName) == 0,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "prop %s not registered!\n",
           propName.c_str());
  return propCounts.at(propName);
}

void lpm_t::registerInterpField(std::string interpFieldName, int Nfields, occa::memory o_fld, bool output)
{
  lowerCase(interpFieldName);
  nrsCheck(this->initialized(),
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "cannot register interpField %s after calling initialize!\n",
           interpFieldName.c_str());

  const auto nInterpFields = interpFieldIds.size();
  if (interpFieldIds.count(interpFieldName) == 0) {
    interpFieldIds[interpFieldName] = nInterpFields;
    interpFieldCounts[interpFieldName] = Nfields;
    outputInterpFields[interpFieldName] = output;
    interpFieldInputs[interpFieldName] = o_fld;
    nInterpFields_ += Nfields;
    fieldType[interpFieldName] = FieldType::INTERP_FIELD;
  }
}

void lpm_t::registerInterpField(std::string interpFieldName, occa::memory o_fld, bool output)
{
  registerInterpField(interpFieldName, 1, o_fld, output);
}

int lpm_t::interpFieldId(std::string interpFieldName) const
{
  lowerCase(interpFieldName);
  nrsCheck(interpFieldIds.count(interpFieldName) == 0,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "interpField %s not registered!\n",
           interpFieldName.c_str());
  return interpFieldIds.at(interpFieldName);
}

int lpm_t::numFieldsInterp(std::string interpFieldName) const
{
  lowerCase(interpFieldName);
  nrsCheck(interpFieldIds.count(interpFieldName) == 0,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "interpField %s not registered!\n",
           interpFieldName.c_str());
  return interpFieldCounts.at(interpFieldName);
}

void lpm_t::setUserRHS(lpm_t::rhsFunc_t userRHS) { userRHS_ = userRHS; }

void lpm_t::addUserData(void *userdata) { userdata_ = userdata; }

occa::memory lpm_t::getDOF(int dofId)
{
  if (fieldOffset_ == 0)
    return o_y;
  return o_y + dofId * fieldOffset_ * sizeof(dfloat);
}
occa::memory lpm_t::getDOF(std::string dofName) { return getDOF(dofId(dofName)); }

std::vector<dfloat> lpm_t::getDOFHost(std::string dofName)
{
  auto o_dof = getDOF(dofName);
  auto Nfields = numDOFs(dofName);

  std::vector<dfloat> h_dof(Nfields * fieldOffset_);
  if (fieldOffset_ == 0)
    return h_dof;
  o_dof.copyTo(h_dof.data(), Nfields * fieldOffset_ * sizeof(dfloat));
  return h_dof;
}

occa::memory lpm_t::getProp(int propId)
{
  if (fieldOffset_ == 0)
    return o_prop;
  return o_prop + propId * fieldOffset_ * sizeof(dfloat);
}
occa::memory lpm_t::getProp(std::string propName) { return getProp(propId(propName)); }

std::vector<dfloat> lpm_t::getPropHost(std::string propName)
{
  auto o_propEntry = getProp(propName);
  auto Nfields = numProps(propName);

  std::vector<dfloat> h_prop(Nfields * fieldOffset_);
  if (fieldOffset_ == 0)
    return h_prop;
  o_propEntry.copyTo(h_prop.data(), Nfields * fieldOffset_ * sizeof(dfloat));
  return h_prop;
}

occa::memory lpm_t::getInterpField(int interpFieldId)
{
  if (fieldOffset_ == 0)
    return o_interpFld;
  return o_interpFld + interpFieldId * fieldOffset_ * sizeof(dfloat);
}
occa::memory lpm_t::getInterpField(std::string interpFieldName)
{
  return getInterpField(interpFieldId(interpFieldName));
}

std::vector<dfloat> lpm_t::getInterpFieldHost(std::string interpFieldName)
{
  auto o_interpFldEntry = getInterpField(interpFieldName);
  auto Nfields = numFieldsInterp(interpFieldName);

  std::vector<dfloat> h_interpField(Nfields * fieldOffset_);
  if (fieldOffset_ == 0)
    return h_interpField;
  o_interpFldEntry.copyTo(h_interpField.data(), Nfields * fieldOffset_ * sizeof(dfloat));
  return h_interpField;
}

void lpm_t::handleAllocation(int offset)
{
  o_y = platform->device.malloc(offset * nDOFs_ * sizeof(dfloat));
  o_ytmp = platform->device.malloc(offset * nDOFs_ * sizeof(dfloat));
  o_ydot = platform->device.malloc(solverOrder * offset * nDOFs_ * sizeof(dfloat));
  o_k = platform->device.malloc(std::max(solverOrder, bootstrapRKOrder) * offset * nDOFs_ * sizeof(dfloat));

  if (nProps_) {
    o_prop = platform->device.malloc(offset * nProps_ * sizeof(dfloat));
  }
  if (nInterpFields_) {
    o_interpFld = platform->device.malloc(offset * nInterpFields_ * sizeof(dfloat));
  }
}

void lpm_t::initialize(int nParticles, dfloat t0, std::vector<dfloat> &y0)
{
  nrsCheck(initialized_, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "lpm_t already initialized!");
  nrsCheck(y0.size() != nParticles * nDOFs_,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "y0.size() = %ld, while expecting %d entries!\n",
           y0.size(),
           nParticles * nDOFs_);

  auto o_y0 = platform->device.malloc(y0.size() * sizeof(dfloat), y0.data());
  this->initialize(nParticles, t0, o_y0);
}

void lpm_t::initialize(int nParticles, dfloat t0, occa::memory o_y0)
{
  nrsCheck(initialized_, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "lpm_t already initialized!");
  nrsCheck(o_y0.size() != nParticles * nDOFs_ * sizeof(dfloat),
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "o_y0.size() = %" PRId64 ", while expecting %ld bytes!\n",
           o_y0.size(),
           nParticles * nDOFs_ * sizeof(dfloat));

  time = t0;

  nParticles_ = nParticles;
  fieldOffset_ = computeFieldOffset(nParticles);

  handleAllocation(fieldOffset_);

  for (auto [fieldName, nFields] : interpFieldCounts) {
    laggedInterpFields[fieldName] =
        platform->device.malloc(nEXT * nFields * nrs->fieldOffset * sizeof(dfloat));
    extrapolatedInterpFields[fieldName] =
        platform->device.malloc(nFields * nrs->fieldOffset * sizeof(dfloat));
  }

  // set initial condition based on user-input
  if (nParticles_ > 0) {
    for (int dofId = 0; dofId < this->nDOFs(); ++dofId) {
      auto o_y_dof = getDOF(dofId);
      auto o_y0_dof = o_y0 + dofId * nParticles_ * sizeof(dfloat);
      o_y_dof.copyFrom(o_y0_dof, nParticles_ * sizeof(dfloat));
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
  for (int i = 0; i < order; ++i)
    coeffAB[i] *= dt[0];
  for (int i = order; i > order; i--)
    coeffAB[i - 1] = 0.0;
  o_coeffAB.copyFrom(coeffAB.data(), solverOrder * sizeof(dfloat));
}

void lpm_t::interpolate()
{
  for (auto [interpFieldName, interpFieldId] : interpFieldInputs) {
    interpolate(interpFieldName);
  }
}

void lpm_t::interpolate(std::string interpFieldName)
{
  if (timerLevel != TimerLevel::None) {
    platform->timer.tic(timerName + "integrate::userRHS::interpolate", 1);
  }
  auto o_fld = extrapolatedInterpFields.at(interpFieldName);
  auto o_interpFld = getInterpField(interpFieldName);
  const auto Nfields = numFieldsInterp(interpFieldName);

  interp->setTimerName(timerName + "integrate::userRHS::interpolate::");
  interp->eval(Nfields, nrs->fieldOffset, o_fld, fieldOffset_, o_interpFld);
  if (timerLevel != TimerLevel::None) {
    platform->timer.toc(timerName + "integrate::userRHS::interpolate");
  }
}

void lpm_t::integrate(dfloat tf)
{
  nrsCheck(!initialized_,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "%s\n",
           "cannot call integrate before calling initialize!");
  nrsCheck(!userRHS_,
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
      o_field.copyFrom(o_currentField, Nfields * nrs->fieldOffset * sizeof(dfloat));
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
    interp = std::make_unique<pointInterpolation_t>(nrs, bb_tol, newton_tol);
  }

  // set extrapolated state to t^n (copy from laggedInterpFields)
  for (auto [fieldName, o_field] : laggedInterpFields) {
    const auto Nfields = numFieldsInterp(fieldName);
    const auto Nbyte = (Nfields * sizeof(dfloat)) * nrs->fieldOffset;
    auto o_extField = extrapolatedInterpFields.at(fieldName);
    o_extField.copyFrom(o_field, Nfields * nrs->fieldOffset * sizeof(dfloat));
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
    userODESolver_(nrs, this, time, tf, tstep, o_y, userdata_, o_ydot);
  }
  else {
    if (solverType == SolverType::AB) {
      integrateAB();
    }
    else if (solverType == SolverType::RK) {
      integrateRK();
    }
  }

  dtEXT[1] = tf - time;
  time = tf;

  // lag previous time states in laggedInterpFields
  for (auto [fieldName, o_field] : laggedInterpFields) {
    const auto Nfields = numFieldsInterp(fieldName);
    const auto Nbyte = (Nfields * sizeof(dfloat)) * nrs->fieldOffset;
    for (int s = nEXT; s > 1; s--) {
      o_field.copyFrom(o_field, Nbyte, (s - 1) * Nbyte, (s - 2) * Nbyte);
    }

    auto o_currentField = interpFieldInputs.at(fieldName);

    // update most recent time state
    o_field.copyFrom(o_currentField, Nfields * nrs->fieldOffset * sizeof(dfloat));
  }

  // always provide (t^n,y^n) for next step
  this->find(this->o_y);

  if (timerLevel != TimerLevel::None) {
    platform->timer.toc(timerName + "integrate");
  }
}

// setup new particle coordinates
void lpm_t::find(occa::memory o_yNew)
{
  occa::memory o_xCoord, o_yCoord, o_zCoord;
  if (fieldOffset_) {
    o_xCoord = o_yNew + 0 * fieldOffset_ * sizeof(dfloat);
    o_yCoord = o_yNew + 1 * fieldOffset_ * sizeof(dfloat);
    o_zCoord = o_yNew + 2 * fieldOffset_ * sizeof(dfloat);
  }
  interp->setPoints(numParticles(), o_xCoord, o_yCoord, o_zCoord);
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

  auto mesh = nrs->meshV;

  std::copy(nrs->dt, nrs->dt + 3, dtEXT.begin());
  dtEXT[0] = tEXT - time;

  nek::extCoeff(coeffEXT.data(), dtEXT.data(), extOrder, bdfOrder);
  for (int i = nEXT; i > extOrder; i--) {
    coeffEXT[i - 1] = 0.0;
  }
  o_coeffEXT.copyFrom(coeffEXT.data(), nEXT * sizeof(dfloat));

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
      const auto Nbyte = (nDOFs_ * sizeof(dfloat)) * fieldOffset_;
      o_ydot.copyFrom(o_ydot, Nbyte, (s - 1) * Nbyte, (s - 2) * Nbyte);
    }
  }

  // boostrap using RK method
  if (tstep <= solverOrder) {
    const auto saveSolverOrder = solverOrder;
    this->solverOrder = bootstrapRKOrder;
    integrateRK4();
    this->solverOrder = saveSolverOrder;
    if (fieldOffset_ > 0) {
      o_ydot.copyFrom(o_k, nDOFs_ * fieldOffset_ * sizeof(dfloat)); // for later lagging
    }
  }
  else {
    platform->timer.tic(timerName + "integrate::userRHS", 1);
    userRHS_(nrs, this, time, o_y, userdata_, o_ydot);
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

  nrsCheck(false,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "RK solver order %d not supported!\n",
           solverOrder);
}

void lpm_t::integrateRK1()
{
  platform->timer.tic(timerName + "integrate::userRHS", 1);
  userRHS_(nrs, this, time, o_y, userdata_, o_k);
  platform->timer.toc(timerName + "integrate::userRHS");

  dfloat rkCoeff = dt[0];

  // o_y += dt[0] * o_ydot
  if (nParticles_)
    platform->linAlg->axpbyMany(nParticles_, nDOFs_, fieldOffset_, rkCoeff, o_k, 1.0, o_y);
}

void lpm_t::integrateRK2()
{
  occa::memory o_k1, o_k2;
  if(fieldOffset_ > 0){
    o_k1 = o_k + 0 * nDOFs_ * fieldOffset_ * sizeof(dfloat);
    o_k2 = o_k + 1 * nDOFs_ * fieldOffset_ * sizeof(dfloat);
  }

  platform->timer.tic(timerName + "integrate::userRHS", 1);
  userRHS_(nrs, this, time, o_y, userdata_, o_k1);
  platform->timer.toc(timerName + "integrate::userRHS");

  // o_ytmp = o_y + dt[0] * o_k1
  if (nParticles_)
    platform->linAlg->axpbyzMany(nParticles_, nDOFs_, fieldOffset_, 1.0, o_y, dt[0], o_k1, o_ytmp);

  this->find(o_ytmp);

  // extrapolate to t^{n+1} using t^{n}, t^{n-1}, ...
  extrapolateFluidState(time + dt[0]);

  platform->timer.tic(timerName + "integrate::userRHS", 1);
  userRHS_(nrs, this, time + dt[0], o_ytmp, userdata_, o_k2);
  platform->timer.toc(timerName + "integrate::userRHS");

  // o_y = o_y + 0.5 * dt[0] * (o_k1 + o_k2)
  coeffRK[0] = 0.5 * dt[0], coeffRK[1] = 0.5 * dt[0];
  o_coeffRK.copyFrom(coeffRK.data(), coeffRK.size() * sizeof(dfloat));

  if (nParticles_)
    this->nStagesSumManyKernel(nParticles_, fieldOffset_, solverOrder, nDOFs_, o_coeffRK, o_k, o_y);
}

void lpm_t::integrateRK3()
{
  occa::memory o_k1, o_k2, o_k3;
  if(fieldOffset_ > 0){
    o_k1 = o_k + 0 * nDOFs_ * fieldOffset_ * sizeof(dfloat);
    o_k2 = o_k + 1 * nDOFs_ * fieldOffset_ * sizeof(dfloat);
    o_k3 = o_k + 2 * nDOFs_ * fieldOffset_ * sizeof(dfloat);
  }

  platform->timer.tic(timerName + "integrate::userRHS", 1);
  userRHS_(nrs, this, time, o_y, userdata_, o_k1);
  platform->timer.toc(timerName + "integrate::userRHS");

  // o_ytmp = o_y + dt[0]/2 * o_k1
  if (nParticles_)
    platform->linAlg->axpbyzMany(nParticles_, nDOFs_, fieldOffset_, 1.0, o_y, 0.5 * dt[0], o_k1, o_ytmp);

  this->find(o_ytmp);

  // extrapolate to t^{n+1/2} using t^{n}, t^{n-1}, ...
  extrapolateFluidState(time + 0.5 * dt[0]);

  platform->timer.tic(timerName + "integrate::userRHS", 1);
  userRHS_(nrs, this, time + 0.5 * dt[0], o_y, userdata_, o_k2);
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
  userRHS_(nrs, this, time + dt[0], o_ytmp, userdata_, o_k3);
  platform->timer.toc(timerName + "integrate::userRHS");

  // o_y = o_y + dt[0] * (1/6 * o_k1 + 2/3 * o_k2 + 1/6 * o_k3)
  coeffRK[0] = 1.0 / 6.0 * dt[0], coeffRK[1] = 2.0 / 3.0 * dt[0], coeffRK[2] = 1.0 / 6.0 * dt[0];
  o_coeffRK.copyFrom(coeffRK.data(), coeffRK.size() * sizeof(dfloat));

  if (nParticles_)
    this->nStagesSumManyKernel(nParticles_, fieldOffset_, solverOrder, nDOFs_, o_coeffRK, o_k, o_y);
}

void lpm_t::integrateRK4()
{
  occa::memory o_k1, o_k2, o_k3, o_k4;
  if(fieldOffset_ > 0){
    o_k1 = o_k + 0 * nDOFs_ * fieldOffset_ * sizeof(dfloat);
    o_k2 = o_k + 1 * nDOFs_ * fieldOffset_ * sizeof(dfloat);
    o_k3 = o_k + 2 * nDOFs_ * fieldOffset_ * sizeof(dfloat);
    o_k4 = o_k + 3 * nDOFs_ * fieldOffset_ * sizeof(dfloat);
  }

  platform->timer.tic(timerName + "integrate::userRHS", 1);
  userRHS_(nrs, this, time, o_y, userdata_, o_k1);
  platform->timer.toc(timerName + "integrate::userRHS");

  // o_ytmp = o_y + dt[0]/2 * o_k1
  if (nParticles_)
    platform->linAlg->axpbyzMany(nParticles_, nDOFs_, fieldOffset_, 1.0, o_y, 0.5 * dt[0], o_k1, o_ytmp);

  this->find(o_ytmp);

  // extrapolate to t^{n+1/2} using t^{n}, t^{n-1}, ...
  extrapolateFluidState(time + 0.5 * dt[0]);

  platform->timer.tic(timerName + "integrate::userRHS", 1);
  userRHS_(nrs, this, time + 0.5 * dt[0], o_y, userdata_, o_k2);
  platform->timer.toc(timerName + "integrate::userRHS");

  // o_ytmp = o_y + dt[0]/2 * o_k2
  if (nParticles_)
    platform->linAlg->axpbyzMany(nParticles_, nDOFs_, fieldOffset_, 1.0, o_y, 0.5 * dt[0], o_k2, o_ytmp);

  this->find(o_ytmp);

  platform->timer.tic(timerName + "integrate::userRHS", 1);
  userRHS_(nrs, this, time + 0.5 * dt[0], o_y, userdata_, o_k3);
  platform->timer.toc(timerName + "integrate::userRHS");

  // o_ytmp = o_y + dt[0] * o_k3
  if (nParticles_)
    platform->linAlg->axpbyzMany(nParticles_, nDOFs_, fieldOffset_, 1.0, o_y, dt[0], o_k3, o_ytmp);

  this->find(o_ytmp);

  // extrapolate to t^{n+1} using t^{n}, t^{n-1}, ...
  extrapolateFluidState(time + dt[0]);

  platform->timer.tic(timerName + "integrate::userRHS", 1);
  userRHS_(nrs, this, time + dt[0], o_ytmp, userdata_, o_k4);
  platform->timer.toc(timerName + "integrate::userRHS");

  // o_y = o_y + dt[0] * (1/6 * o_k1 + 1/3 * o_k2 + 1/3 * o_k3 + 1/6 * o_k4)
  coeffRK[0] = 1.0 / 6.0 * dt[0];
  coeffRK[1] = 1.0 / 3.0 * dt[0];
  coeffRK[2] = 1.0 / 3.0 * dt[0];
  coeffRK[3] = 1.0 / 6.0 * dt[0];
  o_coeffRK.copyFrom(coeffRK.data(), coeffRK.size() * sizeof(dfloat));

  if (nParticles_)
    this->nStagesSumManyKernel(nParticles_, fieldOffset_, solverOrder, nDOFs_, o_coeffRK, o_k, o_y);
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

int lpm_t::fieldOffset(int n) { return computeFieldOffset(n); }

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

namespace {
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
  nrsCheck(entriesPerParticle > lpm_t::maxEntriesPerParticleMigration,
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
    interp->setPoints(this->numParticles(), o_xcoord, o_ycoord, o_zcoord);
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
    auto o_ySend = platform->device.malloc(nNonLocal * nDOFs_ * sizeof(dfloat));
    auto o_ydotSend = platform->device.malloc(solverOrder * nNonLocal * nDOFs_ * sizeof(dfloat));

    occa::memory o_propSend, o_interpFldSend;
    if (nProps_) {
      o_propSend = platform->device.malloc(nNonLocal * nProps_ * sizeof(dfloat));
    }
    if (nInterpFields_) {
      o_interpFldSend = platform->device.malloc(nNonLocal * nInterpFields_ * sizeof(dfloat));
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

    if (o_sendRankMap.size() < sendRankMap.size() * sizeof(dlong)) {
      if (o_sendRankMap.size())
        o_sendRankMap.free();
      o_sendRankMap = platform->device.malloc(sendRankMap.size() * sizeof(dlong));
    }

    o_sendRankMap.copyFrom(sendRankMap.data(), sendRankMap.size() * sizeof(dlong));

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
    o_ySend.copyTo(ySend.data(), nNonLocal * nDOFs_ * sizeof(dfloat));
    o_ydotSend.copyTo(ydotSend.data(), solverOrder * nNonLocal * nDOFs_ * sizeof(dfloat));
    if (nProps_) {
      o_propSend.copyTo(propSend.data(), nNonLocal * nProps_ * sizeof(dfloat));
    }
    if (nInterpFields_) {
      o_interpFldSend.copyTo(interpFldSend.data(), nNonLocal * nInterpFields_ * sizeof(dfloat));
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
  dlong newNParticles = this->numParticles() - nUnfound - nNonLocal + nReceived;
  dlong newFieldOffset = computeFieldOffset(newNParticles);

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

      if (o_migrateMap.size() < migrateMap.size() * sizeof(dlong)) {
        if (o_migrateMap.size())
          o_migrateMap.free();
        o_migrateMap = platform->device.malloc(migrateMap.size() * sizeof(dlong));
      }
      if (this->numParticles() > 0) {
        o_migrateMap.copyFrom(migrateMap.data(), migrateMap.size() * sizeof(dlong));
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
    auto o_yRecv = platform->device.malloc(nReceived * nDOFs_ * sizeof(dfloat), yRecv.data());
    auto o_ydotRecv =
        platform->device.malloc(solverOrder * nReceived * nDOFs_ * sizeof(dfloat), ydotRecv.data());

    occa::memory o_propRecv, o_interpFldRecv;
    if (nProps_) {
      o_propRecv = platform->device.malloc(nReceived * nProps_ * sizeof(dfloat), propRecv.data());
    }
    if (nInterpFields_) {
      o_interpFldRecv =
          platform->device.malloc(nReceived * nInterpFields_ * sizeof(dfloat), interpFldRecv.data());
    }

    std::vector<dlong> recvRankMap(nReceived, -1);
    for (int pid = 0; pid < nReceived; ++pid) {
      recvRankMap[pid] = ctr;
      ctr++;
    }

    if (o_recvRankMap.size() < recvRankMap.size() * sizeof(dlong)) {
      if (o_recvRankMap.size())
        o_recvRankMap.free();
      o_recvRankMap = platform->device.malloc(recvRankMap.size() * sizeof(dlong));
    }

    o_recvRankMap.copyFrom(recvRankMap.data(), recvRankMap.size() * sizeof(dlong));

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

  if (o_yOld.size())
    o_yOld.free();
  if (o_ydotOld.size())
    o_ydotOld.free();
  if (o_propOld.size())
    o_propOld.free();
  if (o_interpFldOld.size())
    o_interpFldOld.free();

  // do an additional findpts call
  if (timerLevel != TimerLevel::None) {
    platform->timer.tic(timerName + "migrate::find", 1);
  }

  {
    auto o_xcoord = getDOF("x");
    auto o_ycoord = getDOF("y");
    auto o_zcoord = getDOF("z");
    interp->setPoints(this->numParticles(), o_xcoord, o_ycoord, o_zcoord);
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
  nrsCheck(yNewPart.size() < expectedYSize,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "yNewPart size is %ld but expected %d words!\n",
           yNewPart.size(),
           expectedYSize);
  nrsCheck(propNewPart.size() < expectedPropSize,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "propNewPart size is %ld but expected %d words!\n",
           propNewPart.size(),
           expectedPropSize);
  nrsCheck(ydotNewPart.size() < expectedYdotSize,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "ydotNewPart size is %ld but expected %d words!\n",
           ydotNewPart.size(),
           expectedYdotSize);

  auto o_yNewPart = platform->device.malloc(expectedYSize * sizeof(dfloat), yNewPart.data());
  auto o_propNewPart = platform->device.malloc(expectedPropSize * sizeof(dfloat), propNewPart.data());
  auto o_ydotNewPart = platform->device.malloc(expectedYdotSize * sizeof(dfloat));

  addParticles(newNParticles, o_yNewPart, o_propNewPart, o_ydotNewPart);

  o_yNewPart.free();
  o_propNewPart.free();
  o_ydotNewPart.free();
}

void lpm_t::addParticles(int newNParticles, occa::memory o_yNewPart, occa::memory o_propNewPart)
{
  const auto expectedYdotSize = solverOrder * newNParticles * nDOFs_;
  auto o_ydotNewPart = platform->device.malloc(expectedYdotSize * sizeof(dfloat));
  addParticles(newNParticles, o_yNewPart, o_propNewPart, o_ydotNewPart);
}

void lpm_t::addParticles(int newNParticles,
                         occa::memory o_yNewPart,
                         occa::memory o_propNewPart,
                         occa::memory o_ydotNewPart)
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
  int newOffset = computeFieldOffset(this->nParticles_ + newNParticles);

  // check that the sizes of o_yNewPart, o_propNewPart are correct
  auto expectedYSize = incomingOffset * nDOFs_ * sizeof(dfloat);
  auto expectedPropSize = incomingOffset * nProps_ * sizeof(dfloat);
  auto expectedYdotSize = solverOrder * incomingOffset * nDOFs_ * sizeof(dfloat);
  nrsCheck(o_yNewPart.size() < expectedYSize,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "o_yNewPart size is %" PRId64 " but expected %ld bytes!\n",
           o_yNewPart.size(),
           expectedYSize);
  nrsCheck(o_propNewPart.size() < expectedPropSize,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "o_propNewPart size is %" PRId64 " but expected %ld bytes!\n",
           o_propNewPart.size(),
           expectedPropSize);
  nrsCheck(o_ydotNewPart.size() < expectedYdotSize,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "o_ydotNewPart size is %" PRId64 " but expected %ld bytes!\n",
           o_ydotNewPart.size(),
           expectedYdotSize);

  std::vector<dlong> remainingMap(this->nParticles_, 0);
  std::vector<dlong> insertMap(newNParticles, 0);

  // remainingMap[id] = id for existing particles
  if (this->nParticles_) {
    std::iota(remainingMap.begin(), remainingMap.end(), 0);
    if (o_remainingMap.size() < this->nParticles_ * sizeof(dlong)) {
      if (o_remainingMap.size())
        o_remainingMap.free();
      o_remainingMap = platform->device.malloc(this->nParticles_ * sizeof(dlong));
    }
    o_remainingMap.copyFrom(remainingMap.data(), this->nParticles_ * sizeof(dlong));
  }

  // insertMap[id] = id + nParticles_ for incoming particles
  if (newNParticles) {
    std::iota(insertMap.begin(), insertMap.end(), this->nParticles_);
    if (o_insertMap.size() < newNParticles * sizeof(dlong)) {
      if (o_insertMap.size())
        o_insertMap.free();
      o_insertMap = platform->device.malloc(newNParticles * sizeof(dlong));
    }
    o_insertMap.copyFrom(insertMap.data(), newNParticles * sizeof(dlong));
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
      o_interpFldDummy = platform->device.malloc(incomingOffset * nInterpFields_ * sizeof(dfloat));
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

  if (o_propOld.size())
    o_propOld.free();
  if (o_interpFldOld.size())
    o_interpFldOld.free();
  if (o_yOld.size())
    o_yOld.free();
  if (o_ydotOld.size())
    o_ydotOld.free();

  nParticles_ += newNParticles;
  fieldOffset_ = newOffset;

  // update findpts results for newly added particles
  {
    auto o_xcoord = getDOF("x");
    auto o_ycoord = getDOF("y");
    auto o_zcoord = getDOF("z");
    interp->setPoints(this->numParticles(), o_xcoord, o_ycoord, o_zcoord);

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

  const auto Nbytes = this->numParticles() * sizeof(dlong);
  if (Nbytes > 0) {
    if (o_remainingMap.size() < Nbytes) {
      if (o_remainingMap.size())
        o_remainingMap.free();
      o_remainingMap = platform->device.malloc(Nbytes);
    }
    o_remainingMap.copyFrom(remainingMap.data(), Nbytes);
  }

  auto o_propOld = this->o_prop;
  auto o_interpFldOld = this->o_interpFld;
  auto o_yOld = this->o_y;
  auto o_ydotOld = this->o_ydot;

  const auto newNParticles = this->numParticles() - nDelete;
  const auto newOffset = computeFieldOffset(newNParticles);

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

  if (o_propOld.size())
    o_propOld.free();
  if (o_interpFldOld.size())
    o_interpFldOld.free();
  if (o_yOld.size())
    o_yOld.free();
  if (o_ydotOld.size())
    o_ydotOld.free();

  nParticles_ = newNParticles;
  fieldOffset_ = newOffset;

  // update findpts results
  {
    auto o_xcoord = getDOF("x");
    auto o_ycoord = getDOF("y");
    auto o_zcoord = getDOF("z");
    interp->setPoints(this->numParticles(), o_xcoord, o_ycoord, o_zcoord);

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

namespace {
std::string lpm_vtu_data(std::string fieldName, long long int nComponent, long long int distance)
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
    interp->setPoints(this->numParticles(), o_xcoord, o_ycoord, o_zcoord);

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
  int mpi_rank = platform->comm.mpiRank;
  int mpi_size = platform->comm.mpiCommSize;

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

  std::string message = "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
  message += "\t<UnstructuredGrid>\n";
  message += "\t\t<FieldData>\n";
  message += "\t\t\t<DataArray type=\"Float32\" Name=\"TIME\" NumberOfTuples=\"1\" format=\"ascii\"> " +
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
    if (!isOutput)
      continue;
    const auto Nfields = propCounts.at(propName);
    message += "\t\t\t\t" + lpm_vtu_data(propName, Nfields, offset);
    offset += (Nfields * globalNPartOutput) * sizeof(float) + 1 * sizeof(long long int);
  }

  // output interpolated fields
  for (auto [interpFieldName, isOutput] : outputInterpFields) {
    if (!isOutput)
      continue;
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
    if (!isOutput)
      continue;
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
    if (!isOutput)
      continue;
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

  kernelName = "remapParticles";
  fileName = oklpath + "/plugins/" + kernelName + ".okl";
  platform->kernels.add(kernelName, fileName, kernelInfo);
}

void lpm_t::setTimerLevel(TimerLevel level)
{
  timerLevel = level;
  interp->setTimerLevel(level);
}

TimerLevel lpm_t::getTimerLevel() const { return timerLevel; }

void lpm_t::setTimerName(std::string name)
{
  timerName = name;
  interp->setTimerName(name);
}

namespace {
long long int parseNumParticles(std::string restartfile, const std::string &header)
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
  }
  catch (std::invalid_argument e) {
    errorLogger << "Could not read number of particles while reading " << restartfile << "!\n";
    errorLogger << "Exception said:\n" << e.what() << std::endl;
  }

  auto errorString = errorLogger.str();
  int errorLength = errorString.length();
  MPI_Allreduce(MPI_IN_PLACE, &errorLength, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);
  
  nrsCheck(errorLength > 0,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "Error in parseNumParticles:\n%s",
           errorString.c_str());

  return nparticles;
}

auto parsePointData(std::string restartfile, std::string pointData)
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
    errorLogger << "Could not parse " << fieldName << " number of components while reading " << restartfile << "\n";
  }

  try {
    numComponents = std::stoll(match[1].str());
  }
  catch (std::invalid_argument &e) {
    errorLogger << "Could not parse " << fieldName << " number of components while reading " << restartfile << "!\n";
    errorLogger << "Exception said:\n" << e.what() << std::endl;
  }

  found = std::regex_search(pointData, match, std::regex(R"(\s*offset=\"(\d+)\")"));
  if (!found) {
    errorLogger << "Could not parse " << fieldName << " offset while reading " << restartfile << "!\n";
  }

  try {
    offset = std::stoll(match[1].str());
  }
  catch (std::invalid_argument &e) {
    errorLogger << "Could not parse " << fieldName << " offset while reading " << restartfile << "!\n";
    errorLogger << "Exception said:\n" << e.what() << std::endl;
  }
  
  auto errorString = errorLogger.str();
  int errorLength = errorString.length();
  MPI_Allreduce(MPI_IN_PLACE, &errorLength, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);
  
  nrsCheck(errorLength > 0,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "Error in parsePointData:\n%s",
           errorString.c_str());

  return std::make_tuple(fieldName, numComponents, offset);
}

auto readHeader(std::string restartFile)
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

void lpm_t::restart(std::string restartFile)
{
  bool fileExists = std::filesystem::exists(restartFile);
  nrsCheck(!fileExists, platform->comm.mpiComm, EXIT_FAILURE, "Restart file %s does not exist!\n", restartFile.c_str());

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
    dfloat t0;
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

  nrsCheck(nPointData != nPartGlobal,
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

  o_xCoord.copyFrom(xCoord.data(), nPartLocal * sizeof(dfloat));
  o_yCoord.copyFrom(yCoord.data(), nPartLocal * sizeof(dfloat));
  o_zCoord.copyFrom(zCoord.data(), nPartLocal * sizeof(dfloat));

  auto readField = [&, &header = header](std::string fieldName, long long int expectedNumComponents, long long int offset) {
    nrsCheck(fieldType.count(fieldName) == 0,
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
    }
    else if (type == FieldType::PROP) {
      nComponents = propCounts.at(fieldName);
      o_fld = getProp(fieldName);
    }
    else if (type == FieldType::INTERP_FIELD) {
      nComponents = interpFieldCounts.at(fieldName);
      o_fld = getInterpField(fieldName);
    }

    nrsCheck(nComponents != expectedNumComponents,
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

    nrsCheck(nPointData != nPartGlobal,
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

    o_fld.copyFrom(fldHost.data(), this->fieldOffset() * nComponents * sizeof(dfloat));
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
    }
    else {
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

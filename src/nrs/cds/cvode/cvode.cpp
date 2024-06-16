#include "platform.hpp"
#include "cvode.hpp"
#include "elliptic.hpp"
#include "maskedFaceIds.hpp"
#include "inipp.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"

#include "plugins/lowMach.hpp"
#include "bdry.hpp"
#include "tabularPrinter.hpp"

#ifdef ENABLE_CVODE

// cvode includes
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sundials/sundials_types.h"
#include "sundials/sundials_math.h"
#include "cvode/cvode.h"

#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#include "nvector/nvector_serial.h"
#include "nvector/nvector_mpiplusx.h"

#ifdef ENABLE_CUDA
#include "nvector/nvector_cuda.h"
#endif
#ifdef ENABLE_HIP
#include "nvector/nvector_hip.h"
#endif

#include "getN_VectorMemory.hpp"

#endif

namespace
{

void check_retval(void *returnvalue, const char *funcname, int opt)
{
}

} // namespace

cvode_t::cvode_t(cds_t *_cds)
{
}

void cvode_t::initialize()
{
}

cvode_t::~cvode_t()
{
}

void cvode_t::setupDirichletMask()
{
}

void cvode_t::applyDirichlet(double time)
{
}

void cvode_t::rhs(double time, const LVector_t<dfloat> &o_y, LVector_t<dfloat> &o_ydot)
{
}

int cvode_t::jtvRhs(double time, const occa::memory &o_y, occa::memory &o_ydot)
{
  return 0;
}

void cvode_t::jtvRhs(double time, const LVector_t<dfloat> &o_y, LVector_t<dfloat> &o_ydot)
{
}

void cvode_t::defaultRHS(double time, double t0, const LVector_t<dfloat> &o_y, LVector_t<dfloat> &o_ydot)
{
}

void cvode_t::makeq(double time)
{
}

void cvode_t::nrsToCv(occa::memory o_EField, LVector_t<dfloat> &o_LField, bool isYdot)
{
}

void cvode_t::cvToNrs(const LVector_t<dfloat> &o_LField, occa::memory o_EField, bool isYdot)
{
}

void cvode_t::solve(double t0, double t1, int tstep)
{
}

long cvode_t::numSteps() const
{
  return this->nsteps;
}

long cvode_t::numRHSEvals() const
{
  return this->nrhs;
}

long cvode_t::numNonlinSolveIters() const
{
  return this->nni;
}

long cvode_t::numLinIters() const
{
  return this->nli;
}

void cvode_t::updateCounters() {}

void cvode_t::resetCounters() {}

void cvode_t::printInfo(bool printVerboseInfo)
{
}

void cvode_t::printTimers()
{
}

void cvode_t::resetTimers()
{
}

std::string cvode_t::rhsTagName() const
{
  return 0;
}

void cvode_t::setLocalPointSource(userLocalPointSourceE_t _userLocalPointSource)
{
}

void cvode_t::setLocalPointSource(userLocalPointSourceL_t _userLocalPointSource)
{
}

void cvode_t::computeUrst()
{
}

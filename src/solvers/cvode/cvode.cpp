#include "cvode.hpp"
#include "elliptic.h"
#include "inipp.hpp"
#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "Urst.hpp"
#include <cstdlib>
#include <limits>
#include <array>
#include <numeric>
#include "nrssys.hpp"
#include "ogs.hpp"
#include "udf.hpp"

#include "timeStepper.hpp"
#include "plugins/lowMach.hpp"
#include "bdry.hpp"
#include "tabularPrinter.hpp"

#ifdef ENABLE_CVODE
// cvode includes
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <nvector/nvector_mpiplusx.h>
#ifdef ENABLE_CUDA
#include <nvector/nvector_cuda.h>
#endif
#ifdef ENABLE_HIP
#include <nvector/nvector_hip.h>
#endif
#endif

//#define USE_E_VECTOR_LAYOUT 1

namespace {

#ifdef ENABLE_CVODE
sunrealtype *__N_VGetDeviceArrayPointer(N_Vector u)
{
}
#endif

void check_retval(void *returnvalue, const char *funcname, int opt)
{
}

} // namespace

cvode_t::cvode_t(nrs_t *_nrs)
{
}

cvode_t::~cvode_t()
{
}

void cvode_t::initialize()
{
}

void cvode_t::setupEToLMapping()
{
}

void cvode_t::setupDirichletMask()
{
}

void cvode_t::applyDirichlet(dfloat time)
{
}

void cvode_t::computeErrorWeight(occa::memory o_y, occa::memory o_ewt)
{
}

void cvode_t::rhs(dfloat time, occa::memory o_y, occa::memory o_ydot)
{
}

void cvode_t::jtvRHS(dfloat time, occa::memory o_y, occa::memory o_ydot)
{
}

void cvode_t::defaultRHS(dfloat time, dfloat t0, occa::memory o_y, occa::memory o_ydot)
{
}

void cvode_t::makeq(dfloat time)
{
}

void cvode_t::nrsToCv(occa::memory o_EField, occa::memory o_LField)
{
}

void cvode_t::cvToNrs(occa::memory o_LField, occa::memory o_EField)
{
}

void cvode_t::solve(double t0, double t1, int tstep)
{
}

#ifdef ENABLE_CVODE
long cvode_t::numSteps() const
{
}

long cvode_t::numRHSEvals() const
{
}

long cvode_t::numNonlinSolveIters() const
{
}

long cvode_t::numLinIters() const
{
}
#else
long cvode_t::numSteps() const { return 0; }

long cvode_t::numRHSEvals() const { return 0; }

long cvode_t::numNonlinSolveIters() const { return 0; }

long cvode_t::numLinIters() const { return 0; }
#endif

void cvode_t::printInfo(bool printVerboseInfo) const
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
}


void cvode_t::setLocalPointSource(userLocalPointSource_t _userLocalPointSource)
{
}

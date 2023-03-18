#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"

namespace lowMach
{
void setup(nrs_t* nrs, dfloat alpha_, occa::memory o_beta_, occa::memory o_kappa_);
void buildKernel(occa::properties kernelInfo);
void qThermalSingleComponent(dfloat time, occa::memory o_div);
void dpdt(occa::memory o_FU);
}

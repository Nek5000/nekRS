#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"

namespace lowMach
{
void setup(nrs_t* nrs, dfloat alpha);
void buildKernel(occa::properties kernelInfo);
void qThermalSingleComponent(dfloat time, occa::memory o_div,occa::memory o_beta, occa::memory o_kappa);
void dpdt(occa::memory o_FU);
}

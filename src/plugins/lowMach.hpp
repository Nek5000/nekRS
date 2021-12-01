#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"

namespace lowMach
{
void setup(nrs_t* nrs);
void buildKernel(occa::properties kernelInfo);
void qThermalIdealGasSingleComponent(dfloat time, occa::memory o_div,dfloat gamma);
void qThermalRealGasSingleComponent(dfloat time, occa::memory o_div,occa::memory o_beta,occa::memory o_kappa);
void dpdt(occa::memory o_FU);
}

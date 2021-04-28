#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"

namespace lowMach
{
void setup(nrs_t* nrs, dfloat gamma0);
void buildKernel(nrs_t* nrs);
void qThermalIdealGasSingleComponent(dfloat time, occa::memory o_div);
void dpdt(occa::memory o_FU);
}

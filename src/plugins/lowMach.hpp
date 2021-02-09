#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"

namespace lowMach
{
void setup(nrs_t* nrs);
void qThermalIdealGasSingleComponent(dfloat time, dfloat gamma, occa::memory o_div);
void dpdt(dfloat gamma, occa::memory o_FU);
}

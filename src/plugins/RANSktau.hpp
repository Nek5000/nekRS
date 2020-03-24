#include <nekrs.hpp>
#include <nekInterfaceAdapter.hpp>

namespace RANSktau {

void buildKernel(ins_t *ins);
void updateSourceTerms();
void setup(ins_t *insIn, dfloat mue, dfloat rho, int startIndex);
void setup(ins_t *insIn, dfloat mue, dfloat rho, int startIndex, const dfloat *coeffIn);
void updateProperties();
occa::memory o_mue_t();

}

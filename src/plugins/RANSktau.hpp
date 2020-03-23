#include <nekrs.hpp>
#include <nekInterfaceAdapter.hpp>

namespace RANSktau {

void buildKernel(ins_t *ins);
void sourceTerms();
void setup(ins_t *insIn, dfloat mue, dfloat rho, int startIndex);
void setup(ins_t *insIn, dfloat mue, dfloat rho, int startIndex, const dfloat *coeffIn);
void mue(dfloat C, occa::memory o_mue);
void mue(dfloat mueLam, dfloat C, occa::memory o_mue);
dfloat sigma_k();
dfloat sigma_tau();

}

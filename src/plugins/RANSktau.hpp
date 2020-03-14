#include <nekrs.hpp>
#include <nekInterfaceAdapter.hpp>

namespace RANSktau {

void buildKernel(ins_t *ins);
void sourceTerms(occa::memory o_BF, occa::memory o_BFDiag);
void setup(ins_t *insIn, occa::memory o_kIn, occa::memory o_tauIn);
void setup(ins_t *insIn, occa::memory o_kIn, occa::memory o_tauIn, dfloat *coeffIn);
void mue(dfloat C, occa::memory o_mue);
void mue(dfloat mueLam, dfloat C, occa::memory o_mue);
dfloat sigma_k(void);
dfloat sigma_tau(void);

}

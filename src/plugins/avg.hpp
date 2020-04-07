#include <nekrs.hpp>
#include <nekInterfaceAdapter.hpp>

namespace avg {

void buildKernel(ins_t *ins);
void run(dfloat time);
void setup(ins_t *ins_);
void outfld();
void reset();
void EX (dlong N, dfloat a, dfloat b, int nflds, occa::memory o_x, occa::memory o_EX);
void EXX(dlong N, dfloat a, dfloat b, int nflds, occa::memory o_x, occa::memory o_EXX);
void EXY(dlong N, dfloat a, dfloat b, int nflds, occa::memory o_x, occa::memory o_y, occa::memory o_EXY);
}

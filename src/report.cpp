#include "nekrs.hpp"
#include "nekInterfaceAdapter.hpp"

void report(ins_t *ins, dfloat time, int tstep){
  nek_ocopyFrom(ins, time, tstep);
  nek_outfld();
}

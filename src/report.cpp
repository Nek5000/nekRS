#include "nekrs.hpp"
#include "nekInterfaceAdapter.hpp"

void report(ins_t *ins, dfloat time, int tstep){
  // copy data back to host
  ins->o_U.copyTo(ins->U);
  ins->o_P.copyTo(ins->P);

  nek_copyFrom(ins, time, tstep);
  nek_userchk();
  nek_outfld();
}

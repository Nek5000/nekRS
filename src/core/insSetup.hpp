#if !defined(nekrs_inssetup_hpp_)
#define nekrs_inssetup_hpp_

#include "nrs.hpp"
ins_t* insSetup(MPI_Comm comm, occa::device device, setupAide &options, int buildOnly);

#endif

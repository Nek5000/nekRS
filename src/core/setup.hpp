#if !defined(nekrs_inssetup_hpp_)
#define nekrs_inssetup_hpp_

#include "nrs.hpp"
nrs_t* nrsSetup(MPI_Comm comm, occa::device device, setupAide &options, int buildOnly);

#endif

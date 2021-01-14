#if !defined(nekrs_inssetup_hpp_)
#define nekrs_inssetup_hpp_

#include "nrs.hpp"
void nrsSetup(MPI_Comm comm, occa::device device, setupAide &options, int buildOnly, nrs_t *nrs);

#endif

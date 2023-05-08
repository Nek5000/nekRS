#if !defined(nekrs_parreader_hpp_)
#define nekrs_parreader_hpp_

#include "nrs.hpp"
#include "inipp.hpp"

void parRead(inipp::Ini *par, const std::string& setupFile, MPI_Comm comm, setupAide &options);
void parEcho();

#endif

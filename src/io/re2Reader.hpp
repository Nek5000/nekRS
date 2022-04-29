#if !defined(nekrs_re2reader_hpp_)
#define nekrs_re2reader_hpp_

#include "nrs.hpp"

namespace re2 
{
void nelg(const std::string& meshFile, int& nelgt, int& nelgv, MPI_Comm comm);
}

#endif

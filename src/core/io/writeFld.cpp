#include "platform.hpp"
#include "nekInterfaceAdapter.hpp"
#include "fld.hpp"

namespace fld
{
ElementFilter elementFilter;
}

void fld::write(std::string suffix, double t, int step,
                const std::vector<occa::memory>& o_u, const occa::memory& o_p, const std::vector<occa::memory>& o_s,
                bool outXYZ, bool FP64, int Nout, bool uniform)
{
  nek::outfld(suffix.c_str(), t, step, outXYZ, FP64, o_u, o_p, o_s, Nout, uniform); 
}

void fld::write(std::string suffix, double t, int step,
                const std::vector<occa::memory>& o_s, 
                bool outXYZ, bool FP64, int Nout, bool uniform)
{
  write(suffix, t, step, std::vector<occa::memory>(), o_NULL, o_s, outXYZ, FP64, Nout, uniform); 
}



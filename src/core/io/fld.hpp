#if !defined(nekrs_io_hpp_)
#define nekrs_io_hpp_

#include "platform.hpp"

class ElementFilter
{
  public:
    // mask: local element indices to include 
    void set(const std::vector<int>& mask) 
    { 
      _mask = mask;
    };
    void clear() { _mask.clear(); };

    std::vector<int> mask() const { return _mask; };

  private:
    std::vector<int> _mask;
};

namespace fld {

extern ElementFilter elementFilter;

void write(std::string suffix, double t, int step,
           const std::vector<occa::memory>& o_s, 
           bool outXYZ = true, bool FP64 = false, int Nout = 0, bool uniform = false);

void write(std::string suffix, double t, int step,
           const std::vector<occa::memory>& o_u, const occa::memory& o_p, const std::vector<occa::memory>& o_s,
           bool outXYZ = true, bool FP64 = false, int Nout = 0, bool uniform = false);

}

#endif

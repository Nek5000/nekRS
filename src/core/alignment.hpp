#if !defined(alignment_hpp_)
#define alignment_hpp_
#include "nrssys.hpp"
#include <array>
#include <vector>
enum class boundaryAlignment_t { X, Y, Z, UNALIGNED };
std::string to_string(boundaryAlignment_t a);
boundaryAlignment_t computeAlignment(const std::array<dfloat, 3> &n);
#endif
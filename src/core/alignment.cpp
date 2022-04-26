#include "alignment.hpp"
std::string to_string(boundaryAlignment_t a)
{
  switch (a) {
  case boundaryAlignment_t::X:
    return "X";
  case boundaryAlignment_t::Y:
    return "Y";
  case boundaryAlignment_t::Z:
    return "Z";
  case boundaryAlignment_t::UNALIGNED:
    return "UNALIGNED";
  }

  return "";
}
boundaryAlignment_t computeAlignment(const std::array<dfloat, 3> &n)
{
  const dfloat alignmentTol = 1e-4;
  const dfloat nxDiff = std::abs(std::abs(n[0]) - 1.0);
  const dfloat nyDiff = std::abs(std::abs(n[1]) - 1.0);
  const dfloat nzDiff = std::abs(std::abs(n[2]) - 1.0);

  if (nxDiff < alignmentTol)
    return boundaryAlignment_t::X;
  if (nyDiff < alignmentTol)
    return boundaryAlignment_t::Y;
  if (nzDiff < alignmentTol)
    return boundaryAlignment_t::Z;

  return boundaryAlignment_t::UNALIGNED;
}
#include "alignment.hpp"
std::string to_string(alignment_t a)
{
  switch (a) {
  case alignment_t::X:
    return "X";
  case alignment_t::Y:
    return "Y";
  case alignment_t::Z:
    return "Z";
  case alignment_t::UNALIGNED:
    return "UNALIGNED";
  }

  return "";
}
alignment_t computeAlignment(const std::array<dfloat,3>& n)
{
  const dfloat alignmentTol = 1e-4;
  const dfloat nxDiff = std::abs(std::abs(n[0]) - 1.0);
  const dfloat nyDiff = std::abs(std::abs(n[1]) - 1.0);
  const dfloat nzDiff = std::abs(std::abs(n[2]) - 1.0);

  if (nxDiff < alignmentTol)
    return alignment_t::X;
  if (nyDiff < alignmentTol)
    return alignment_t::Y;
  if (nzDiff < alignmentTol)
    return alignment_t::Z;

  return alignment_t::UNALIGNED;
}
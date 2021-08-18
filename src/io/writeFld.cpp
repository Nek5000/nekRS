#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"

void writeFld(const char* suffix, dfloat t, int outXYZ, int FP64,
              void* o_u, void* o_p, void* o_s,
              int NSfields)
{
  nek::outfld(suffix, t, outXYZ, FP64, o_u, o_p, o_s, NSfields); 
}

void writeFld(nrs_t *nrs, dfloat t, int outXYZ, int FP64) 
{
  int Nscalar = 0;
  occa::memory o_s;
  if(nrs->Nscalar) {
    o_s = nrs->cds->o_S;
    Nscalar = nrs->Nscalar;
  }
  nek::outfld("   ", t, outXYZ, FP64, &nrs->o_U, &nrs->o_P, &o_s, Nscalar); 
}

void writeFld(nrs_t *nrs, dfloat t) 
{
  string precision;
  platform->options.getArgs("CHECKPOINT PRECISION", precision);
  int FP64 = 0;
  if(precision == "DP") FP64 = 1;

  int outXYZ = 1;
  if(platform->options.compareArgs("CHECKPOINT OUTPUT MESH", "FALSE")) outXYZ = 0;

  writeFld(nrs, t, outXYZ, FP64); 
}

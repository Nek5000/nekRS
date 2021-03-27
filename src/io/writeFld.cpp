#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"

void writeFld(const char* suffix, dfloat t, int coords, int FP64,
              void* o_u, void* o_p, void* o_s,
              int NSfields)
{
  nek::outfld(suffix, t, coords, FP64, o_u, o_p, o_s, NSfields); 
}

void writeFld(nrs_t *nrs, dfloat t, int FP64) 
{
  int coords = 1;
  int Nscalar = 0;
  occa::memory o_s;
  if(nrs->Nscalar) {
    o_s = nrs->cds->o_S;
    Nscalar = nrs->Nscalar;
  }
  nek::outfld("   ", t, coords, FP64, &nrs->o_U, &nrs->o_P, &o_s, Nscalar); 
}

void writeFld(nrs_t *nrs, dfloat t) 
{
  writeFld(nrs, t, 0); 
}

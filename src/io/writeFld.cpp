#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"

void writeFld(std::string suffix, dfloat t, int outXYZ, int FP64,
              void* o_u, void* o_p, void* o_s,
              int NSfields)
{
  std::string casename;
  platform->options.getArgs("CASENAME", casename);

  std::string filename; 
  if(suffix.length())
    filename = suffix + casename; 
  else
    filename = casename; 

  nek::outfld(filename.c_str(), t, outXYZ, FP64, o_u, o_p, o_s, NSfields); 
}

void writeFld(nrs_t *nrs, dfloat t, int outXYZ, int FP64, std::string suffix) 
{
  int Nscalar = 0;
  occa::memory o_s;
  if(nrs->Nscalar) {
    o_s = nrs->cds->o_S;
    Nscalar = nrs->Nscalar;
  }
  writeFld(suffix, t, outXYZ, FP64, &nrs->o_U, &nrs->o_P, &o_s, Nscalar); 
}

void writeFld(nrs_t *nrs, dfloat t, int outXYZ, int FP64) 
{
  writeFld(nrs, t, outXYZ, FP64, "");
}

void writeFld(nrs_t *nrs, dfloat t, std::string suffix) 
{
  std::string precision;
  platform->options.getArgs("CHECKPOINT PRECISION", precision);
  int FP64 = 0;
  if(precision == "DP") FP64 = 1;

  int outXYZ = 1;
  if(platform->options.compareArgs("CHECKPOINT OUTPUT MESH", "FALSE")) outXYZ = 0;

  writeFld(nrs, t, outXYZ, FP64, suffix); 
}

void writeFld(nrs_t *nrs, dfloat t) 
{
  writeFld(nrs, t, "");
}

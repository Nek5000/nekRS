#include "nrs.hpp"

void writeFld(nrs_t *nrs, dfloat t);
void writeFld(nrs_t *nrs, dfloat t, int FP64);
void writeFld(const char* suffix, dfloat t, int coords, int FP64,
              occa::memory &o_u, occa::memory &o_p, occa::memory &o_s,
              int NSfields);

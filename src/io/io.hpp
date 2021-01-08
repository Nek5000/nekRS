#include "nrs.hpp"

void writeFld(nrs_t *nrs, dfloat t);
void writeFld(nrs_t *nrs, dfloat t, int FP64);
void writeFld(const char* suffix, dfloat t, int coords, int FP64,
              void* o_u, void *o_p,  void *o_s,
              int NSfields);

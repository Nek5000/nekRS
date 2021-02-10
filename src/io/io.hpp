#include "nrs.hpp"

void copyFile(const char *srcName, const char* destName);
bool isFileNewer(const char *file1, const char* file2);
void writeFld(nrs_t *nrs, dfloat t);
void writeFld(nrs_t *nrs, dfloat t, int FP64);
void writeFld(const char* suffix, dfloat t, int coords, int FP64,
              void* o_u, void *o_p,  void *o_s,
              int NSfields);

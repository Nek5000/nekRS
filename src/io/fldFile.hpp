#if !defined(nekrs_io_hpp_)
#define nekrs_io_hpp_

#include "nrs.hpp"

void fileSync(const char *file);
void copyFile(const char *srcName, const char* destName);
bool isFileEmpty(const char *file);
bool isFileNewer(const char *file1, const char* file2);
bool fileExists(const char *file);
void writeFld(nrs_t *nrs, dfloat t, int step);
void writeFld(nrs_t *nrs, dfloat t, int step, std::string suffix);
void writeFld(nrs_t *nrs, dfloat t, int step, int outXYZ, int FP64);
void writeFld(nrs_t *nrs, dfloat t, int step, int outXYZ, int FP64, std::string suffix);

void writeFld(std::string suffix, dfloat t, int step, int outXYZ, int FP64,
              void* o_s, int NSfields);

void writeFld(std::string suffix, dfloat t, int step, int outXYZ, int FP64,
              void* o_u, void *o_p,  void *o_s,
              int NSfields);

#endif

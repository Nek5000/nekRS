#if !defined(nekrs_nrs_hpp_)
#define nekrs_nrs_hpp_

#include <string>
#include <mpi.h>

namespace nekrs
{
void setup(MPI_Comm comm, int buildOnly, int sizeTarget,
           int ciMode, std::string cacheDir, std::string setupFile,
           std::string backend, std::string deviceID);

void runStep(double time, double dt, int tstep);
void copyToNek(double time, int tstep);
void udfExecuteStep(double time, int tstep, int isOutputStep);
void nekOutfld(void);
void nekUserchk(void);
void printRuntimeStatistics(void);

const double dt(void);
const int outputStep(void);
const int NtimeSteps(void);
const double startTime(void);
const double finalTime(void);

void* nekPtr(const char* id);
}

#endif

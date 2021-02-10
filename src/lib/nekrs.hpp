#if !defined(nekrs_nrs_hpp_)
#define nekrs_nrs_hpp_

#include <string>
#include <mpi.h>

namespace nekrs
{
void setup(MPI_Comm comm, int buildOnly, int targetSize,
           int ciMode, std::string cacheDir, std::string setupFile,
           std::string backend, std::string deviceID);

void runStep(double time, double dt, int tstep);
void copyToNek(double time, int tstep);
void udfExecuteStep(double time, int tstep, int isOutputStep);
void outfld(double time);
const int isOutputStep(double time, int tStep);
void nekUserchk(void);
void printRuntimeStatistics(void);
const double writeInterval(void);
const double dt(void);
const double startTime(void);
const double endTime(void);
const int numSteps(void);
const int lastStep(double time, int tstep, double elapsedTime);
const int writeControlRunTime(void);

void* nrsPtr(void);
void* nekPtr(const char* id);
}

#endif

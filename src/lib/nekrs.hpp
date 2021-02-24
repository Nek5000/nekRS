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
int isOutputStep(double time, int tStep);
void nekUserchk(void);
void printRuntimeStatistics(void);
double writeInterval(void);
double dt(void);
double startTime(void);
double endTime(void);
int numSteps(void);
int lastStep(double time, int tstep, double elapsedTime);
int writeControlRunTime(void);

void* nrsPtr(void);
void* nekPtr(const char* id);
}

#endif

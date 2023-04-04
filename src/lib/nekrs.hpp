#if !defined(nekrs_nrs_hpp_)
#define nekrs_nrs_hpp_

#include <mpi.h>
#include <functional> 
#include <string>

namespace nekrs
{
void setup(MPI_Comm commg_in, MPI_Comm comm_in,
           int buildOnly, int commSizeTarget,
           int ciMode, std::string _setupFile,
           std::string _backend, std::string _deviceID,
           int _nSessions, int _sessionID,
           int debug);
void copyFromNek(double time, int tstep);
void udfExecuteStep(double time, int tstep, int isOutputStep);
void outfld(double time, int step);
void outfld(double time, int step, std::string suffix);
int outputStep(double time, int tStep);
void outputStep(int val);
int finalize();
void nekUserchk(void);
int runTimeStatFreq();
int printInfoFreq();
int updateFileCheckFreq();
void printRuntimeStatistics(int step);
double writeInterval(void);
double dt(int tStep);
double startTime(void);
double endTime(void);
int numSteps(void);
void lastStep(int val);
int lastStep(double time, int tstep, double elapsedTime);
int writeControlRunTime(void);
int exitValue(void);
bool stepConverged(void);
void processUpdFile();
void printInfo(double time, int tstep, bool printStepInfo, bool printVerboseInfo);
void verboseInfo(bool enabled);
void updateTimer(const std::string &key, double time);
void resetTimer(const std::string &key); 
void* nrsPtr(void);
void* nekPtr(const char* id);
void initStep(double time, double dt, int tstep);
bool runStep(std::function<bool(int)> convergenceCheck, int corrector);
bool runStep(int corrector);
double finishStep();
bool stepConverged();

}

#endif

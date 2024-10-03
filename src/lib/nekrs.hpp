#if !defined(nekrs_nrs_hpp_)
#define nekrs_nrs_hpp_

#include <mpi.h>
#include <functional>
#include <string>

/*
  Basic high-level API
*/

namespace nekrs
{
void setup(MPI_Comm commg_in,
           MPI_Comm comm_in,
           int buildOnly,
           int commSizeTarget,
           int ciMode,
           const std::map<std::string, std::map<std::string, std::string>> &parKeyValuePairs,
           std::string casename,
           std::string _backend,
           std::string _deviceID,
           int nSessions,
           int sessionID,
           int debug);
void udfExecuteStep(double time, int tstep, int checkpointStep);
void writeCheckpoint(double time, int step);
int checkpointStep(double time, int tStep);
void checkpointStep(int val);
int finalize();
void nekUserchk(void);
int runTimeStatFreq();
int printStepInfoFreq();
int updateFileCheckFreq();
void printRuntimeStatistics(int step);
double writeInterval(void);
std::tuple<double, double> dt(int tStep);
double startTime(void);
double endTime(void);
int numSteps(void);
void lastStep(int val);
int lastStep(double time, int tstep, double elapsedTime);
int writeControlRunTime(void);
int exitValue(void);
bool stepConverged(void);
void processUpdFile();
void printStepInfo(double time, int tstep, bool printStepInfo, bool printVerboseInfo);
void updateTimer(const std::string &key, double time);
void resetTimer(const std::string &key);

void initStep(double time, double dt, int tstep);
bool runStep(std::function<bool(int)> convergenceCheck, int corrector);
bool runStep(int corrector);
void finishStep();
bool stepConverged();
int timeStep();
double finalTimeStepSize(double time);
} // namespace nekrs

#endif

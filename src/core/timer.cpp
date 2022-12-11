#include <iostream>
#include <string>
#include <map>
#include <algorithm>

#include "timer.hpp"
#include "platform.hpp"
#include "ogs.hpp"

namespace timer {
namespace {
struct tagData {
  long long int count;
  double hostElapsed;
  double deviceElapsed;
  double startTime;
  occa::streamTag startTag;
};
std::map<std::string, tagData> m_;

const int NEKRS_TIMER_INVALID_KEY = -1;
const int NEKRS_TIMER_INVALID_METRIC = -2;

int ifSync_;
inline int ifSync() { return ifSync_; }

int enable_sync_;

int enabled;

occa::device device_;
MPI_Comm comm_;

inline void sync()
{
  if (enable_sync_)
    MPI_Barrier(comm_);
}

double tElapsedTime = 0;
} // namespace

timer_t::timer_t(MPI_Comm comm, occa::device device, int ifSyncDefault, int enableSync)
{
  init(comm, device, ifSyncDefault, enableSync);
}

void timer_t::init(MPI_Comm comm, occa::device device, int ifSyncDefault, int enableSync)
{
  device_ = device;
  ifSync_ = ifSyncDefault;
  comm_ = comm;
  enable_sync_ = enableSync;
  enabled = 1;
}

void timer_t::set(const std::string tag, double time, long long int count)
{
  m_[tag].startTime = time;
  auto it = m_.find(tag);
  if (it == m_.end()) {
    printf("Error in set: Invalid tag name %s\n", tag.c_str());
    MPI_Abort(comm_, 1);
  }

  it->second.hostElapsed = time;
  it->second.deviceElapsed = it->second.hostElapsed;
  it->second.count = count;
}

void timer_t::enable() { enabled = 1; }

void timer_t::disable() { enabled = 0; }

void timer_t::reset()
{
  for (auto &it : m_) {
    it.second.startTime = 0;
    it.second.hostElapsed = 0;
    it.second.deviceElapsed = 0;
    it.second.count = 0;
  }
  ogsResetTime();
}

void timer_t::enableSync() { enable_sync_ = 1; }

void timer_t::disableSync() { enable_sync_ = 0; }

void timer_t::reset(const std::string tag)
{
  std::map<std::string, tagData>::iterator it = m_.find(tag);
  it->second.startTime = 0;
  it->second.hostElapsed = 0;
  it->second.deviceElapsed = 0;
  it->second.count = 0;
}

void timer_t::finalize() { reset(); }

void timer_t::deviceTic(const std::string tag, int ifSync)
{
  if (!enabled)
    return;
  if (ifSync)
    sync();
  m_[tag].startTag = device_.tagStream();
}

void timer_t::deviceTic(const std::string tag)
{
  if (!enabled)
    return;
  if (ifSync())
    sync();
  m_[tag].startTag = device_.tagStream();
}

void timer_t::deviceToc(const std::string tag)
{
  if (!enabled)
    return;
  occa::streamTag stopTag = device_.tagStream();

  std::map<std::string, tagData>::iterator it = m_.find(tag);
  if (it == m_.end()) {
    printf("Error in deviceToc: Invalid tag name %s\n", tag.c_str());
    MPI_Abort(comm_, 1);
  }

  it->second.deviceElapsed += device_.timeBetween(it->second.startTag, stopTag);
  it->second.count++;
}

void timer_t::hostTic(const std::string tag, int ifSync)
{
  if (!enabled)
    return;
  if (ifSync)
    sync();
  m_[tag].startTime = MPI_Wtime();
}

void timer_t::hostTic(const std::string tag)
{
  if (!enabled)
    return;
  if (ifSync())
    sync();
  m_[tag].startTime = MPI_Wtime();
}

void timer_t::hostToc(const std::string tag)
{
  if (!enabled)
    return;
  double stopTime = MPI_Wtime();

  auto it = m_.find(tag);
  if (it == m_.end()) {
    printf("Error in deviceToc: Invalid tag name %s\n", tag.c_str());
    MPI_Abort(comm_, 1);
  }

  it->second.hostElapsed += (stopTime - it->second.startTime);
  it->second.count++;
}

void timer_t::tic(const std::string tag, int ifSync)
{
  if (!enabled)
    return;
  if (ifSync)
    sync();
  m_[tag].startTime = MPI_Wtime();
  m_[tag].startTag = device_.tagStream();
}

void timer_t::tic(const std::string tag)
{
  if (!enabled)
    return;
  if (ifSync())
    sync();
  m_[tag].startTime = MPI_Wtime();
  m_[tag].startTag = device_.tagStream();
}

void timer_t::toc(const std::string tag)
{
  if (!enabled)
    return;
  auto stopTime = MPI_Wtime();
  auto stopTag = device_.tagStream();

  auto it = m_.find(tag);
  if (it == m_.end()) {
    printf("Error in deviceToc: Invalid tag name %s\n", tag.c_str());
    MPI_Abort(comm_, 1);
  }

  it->second.hostElapsed += (stopTime - it->second.startTime);
  it->second.deviceElapsed += device_.timeBetween(it->second.startTag, stopTag);
  it->second.count++;
}

double timer_t::hostElapsed(const std::string tag)
{
  auto it = m_.find(tag);
  if (it == m_.end())
    return NEKRS_TIMER_INVALID_KEY;
  return it->second.hostElapsed;
}

double timer_t::deviceElapsed(const std::string tag)
{
  auto it = m_.find(tag);
  if (it == m_.end())
    return NEKRS_TIMER_INVALID_KEY;
  return it->second.deviceElapsed;
}

long long int timer_t::count(const std::string tag)
{
  auto it = m_.find(tag);
  if (it == m_.end())
    return NEKRS_TIMER_INVALID_KEY;
  return it->second.count;
}

double timer_t::query(const std::string tag, const std::string metric)
{
  int size;
  MPI_Comm_size(comm_, &size);

  auto it = m_.find(tag);
  if (it == m_.end())
    return NEKRS_TIMER_INVALID_KEY;
  auto hostElapsed = it->second.hostElapsed;
  auto deviceElapsed = it->second.deviceElapsed;
  auto count = it->second.count;

  double retVal;

  std::string upperMetric = metric;
  std::transform(upperMetric.begin(), upperMetric.end(), upperMetric.begin(), ::toupper);

  if (upperMetric.compare("HOST:MIN") == 0) {
    MPI_Allreduce(&hostElapsed, &retVal, 1, MPI_DOUBLE, MPI_MIN, comm_);
    return retVal;
  }
  if (upperMetric.compare("HOST:MAX") == 0) {
    MPI_Allreduce(&hostElapsed, &retVal, 1, MPI_DOUBLE, MPI_MAX, comm_);
    return retVal;
  }
  if (upperMetric.compare("HOST:SUM") == 0) {
    MPI_Allreduce(&hostElapsed, &retVal, 1, MPI_DOUBLE, MPI_SUM, comm_);
    return retVal;
  }
  if (upperMetric.compare("HOST:AVG") == 0) {
    MPI_Allreduce(&hostElapsed, &retVal, 1, MPI_DOUBLE, MPI_SUM, comm_);
    return retVal / (size * count);
  }
  if (upperMetric.compare("DEVICE:MIN") == 0) {
    MPI_Allreduce(&deviceElapsed, &retVal, 1, MPI_DOUBLE, MPI_MIN, comm_);
    return retVal;
  }
  if (upperMetric.compare("DEVICE:MAX") == 0) {
    MPI_Allreduce(&deviceElapsed, &retVal, 1, MPI_DOUBLE, MPI_MAX, comm_);
    return retVal;
  }
  if (upperMetric.compare("DEVICE:SUM") == 0) {
    MPI_Allreduce(&deviceElapsed, &retVal, 1, MPI_DOUBLE, MPI_SUM, comm_);
    return retVal;
  }
  if (upperMetric.compare("DEVICE:AVG") == 0) {
    MPI_Allreduce(&deviceElapsed, &retVal, 1, MPI_DOUBLE, MPI_SUM, comm_);
    return retVal / (size * count);
  }
  return NEKRS_TIMER_INVALID_METRIC;
}

std::string printPercentage(double num, double dom)
{
  char buf[4096];
  double frac = num / dom;
  snprintf(buf, sizeof(buf), "%4.1f", 100 * frac);
  return std::string(buf);
}

void timer_t::printStatEntry(std::string name, std::string tag, std::string type, double tNorm)
{
  int rank;
  MPI_Comm_rank(comm_, &rank);
  const long long int nCalls = count(tag);
  const double tTag = query(tag, type);
  const bool child = (tNorm != tElapsedTime);
  if (tTag > 0) {
    if (rank == 0) {
      std::cout << name << tTag << "s"
                << "  " << printPercentage(tTag, tElapsedTime);
      if (child)
        std::cout << "  " << printPercentage(tTag, tNorm);
      else
        std::cout << "      ";
      std::cout << "  " << nCalls << "\n";
    }
  }
}

void timer_t::printStatEntry(std::string name, double time, double tNorm)
{
  int rank;
  MPI_Comm_rank(comm_, &rank);
  const bool child = (tNorm != tElapsedTime);
  if (time > 0) {
    if (rank == 0) {
      std::cout << name << time << "s"
                << "  " << printPercentage(time, tElapsedTime);
      if (child)
        std::cout << "  " << printPercentage(time, tNorm);
      else
        std::cout << "      ";
      std::cout << "\n";
    }
  }
}

void timer_t::printRunStat(int step)
{
  int rank;
  MPI_Comm_rank(comm_, &rank);

  set("velocity proj",
      query("velocity proj pre", "DEVICE:MAX") + query("velocity proj post", "DEVICE:MAX"),
      count("velocity proj pre"));

  set("pressure proj",
      query("pressure proj pre", "DEVICE:MAX") + query("pressure proj post", "DEVICE:MAX"),
      count("pressure proj pre"));

  set("scalar proj",
      query("scalar proj pre", "DEVICE:MAX") + query("scalar proj post", "DEVICE:MAX"),
      count("scalar proj pre"));

  double gsTime = ogsTime(/* reportHostTime */ true);
  MPI_Allreduce(MPI_IN_PLACE, &gsTime, 1, MPI_DOUBLE, MPI_MAX, comm_);

  tElapsedTime = query("elapsed", "DEVICE:MAX");

  if (rank == 0)
    std::cout << "\n>>> runtime statistics (step= " << step << "  totalElapsed= " << tElapsedTime << "s"
              << "):\n";

  std::cout.setf(std::ios::scientific);
  int outPrecisionSave = std::cout.precision();
  std::cout.precision(5);

  if (rank == 0)
    std::cout << "name                    "
              << "time          "
              << "abs%  "
              << "rel%  "
              << "calls\n";

  const double tElapsedTimeSolve = query("elapsedStepSum", "DEVICE:MAX");
  const double tSetup = query("setup", "DEVICE:MAX");

  const double tMinSolveStep = query("minSolveStep", "DEVICE:MAX");
  const double tMaxSolveStep = query("maxSolveStep", "DEVICE:MAX");
  const double flops =
      platform->flopCounter->get(platform->comm.mpiComm) / (tElapsedTimeSolve * platform->comm.mpiCommSize);
  bool printFlops = !platform->options.compareArgs("PRESSURE PRECONDITIONER", "SEMFEM");

  printStatEntry("  setup                 ", "setup", "DEVICE:MAX", tElapsedTime);
  printStatEntry("    loadKernels         ", "loadKernels", "HOST:MAX", tSetup);

  printStatEntry("  solve                 ", tElapsedTimeSolve, tElapsedTime);
  if (tElapsedTimeSolve > 0 && rank == 0) {
    std::cout << "    min                 " << tMinSolveStep << "s\n";
    std::cout << "    max                 " << tMaxSolveStep << "s\n";
    if (printFlops)
      std::cout << "    flops/rank          " << flops << "\n";
  }

  printStatEntry("    checkpointing       ", "checkpointing", "DEVICE:MAX", tElapsedTimeSolve);
  printStatEntry("    udfExecuteStep      ", "udfExecuteStep", "DEVICE:MAX", tElapsedTimeSolve);

  const double tMakef = query("makef", "DEVICE:MAX");
  printStatEntry("    makef               ", "makef", "DEVICE:MAX", tElapsedTimeSolve);
  printStatEntry("      udfUEqnSource     ", "udfUEqnSource", "DEVICE:MAX", tMakef);

  const double tMakeq = query("makeq", "DEVICE:MAX");
  printStatEntry("    makeq               ", "makeq", "DEVICE:MAX", tElapsedTime);
  printStatEntry("      udfSEqnSource     ", "udfSEqnSource", "DEVICE:MAX", tMakeq);

  printStatEntry("    udfProperties       ", "udfProperties", "DEVICE:MAX", tElapsedTimeSolve);

  printStatEntry("    meshUpdate          ", "meshUpdate", "DEVICE:MAX", tElapsedTimeSolve);
  const double tMesh = query("meshSolve", "DEVICE:MAX");
  printStatEntry("    meshSolve           ", "meshSolve", "DEVICE:MAX", tElapsedTimeSolve);
  printStatEntry("      initial guess     ", "mesh proj", "DEVICE:MAX", tMesh);

  const double tVelocity = query("velocitySolve", "DEVICE:MAX");
  printStatEntry("    velocitySolve       ", "velocitySolve", "DEVICE:MAX", tElapsedTimeSolve);
  printStatEntry("      rhs               ", "velocity rhs", "DEVICE:MAX", tVelocity);
  printStatEntry("      preconditioner    ", "velocity preconditioner", "DEVICE:MAX", tVelocity);
  printStatEntry("      initial guess     ", "velocity proj", "DEVICE:MAX", tVelocity);

  const double tPressure = query("pressureSolve", "DEVICE:MAX");
  printStatEntry("    pressureSolve       ", "pressureSolve", "DEVICE:MAX", tElapsedTimeSolve);
  printStatEntry("      rhs               ", "pressure rhs", "DEVICE:MAX", tPressure);

  const double tPressurePreco = query("pressure preconditioner", "DEVICE:MAX");
  printStatEntry("      preconditioner    ", "pressure preconditioner", "DEVICE:MAX", tPressure);

  for (int i = 15; i > 0; i--) {
    const std::string tag = "pressure preconditioner smoother N=" + std::to_string(i);
    if (m_.find(tag) == m_.end())
      continue;
    printStatEntry("        pMG smoother    ", tag, "DEVICE:MAX", tPressurePreco);
  }

  printStatEntry("        coarse grid     ", "coarseSolve", "DEVICE:MAX", tPressurePreco);
  printStatEntry("      initial guess     ", "pressure proj", "DEVICE:MAX", tPressure);

  int nScalar = 0;
  platform->options.getArgs("NUMBER OF SCALARS", nScalar);

  const double tScalar = query("scalarSolve", "DEVICE:MAX");
  printStatEntry("    scalarSolve         ", "scalarSolve", "DEVICE:MAX", tElapsedTimeSolve);
  printStatEntry("      rhs               ", "scalar rhs", "DEVICE:MAX", tScalar);

  auto precoTimeScalars = 0.0;
  auto precoCallsScalars = 0.0;
  for (int is = 0; is < nScalar; is++) {
    std::stringstream ss;
    const int scalarWidth = getDigitsRepresentation(NSCALAR_MAX - 1);
    ss << std::setfill('0') << std::setw(scalarWidth) << is;
    std::string sid = ss.str();
    precoTimeScalars += query("scalar" + sid + " preconditioner", "DEVICE:MAX"); 
    precoCallsScalars += count("scalar" + sid + " preconditioner");
  }  
  set("scalar preconditioner", precoTimeScalars, precoCallsScalars);

  printStatEntry("      preconditioner    ", "scalar preconditioner", "DEVICE:MAX", tScalar);
  printStatEntry("      initial guess     ", "scalar proj", "DEVICE:MAX", tScalar);


  printStatEntry("    gsMPI               ", gsTime, tElapsedTimeSolve);

  printStatEntry("    dotp                ", "dotp", "DEVICE:MAX", tElapsedTimeSolve);

  printStatEntry("    dotp multi          ", "dotpMulti", "DEVICE:MAX", tElapsedTimeSolve);

  if (rank == 0)
    std::cout << std::endl;

  std::cout.unsetf(std::ios::scientific);
  std::cout.precision(outPrecisionSave);
}

} // namespace timer

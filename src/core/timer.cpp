#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include <tuple>

#include "timer.hpp"
#include "platform.hpp"
#include "ogs.hpp"
#include "orderedMap.hpp"

namespace timer {
namespace {
struct tagData {
  long long int count;
  double hostElapsed;
  double deviceElapsed;
  double startTime;
  occa::streamTag startTag;
};
orderedMap<std::string, tagData> m_;

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

double tElapsedTimeSolve = 0;

auto sumAllMatchingTags(std::function<bool(std::string)> predicate, const std::string metric)
{
  long long int count = 0;
  double elapsed = 0;
  const auto timerTags = platform->timer.tags();

  // filter out tags that do contain tag in the name
  std::vector<std::string> filteredTags;
  std::copy_if(timerTags.begin(),
               timerTags.end(),
               std::back_inserter(filteredTags),
               [&](const std::string &tag) { return predicate(tag); });

  for (auto &&t : filteredTags) {
    count += platform->timer.count(t);
    elapsed += platform->timer.query(t, metric);
  }

  return std::make_tuple(elapsed, count);
}

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

void timer_t::clear()
{
  m_.clear();
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

void timer_t::printStatEntry(std::string name, double tTag, long long int nCalls, double tNorm)
{
  int rank;
  MPI_Comm_rank(comm_, &rank);
  const bool child = (tNorm != tElapsedTimeSolve);
  if (tTag > 0) {
    if (rank == 0) {
      std::cout << name << tTag << "s"
                << "  " << printPercentage(tTag, tElapsedTimeSolve);
      if (child)
        std::cout << "  " << printPercentage(tTag, tNorm);
      else
        std::cout << "      ";
      std::cout << "  " << nCalls << "\n";
    }
  }
}

void timer_t::printStatEntry(std::string name, std::string tag, std::string type, double tNorm)
{
  int rank;
  MPI_Comm_rank(comm_, &rank);
  const long long int nCalls = count(tag);
  const double tTag = query(tag, type);
  const bool child = (tNorm != tElapsedTimeSolve);
  if (tTag > 0) {
    if (rank == 0) {
      std::cout << name << tTag << "s"
                << "  " << printPercentage(tTag, tElapsedTimeSolve);
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
  const bool child = (tNorm != tElapsedTimeSolve);
  if (time > 0) {
    if (rank == 0) {
      std::cout << name << time << "s"
                << "  " << printPercentage(time, tElapsedTimeSolve);
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

  set("mesh proj",
      query("mesh proj pre", "DEVICE:MAX") + query("mesh proj post", "DEVICE:MAX"),
      count("mesh proj pre"));

  double gsTime = ogsTime(/* reportHostTime */ true);
  MPI_Allreduce(MPI_IN_PLACE, &gsTime, 1, MPI_DOUBLE, MPI_MAX, comm_);

  const double tElapsedTime = query("elapsed", "DEVICE:MAX");

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

  tElapsedTimeSolve = query("elapsedStepSum", "DEVICE:MAX");
  const double tSetup = query("setup", "DEVICE:MAX");

  const double tMinSolveStep = query("minSolveStep", "DEVICE:MAX");
  const double tMaxSolveStep = query("maxSolveStep", "DEVICE:MAX");
  const double flops =
      platform->flopCounter->get(platform->comm.mpiComm) / (tElapsedTimeSolve * platform->comm.mpiCommSize);
  bool printFlops = !platform->options.compareArgs("PRESSURE PRECONDITIONER", "SEMFEM");

  printStatEntry("  solve                 ", tElapsedTimeSolve, tElapsedTimeSolve);
  if (tElapsedTimeSolve > 0 && rank == 0) {
    std::cout << "    min                 " << tMinSolveStep << "s\n";
    std::cout << "    max                 " << tMaxSolveStep << "s\n";
    if (printFlops)
      std::cout << "    flops/rank          " << flops << "\n";
  }

  auto lpmLocalKernelPredicate = [](const std::string &tag) {
    return tag.find("lpm_t::") != std::string::npos && tag.find("localKernel") != std::string::npos;
  };

  auto lpmLocalEvalKernelPredicate = [](const std::string &tag) {
    return tag.find("lpm_t::") != std::string::npos && tag.find("localEvalKernel") != std::string::npos;
  };

  auto neknekLocalKernelPredicate = [](const std::string &tag) {
    return tag.find("neknek_t::") != std::string::npos && tag.find("localKernel") != std::string::npos;
  };

  auto neknekLocalEvalKernelPredicate = [](const std::string &tag) {
    return tag.find("neknek_t::") != std::string::npos && tag.find("localEvalKernel") != std::string::npos;
  };

  printStatEntry("    checkpointing       ", "checkpointing", "DEVICE:MAX", tElapsedTimeSolve);
  printStatEntry("    udfExecuteStep      ", "udfExecuteStep", "DEVICE:MAX", tElapsedTimeSolve);
  const double tudf = query("udfExecuteStep", "DEVICE:MAX");
  printStatEntry("      lpm integrate     ", "lpm_t::integrate", "DEVICE:MAX", tudf);
  const double tlpm = query("lpm_t::integrate", "DEVICE:MAX");
  printStatEntry("        userRHS         ", "lpm_t::integrate::userRHS", "DEVICE:MAX", tlpm);
  const double tParticleRHS = query("lpm_t::integrate::userRHS", "DEVICE:MAX");
  printStatEntry("          interpolate   ",
                 "lpm_t::integrate::userRHS::interpolate",
                 "DEVICE:MAX",
                 tParticleRHS);
  const double tInterpPart = query("lpm_t::integrate::userRHS::interpolate", "DEVICE:MAX");
  auto [tLocalKernel, nLocalKernel] = sumAllMatchingTags(lpmLocalEvalKernelPredicate, "DEVICE:MAX");
  printStatEntry("            eval kernel ", tLocalKernel, nLocalKernel, tInterpPart);
  printStatEntry("        findpts         ", "lpm_t::integrate::find", "DEVICE:MAX", tlpm);
  const double tFindPart = query("lpm_t::integrate::find", "DEVICE:MAX");
  auto [tFindKernel, nFindKernel] = sumAllMatchingTags(lpmLocalKernelPredicate, "DEVICE:MAX");
  printStatEntry("          find kernel   ", tFindKernel, nFindKernel, tFindPart);
  printStatEntry("        delete          ", "lpm_t::deleteParticles", "DEVICE:MAX", tlpm);
  printStatEntry("      lpm add           ", "lpm_t::addParticles", "DEVICE:MAX", tudf);
  printStatEntry("      lpm write         ", "lpm_t::write", "DEVICE:MAX", tudf);
  
  const double tDiv = query("udfDiv", "DEVICE:MAX");
  printStatEntry("    udfDiv              ", "udfDiv", "DEVICE:MAX", tElapsedTimeSolve);
  printStatEntry("      udfSEqnSource     ", "udfDiv::udfSEqnSource", "DEVICE:MAX", tDiv);

  const double tMakef = query("makef", "DEVICE:MAX");
  printStatEntry("    makef               ", "makef", "DEVICE:MAX", tElapsedTimeSolve);
  printStatEntry("      udfUEqnSource     ", "udfUEqnSource", "DEVICE:MAX", tMakef);

  const double tMakeq = query("makeq", "DEVICE:MAX");
  printStatEntry("    makeq               ", "makeq", "DEVICE:MAX", tElapsedTimeSolve);
  printStatEntry("      udfSEqnSource     ", "udfSEqnSource", "DEVICE:MAX", tMakeq);

  printStatEntry("    udfProperties       ", "udfProperties", "DEVICE:MAX", tElapsedTimeSolve);

  printStatEntry("    meshUpdate          ", "meshUpdate", "DEVICE:MAX", tElapsedTimeSolve);
  const double tMesh = query("meshSolve", "DEVICE:MAX");
  printStatEntry("    meshSolve           ", "meshSolve", "DEVICE:MAX", tElapsedTimeSolve);
  printStatEntry("      preconditioner    ", "mesh preconditioner", "DEVICE:MAX", tMesh);
  printStatEntry("      initial guess     ", "mesh proj", "DEVICE:MAX", tMesh);

  const double tNekNek = query("neknek update boundary", "DEVICE:MAX");
  printStatEntry("    neknek              ", "neknek update boundary", "DEVICE:MAX", tElapsedTimeSolve);
  printStatEntry("      sync              ", "neknek sync", "DEVICE:MAX", tNekNek);
  printStatEntry("      exchange          ", "neknek exchange", "DEVICE:MAX", tNekNek);
  const double tExchange = query("neknek exchange", "DEVICE:MAX");
  std::tie(tLocalKernel, nLocalKernel) = sumAllMatchingTags(neknekLocalEvalKernelPredicate, "DEVICE:MAX");
  printStatEntry("        eval kernel     ", tLocalKernel, nLocalKernel, tExchange);
  printStatEntry("      findpts           ", "neknek updateInterpPoints", "DEVICE:MAX", tNekNek);
  const double tFindpts = query("neknek updateInterpPoints", "DEVICE:MAX");

  if (tFindpts > 0.0) {
    std::tie(tFindKernel, nFindKernel) = sumAllMatchingTags(neknekLocalKernelPredicate, "DEVICE:MAX");
    printStatEntry("        find kernel     ", tFindKernel, nFindKernel, tFindpts);
  }

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
  
  const double tScalarCvode = query("cvode_t::solve", "DEVICE:MAX");

  auto cvodeMakeQPredicate = [](const std::string &tag) {
    bool match = tag.find("cvode_t::") != std::string::npos && tag.find("makeq") != std::string::npos;
    // ensure children of the timer aren't doubly counted
    return match && tag.find("makeq::") == std::string::npos;
  };
  auto [tMakeqCvode, nMakeqCvode] = sumAllMatchingTags(cvodeMakeQPredicate, "DEVICE:MAX");
  
  auto cvodeUdfSEqnSourcePredicate = [](const std::string &tag) {
    bool match = tag.find("cvode_t::") != std::string::npos && tag.find("udfSEqnSource") != std::string::npos;
    // ensure children of the timer aren't doubly counted
    return match && tag.find("udfSEqnSource::") == std::string::npos;
  };
  auto [tSEqnSourceCvode, nSEqnSourceCvode] = sumAllMatchingTags(cvodeUdfSEqnSourcePredicate, "DEVICE:MAX");
  
  auto cvodePropertiesPredicate = [](const std::string &tag) {
    bool match = tag.find("cvode_t::") != std::string::npos && tag.find("evaluateProperties") != std::string::npos;
    // ensure children of the timer aren't doubly counted
    return match && tag.find("evaluateProperties::") == std::string::npos;
  };
  auto [tPropCvode, nPropCvode] = sumAllMatchingTags(cvodePropertiesPredicate, "DEVICE:MAX");
  printStatEntry("    scalarSolveCvode    ", "cvode_t::solve", "DEVICE:MAX", tElapsedTimeSolve);
  printStatEntry("      makeq             ", tMakeqCvode, nMakeqCvode, tScalarCvode);
  printStatEntry("        udfSEqnSource   ", tSEqnSourceCvode, nSEqnSourceCvode, tMakeqCvode);
  printStatEntry("      udfProperties     ", tPropCvode, nPropCvode, tScalarCvode);

  auto precoTimeScalars = 0.0;
  auto precoCallsScalars = 0.0;
  for (int is = 0; is < nScalar; is++) {
    std::string sid = scalarDigitStr(is);
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

void timer_t::printAll()
{
  if (platform->comm.mpiRank != 0)
    return;
  std::cout << "Device timers: {\n";
  for (auto &&[name, data] : m_) {
    std::cout << "\t" << name << " " << data.deviceElapsed << ",\n";
  }
  std::cout << "}\n";

  std::cout << "Host timers: {\n";
  for (auto &&[name, data] : m_) {
    std::cout << "\t" << name << " " << data.hostElapsed << ",\n";
  }
  std::cout << "}\n";
}

std::vector<std::string> timer_t::tags()
{
  std::vector<std::string> entries;
  auto &keys = m_.keys();
  std::copy(keys.begin(), keys.end(), std::back_inserter(entries));
  return entries;
}

} // namespace timer

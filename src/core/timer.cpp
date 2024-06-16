#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include <tuple>

#include "platform.hpp"
#include "ogs.hpp"
#include "orderedMap.hpp"

#include "tabularPrinter.hpp"

namespace timer
{
namespace
{
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

std::vector<std::string> userStat;

int ifSync_;

inline int ifSync()
{
  return ifSync_;
}

int enable_sync_;

int enabled;

occa::device device_;
MPI_Comm comm_;

inline void sync()
{
  if (enable_sync_) {
    MPI_Barrier(comm_);
  }
}

double tElapsedTimeSolve = 0;

} // namespace

void timer_t::printStatSetElapsedTimeSolve(double time)
{
  tElapsedTimeSolve = time;
}

std::tuple<double, long long int> timer_t::sumAllMatchingTags(std::function<bool(std::string)> predicate, const std::string metric)
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

  nekrsCheck(it == m_.end(), MPI_COMM_SELF, EXIT_FAILURE, "Invalid tag name %s\n", tag.c_str());

  it->second.hostElapsed = time;
  it->second.deviceElapsed = it->second.hostElapsed;
  it->second.count = count;
}

void timer_t::enable()
{
  enabled = 1;
}

void timer_t::disable()
{
  enabled = 0;
}

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

void timer_t::enableSync()
{
  enable_sync_ = 1;
}

void timer_t::disableSync()
{
  enable_sync_ = 0;
}

void timer_t::reset(const std::string timerName)
{
  const auto timerTags = platform->timer.tags();
  std::vector<std::string> filteredTags;
  std::copy_if(timerTags.begin(),
               timerTags.end(),
               std::back_inserter(filteredTags),
               [&](const std::string &tag) { return tag.find(timerName) == 0; });

  for (auto &&tag : filteredTags) {
    std::map<std::string, tagData>::iterator it = m_.find(tag);
    if (it != m_.end()) {
      it->second.startTime = 0;
      it->second.hostElapsed = 0;
      it->second.deviceElapsed = 0;
      it->second.count = 0;
    }
  }
}

void timer_t::finalize()
{
  reset();
}

void timer_t::deviceTic(const std::string tag, int ifSync)
{
  if (!enabled) {
    return;
  }
  if (ifSync) {
    sync();
  }
  m_[tag].startTag = device_.tagStream();
}

void timer_t::deviceToc(const std::string tag)
{
  if (!enabled) {
    return;
  }
  occa::streamTag stopTag = device_.tagStream();

  std::map<std::string, tagData>::iterator it = m_.find(tag);
  nekrsCheck(it == m_.end(), MPI_COMM_SELF, EXIT_FAILURE, "Invalid tag name %s\n", tag.c_str());

  it->second.deviceElapsed += device_.timeBetween(it->second.startTag, stopTag);
  it->second.count++;
}

void timer_t::hostTic(const std::string tag, int ifSync)
{
  if (!enabled) {
    return;
  }
  if (ifSync) {
    sync();
  }
  m_[tag].startTime = MPI_Wtime();
}

void timer_t::hostToc(const std::string tag)
{
  if (!enabled) {
    return;
  }
  double stopTime = MPI_Wtime();

  auto it = m_.find(tag);
  nekrsCheck(it == m_.end(), MPI_COMM_SELF, EXIT_FAILURE, "Invalid tag name %s\n", tag.c_str());

  it->second.hostElapsed += (stopTime - it->second.startTime);
  it->second.count++;
}

void timer_t::tic(const std::string tag, int ifSync)
{
  if (!enabled) {
    return;
  }
  if (ifSync) {
    sync();
  }
  m_[tag].startTime = MPI_Wtime();
  m_[tag].startTag = device_.tagStream();
}

void timer_t::toc(const std::string tag)
{
  if (!enabled) {
    return;
  }
  auto stopTime = MPI_Wtime();
  auto stopTag = device_.tagStream();

  auto it = m_.find(tag);
  nekrsCheck(it == m_.end(), MPI_COMM_SELF, EXIT_FAILURE, "Invalid tag name %s\n", tag.c_str());

  it->second.hostElapsed += (stopTime - it->second.startTime);
  it->second.deviceElapsed += device_.timeBetween(it->second.startTag, stopTag);
  it->second.count++;
}

double timer_t::hostElapsed(const std::string tag)
{
  auto it = m_.find(tag);
  if (it == m_.end()) {
    return NEKRS_TIMER_INVALID_KEY;
  }
  return it->second.hostElapsed;
}

double timer_t::deviceElapsed(const std::string tag)
{
  auto it = m_.find(tag);
  if (it == m_.end()) {
    return NEKRS_TIMER_INVALID_KEY;
  }
  return it->second.deviceElapsed;
}

long long int timer_t::count(const std::string tag)
{
  auto it = m_.find(tag);
  if (it == m_.end()) {
    return NEKRS_TIMER_INVALID_KEY;
  }
  return it->second.count;
}

double timer_t::query(const std::string tag, const std::string metric)
{
  int size;
  MPI_Comm_size(comm_, &size);

  auto it = m_.find(tag);
  if (it == m_.end()) {
    return NEKRS_TIMER_INVALID_KEY;
  }
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
      if (child) {
        std::cout << "  " << printPercentage(tTag, tNorm);
      } else {
        std::cout << "      ";
      }
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
      if (child) {
        std::cout << "  " << printPercentage(tTag, tNorm);
      } else {
        std::cout << "      ";
      }
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
      if (child) {
        std::cout << "  " << printPercentage(time, tNorm);
      } else {
        std::cout << "      ";
      }
      std::cout << "\n";
    }
  }
}

void timer_t::print(std::string timerName, long long int DOF)
{
  const auto timerTags = tags();

  // filter out tags that do not start with timerName
  std::vector<std::string> filteredTags;
  std::copy_if(timerTags.begin(),
               timerTags.end(),
               std::back_inserter(filteredTags),
               [&](const std::string &tag) { return tag.find(timerName) == 0; });

  // construct tree where parent entries are the portion of the tag left of the last '::'
  std::map<std::string, std::vector<std::string>> tree;

  for (auto &&tag : filteredTags) {
    auto pos = tag.rfind("::");
    if (pos == std::string::npos) {
      tree[""].push_back(tag);
    } else {
      auto parent = tag.substr(0, pos);
      tree[parent].push_back(tag);
    }
  }

  auto pos = timerName.rfind("::");
  const auto start = timerName.substr(0, pos);

  std::ios oldState(nullptr);
  oldState.copyfmt(std::cout);

  // gather timer information from tree
  std::vector<std::string> operations;
  std::vector<std::string> times;
  std::vector<std::string> calls;
  std::vector<std::string> relPercentage;
  std::vector<std::string> absPercentage;
  std::vector<std::string> throughputs;

  std::cout.setf(std::ios::fixed);
  std::function<void(std::string, std::string, std::string, int)> gatherTreeStats;
  gatherTreeStats = [&](std::string tag, std::string rootTag, std::string parentTag, int level) {
    if (level > 0) {
      if (level == 1) {
        rootTag = tag; // set as root of timer tree
      }
      const auto tTag = platform->timer.query(tag, "DEVICE:MAX");
      const auto nCalls = platform->timer.count(tag);

      if (nCalls == 0) {
        return; // nothing to print
      }

      auto tParent = platform->timer.query(parentTag, "DEVICE:MAX");

      if (tParent < 0.0) {
        tParent = tTag;
      }

      const auto tRoot = platform->timer.query(rootTag, "DEVICE:MAX");

      const auto tCall = tTag / nCalls;

      // trim parentTag from the current tag
      auto pos = tag.rfind(parentTag);
      auto trimmedTag = tag.substr(pos + parentTag.length() + 2);

      if (platform->comm.mpiRank == 0) {
        std::ostringstream ss;
        for (int i = 0; i < level; ++i) {
          ss << "> ";
        }
        ss << trimmedTag;
        operations.push_back(ss.str());

        ss.str("");
        ss.clear();
        ss << std::setprecision(3) << std::scientific << tTag;
        times.push_back(ss.str());

        ss.str("");
        ss.clear();
        ss << std::setw(6) << nCalls;
        calls.push_back(ss.str());

        ss.str("");
        ss.clear();
        ss << std::setprecision(1) << std::fixed << 100.0 * tTag / tParent;
        relPercentage.push_back(level == 1 ? "" : ss.str());

        ss.str("");
        ss.clear();
        ss << std::setprecision(1) << std::fixed << 100.0 * tTag / tRoot;
        absPercentage.push_back(ss.str());

        if (DOF) {
          const double GDOFs = DOF/(1e9 * platform->comm.mpiCommSize);
          ss.str("");
          ss.clear();
          ss << std::setprecision(3) << std::scientific << GDOFs / tCall;
          throughputs.push_back(ss.str());
        }
      }
    }

    std::vector<std::string> children;
    for (auto &&child : tree[tag]) {
      children.push_back(child);
    }

    // sort children by max time, from largest to smallest
    std::sort(children.begin(), children.end(), [&](const std::string &a, const std::string &b) {
      const auto ta = platform->timer.query(a, "DEVICE:MAX");
      const auto tb = platform->timer.query(b, "DEVICE:MAX");
      return ta > tb;
    });

    for (auto &&child : children) {
      gatherTreeStats(child, rootTag, tag, level + 1);
    }
  };

  gatherTreeStats(start, "", "", 0);

  std::map<int, std::vector<std::string>> table;
  table[0] = operations;
  table[1] = times;
  table[2] = calls;
  table[3] = relPercentage;
  table[4] = absPercentage;
  if (DOF) table[5] = throughputs;

  std::vector<std::string> headers = {"name", "time", "calls", "rel %", "abs %"};
  if (DOF) headers.push_back("GDOF/s/rank");

  if (platform->comm.mpiRank == 0) {
    std::cout << "\n";
    std::cout << "timers for " << start << ":\n";
    printTable(table, headers, "    ");
    std::cout << "\n";
  }

  std::cout.copyfmt(oldState);
}

std::vector<std::string> timer_t::tags()
{
  std::vector<std::string> entries;
  auto &keys = m_.keys();
  std::copy(keys.begin(), keys.end(), std::back_inserter(entries));
  return entries;
}

void timer_t::addUserStat(const std::string& tag)
{
  userStat.push_back(tag);
}

void timer_t::printUserStat()
{
  for(const auto& entry : userStat) { 
    platform->timer.print(entry);
  }
}

} // namespace timer

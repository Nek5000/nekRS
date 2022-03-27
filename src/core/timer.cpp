#include <iostream>
#include <string>
#include <map>
#include <algorithm>

#include "timer.hpp"
#include "platform.hpp"
#include "ogs.hpp"


namespace timer
{
namespace
{
struct tagData
{
  long long int count;
  double hostElapsed;
  double deviceElapsed;
  double startTime;
  occa::streamTag startTag;
};
std::map<std::string,tagData> m_;

const int NEKRS_TIMER_INVALID_KEY    = -1;
const int NEKRS_TIMER_INVALID_METRIC = -2;

int ifSync_;
inline int ifSync(){ return ifSync_; }

occa::device device_;
MPI_Comm comm_;
}

timer_t::timer_t(MPI_Comm comm,occa::device device,int ifSync)
{
  init(comm, device, ifSync);
}
void timer_t::init(MPI_Comm comm,occa::device device,int ifSync)
{
  device_ = device;
  ifSync_ = ifSync;
  comm_ = comm;
}

void timer_t::set(const std::string tag, double time, long long int count)
{
  m_[tag].startTime = time;	
  auto it = m_.find(tag);
  if(it == m_.end()) {
    printf("Error in set: Invalid tag name. %s:%u\n",__FILE__,__LINE__);
    MPI_Abort(comm_,1);
  }

  it->second.hostElapsed = time; 
  it->second.deviceElapsed = it->second.hostElapsed;
  it->second.count = count;
}

void timer_t::reset()
{
  m_.clear();
  ogsResetTime();
}

void timer_t::reset(const std::string tag)
{
  std::map<std::string,tagData>::iterator it = m_.find(tag);
  if(it != m_.end()) m_.erase(it);
}

void timer_t::finalize()
{
  reset();
}

void timer_t::deviceTic(const std::string tag,int ifSync)
{
  if(ifSync) MPI_Barrier(comm_);
  m_[tag].startTag = device_.tagStream();
}

void timer_t::deviceTic(const std::string tag)
{
  if(ifSync()) MPI_Barrier(comm_);
  m_[tag].startTag = device_.tagStream();
}

void timer_t::deviceToc(const std::string tag)
{
  occa::streamTag stopTag = device_.tagStream();

  std::map<std::string,tagData>::iterator it = m_.find(tag);
  if(it == m_.end()) {
    printf("Error in deviceToc: Invalid tag name. %s:%u\n",__FILE__,__LINE__);
    MPI_Abort(comm_,1);
  }

  it->second.deviceElapsed += device_.timeBetween(it->second.startTag,stopTag);
  it->second.count++;
}

void timer_t::hostTic(const std::string tag,int ifSync)
{
  if(ifSync) MPI_Barrier(comm_);
  m_[tag].startTime = MPI_Wtime();
}

void timer_t::hostTic(const std::string tag)
{
  if(ifSync()) MPI_Barrier(comm_);
  m_[tag].startTime = MPI_Wtime();
}

void timer_t::hostToc(const std::string tag)
{
  double stopTime = MPI_Wtime();

  auto it = m_.find(tag);
  if(it == m_.end()) {
    printf("Error in deviceToc: Invalid tag name. %s:%u\n",__FILE__,__LINE__);
    MPI_Abort(comm_,1);
  }

  it->second.hostElapsed += (stopTime - it->second.startTime);
  it->second.count++;
}

void timer_t::tic(const std::string tag,int ifSync)
{
  if(ifSync) MPI_Barrier(comm_);
  m_[tag].startTime = MPI_Wtime();
  m_[tag].startTag = device_.tagStream();
}

void timer_t::tic(const std::string tag)
{
  if(ifSync()) MPI_Barrier(comm_);
  m_[tag].startTime = MPI_Wtime();
  m_[tag].startTag = device_.tagStream();
}

void timer_t::toc(const std::string tag)
{
  auto stopTime = MPI_Wtime();
  auto stopTag = device_.tagStream();

  auto it = m_.find(tag);
  if(it == m_.end()) {
    printf("Error in deviceToc: Invalid tag name. %s:%u\n",__FILE__,__LINE__);
    MPI_Abort(comm_,1);
  }

  it->second.hostElapsed  += (stopTime - it->second.startTime);
  it->second.deviceElapsed += device_.timeBetween(it->second.startTag,stopTag);
  it->second.count++;
}

double timer_t::hostElapsed(const std::string tag)
{
  auto it = m_.find(tag);
  if(it == m_.end()) return NEKRS_TIMER_INVALID_KEY;
  return it->second.hostElapsed;
}

double timer_t::deviceElapsed(const std::string tag)
{
  auto it = m_.find(tag);
  if(it == m_.end()) return NEKRS_TIMER_INVALID_KEY;
  return it->second.deviceElapsed;
}

long long int timer_t::count(const std::string tag)
{
  auto it = m_.find(tag);
  if(it == m_.end()) return NEKRS_TIMER_INVALID_KEY;
  return it->second.count;
}

double timer_t::query(const std::string tag,const std::string metric)
{
  int size;
  MPI_Comm_size(comm_,&size);

  auto it = m_.find(tag);
  if(it == m_.end()) return NEKRS_TIMER_INVALID_KEY;
  auto hostElapsed = it->second.hostElapsed;
  auto deviceElapsed = it->second.deviceElapsed;
  auto count = it->second.count;

  double retVal;

  std::string upperMetric = metric;
  std::transform(upperMetric.begin(),upperMetric.end(),upperMetric.begin(),::toupper);

  if(upperMetric.compare("HOST:MIN"  ) == 0) {
    MPI_Allreduce(&hostElapsed,&retVal,1,MPI_DOUBLE,MPI_MIN,comm_);
    return retVal;
  }
  if(upperMetric.compare("HOST:MAX"  ) == 0) {
    MPI_Allreduce(&hostElapsed,&retVal,1,MPI_DOUBLE,MPI_MAX,comm_);
    return retVal;
  }
  if(upperMetric.compare("HOST:SUM"  ) == 0) {
    MPI_Allreduce(&hostElapsed,&retVal,1,MPI_DOUBLE,MPI_SUM,comm_);
    return retVal;
  }
  if(upperMetric.compare("HOST:AVG"  ) == 0) {
    MPI_Allreduce(&hostElapsed,&retVal,1,MPI_DOUBLE,MPI_SUM,comm_);
    return retVal / (size * count);
  }
  if(upperMetric.compare("DEVICE:MIN") == 0) {
    MPI_Allreduce(&deviceElapsed,&retVal,1,MPI_DOUBLE,MPI_MIN,comm_);
    return retVal;
  }
  if(upperMetric.compare("DEVICE:MAX") == 0) {
    MPI_Allreduce(&deviceElapsed,&retVal,1,MPI_DOUBLE,MPI_MAX,comm_);
    return retVal;
  }
  if(upperMetric.compare("DEVICE:SUM") == 0) {
    MPI_Allreduce(&deviceElapsed,&retVal,1,MPI_DOUBLE,MPI_SUM,comm_);
    return retVal;
  }
  if(upperMetric.compare("DEVICE:AVG") == 0) {
    MPI_Allreduce(&deviceElapsed,&retVal,1,MPI_DOUBLE,MPI_SUM,comm_);
    return retVal / (size * count);
  }
  return NEKRS_TIMER_INVALID_METRIC;
}

std::string printPercentage(double num, double dom)
{
  char buf[4096];
  double frac = num/dom;
  snprintf(buf, sizeof(buf), "%7.2f", frac);
  return std::string(buf);
}

void timer_t::printStatEntry(std::string name, std::string tag, std::string type) 
{
  int rank;
  MPI_Comm_rank(comm_, &rank);
  const double tTag = query(tag, type);
  const double tSolve = query("solve", type);
  const long long int nCalls = count(tag);
  if(tTag > 0) {
    if(rank == 0){
      std::cout << name 
                << tTag << "s"  
                << "  " << printPercentage(tTag, tSolve)
                << "  " << nCalls << "\n";
    }
  } 
}

void timer_t::printStatEntry(std::string name, double time) 
{
  int rank;
  MPI_Comm_rank(comm_, &rank);
  const double tSolve = query("solve", "DEVICE:MAX");
  if(time > 0) {
    if(rank == 0){
      std::cout << name 
                << time << "s"  
                << "  " << printPercentage(time, tSolve) << "\n";
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

  if(rank == 0) std::cout << "\n>>> runtime statistics (step = " << step << "):\n";

  std::cout.setf(std::ios::scientific);
  int outPrecisionSave = std::cout.precision();
  std::cout.precision(5);

  if(rank == 0) std::cout <<   "name                    " << "time          " << "time(%)  " << "calls\n";

  printStatEntry("  setup                 ", "setup", "DEVICE:MAX");

  printStatEntry("    loadKernels         ", "loadKernels", "HOST:MAX");

  printStatEntry("  checkpointing         ", "checkpointing", "DEVICE:MAX");

  printStatEntry("  udfExecuteStep        ", "udfExecuteStep", "DEVICE:MAX");

  const double tSolve        = query("solve", "DEVICE:MAX");
  const double tMinSolveStep = query("minSolveStep", "HOST:MAX");
  const double tMaxSolveStep = query("maxSolveStep", "HOST:MAX");
  const double flops = platform->flopCounter->get(platform->comm.mpiComm);
  bool printFlops = !platform->options.compareArgs("PRESSURE PRECONDITIONER", "SEMFEM");

  if(tSolve > 0 && rank == 0) {
    std::cout << "  total solve           " << tSolve << "s\n";
    std::cout << "    step min            " << tMinSolveStep << "s\n";
    std::cout << "    step max            " << tMaxSolveStep << "s\n";
    if (flops > 0 && printFlops)
    std::cout << "    FLOPS/s             " << flops/tSolve << "\n";
  }

  printStatEntry("    meshUpdate          ", "meshUpdate", "DEVICE:MAX");
  printStatEntry("    makef               ", "makef", "DEVICE:MAX");
  printStatEntry("      udfUEqnSource     ", "udfUEqnSource", "DEVICE:MAX");

  printStatEntry("    makeq               ", "makeq", "DEVICE:MAX");
  printStatEntry("      udfSEqnSource     ", "udfSEqnSource", "DEVICE:MAX");

  printStatEntry("    udfProperties       ", "udfProperties", "DEVICE:MAX");
 

  printStatEntry("    velocitySolve       ", "velocitySolve", "DEVICE:MAX");
  printStatEntry("      projection        ", "velocity proj", "DEVICE:MAX");

  printStatEntry("    pressureSolve       ", "pressureSolve", "DEVICE:MAX");
  printStatEntry("      preconditioner    ", "pressure preconditioner", "DEVICE:MAX");
  printStatEntry("        pMG smoother    ", "pressure preconditioner smoother", "DEVICE:MAX");
  printStatEntry("        coarse grid     ", "coarseSolve", "DEVICE:MAX");
  printStatEntry("      projection        ", "pressure proj", "DEVICE:MAX");

  printStatEntry("    scalarSolve         ", "scalarSolve", "DEVICE:MAX");
  printStatEntry("      projection        ", "scalar proj", "DEVICE:MAX");

  printStatEntry("    meshSolve           ", "meshSolve", "DEVICE:MAX");

  printStatEntry("    gsMPI               ", gsTime);

  printStatEntry("    dotp                ", "dotp", "DEVICE:MAX");

  if(rank == 0) std::cout << std::endl;

  std::cout.unsetf(std::ios::scientific);
  std::cout.precision(outPrecisionSave);
}

} // namespace



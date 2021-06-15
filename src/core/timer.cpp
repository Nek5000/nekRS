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
  int count;
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

void timer_t::set(const std::string tag, double time)
{
  m_[tag].startTime = time;	
  auto it = m_.find(tag);
  if(it == m_.end()) {
    printf("Error in set: Invalid tag name. %s:%u\n",__FILE__,__LINE__);
    MPI_Abort(comm_,1);
  }

  it->second.hostElapsed = time; 
  it->second.deviceElapsed = it->second.hostElapsed;
  it->second.count++;
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

int timer_t::count(const std::string tag)
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

void timer_t::printRunStat()
{
  int rank;
  MPI_Comm_rank(comm_, &rank);

  double dEtime[20];
  dEtime[0] = query("makef", "DEVICE:MAX");
  dEtime[1] = query("velocitySolve", "DEVICE:MAX");
  dEtime[17] = query("velocity proj pre", "DEVICE:MAX");
  dEtime[17]+= query("velocity proj post", "DEVICE:MAX");
  dEtime[2] = query("pressureSolve", "DEVICE:MAX");
  dEtime[3] = query("makeq", "DEVICE:MAX");
  dEtime[4] = query("scalarSolve", "DEVICE:MAX");
  dEtime[5] = query("pressure preconditioner", "DEVICE:MAX");
  dEtime[16] = query("pressure preconditioner smoother", "DEVICE:MAX");
  dEtime[6] = query("pressure proj pre", "DEVICE:MAX");
  dEtime[6]+= query("pressure proj post", "DEVICE:MAX");

  dEtime[8] = query("dotp", "DEVICE:MAX");

  dEtime[9] = query("solve", "DEVICE:MAX");
  dEtime[10] = query("setup", "DEVICE:MAX");
  dEtime[11] = query("checkpointing", "DEVICE:MAX");

  dEtime[12] = query("udfExecuteStep", "DEVICE:MAX");
  dEtime[13] = query("udfUEqnSource", "DEVICE:MAX");
  dEtime[14] = query("udfSEqnSource", "DEVICE:MAX");
  dEtime[15] = query("udfProperties", "DEVICE:MAX");

  double hEtime[10];
  hEtime[0] = query("BoomerAMGSolve", "HOST:MAX");
  const double amgxTime = query("AmgXSolve", "DEVICE:MAX");
  hEtime[0] = hEtime[0] > amgxTime ? hEtime[0] : amgxTime;
  hEtime[1] = ogsTime(/* reportHostTime */ true);
  MPI_Allreduce(MPI_IN_PLACE, &hEtime[1], 1, MPI_DOUBLE, MPI_MAX, comm_);

  if (rank == 0) {
    std::cout.setf ( std::ios::scientific );

    std::cout << "runtime statistics\n\n"
  	      << "  setup                 " << dEtime[10]<< " s\n";

    if(dEtime[11] > 0)
    std::cout << "  checkpointing         " << dEtime[11]<< " s\n";

    if(dEtime[12] > 0)
    std::cout << "  udfExecuteStep        " << dEtime[12] << " s\n";

    std::cout << "  total solve           " << dEtime[9] << " s\n"
  	      << "    makef               " << dEtime[0] << " s\n";
    if(dEtime[13] > 0)
    std::cout << "      udfUEqnSource     " << dEtime[13] << " s\n";
    std::cout << "    velocitySolve       " << dEtime[1] << " s\n";
    if(dEtime[17] > 0)
    std::cout << "      projection        " << dEtime[17] << " s\n";

    std::cout << "    pressureSolve       " << dEtime[2] << " s\n"
              << "      preconditioner    " << dEtime[5] << " s\n";
    if(dEtime[16] > 0)
    std::cout << "        pMG smoother    " << dEtime[16] << " s\n";
    if(hEtime[0] > 0)
    std::cout << "        coarse grid     " << hEtime[0] << " s\n";
    if(dEtime[6] > 0)
    std::cout << "      projection        " << dEtime[6] << " s\n";

    if(dEtime[14] > 0)
    std::cout << "    scalarSolve         " << dEtime[4] << " s\n"
              << std::endl;
    if(dEtime[4] > 0) {
    std::cout << "    makeq               " << dEtime[3] << " s\n";
     if(dEtime[14] > 0)
    std::cout << "      udfSEqnSource     " << dEtime[14] << " s\n";
    }

    if(dEtime[15] > 0)
    std::cout << "    udfProperties       " << dEtime[15] << " s\n"
              << std::endl;

    if(hEtime[1] > 0)
    std::cout << "    gsMPI               " << hEtime[1] << " s\n";
    if(dEtime[8] > 0)
    std::cout << "    dotp                " << dEtime[8] << " s\n";

    std::cout << std::endl;

    std::cout.unsetf ( std::ios::scientific );
  }
}
} // namespace

#include <iostream>
#include <string>
#include <map>
#include <algorithm>

#include "occa.hpp"
#include "mpi.h"
#include "timer.hpp"

namespace timer {

namespace {

typedef struct tagData_{
  int count;
  double hostElapsed;
  double deviceElapsed;
  double startTime;
  occa::streamTag startTag;
} tagData;
std::map<std::string,tagData> m_;

const int NEKRS_TIMER_INVALID_KEY    = -1;
const int NEKRS_TIMER_INVALID_METRIC = -2;

int ifSync_;
inline int ifSync(){ return ifSync_; }

occa::device device_;
MPI_Comm comm_;
}

void init(MPI_Comm comm,occa::device device,int ifSync){
  device_=device;
  ifSync_=ifSync;
  comm_ = comm;
}

void reset(){
  m_.clear();
}

void reset(const std::string tag){
  std::map<std::string,tagData>::iterator it=m_.find(tag);
  if(it!=m_.end()) m_.erase(it);
}

void finalize(){
  reset();
}

void deviceTic(const std::string tag,int ifSync){
  if(ifSync) MPI_Barrier(comm_);
  m_[tag].startTag =device_.tagStream();
}

void deviceTic(const std::string tag){
  if(ifSync()) MPI_Barrier(comm_);
  m_[tag].startTag =device_.tagStream();
}

void deviceToc(const std::string tag){
  occa::streamTag stopTag=device_.tagStream();

  std::map<std::string,tagData>::iterator it=m_.find(tag);
  if(it==m_.end()){
    printf("Error in deviceToc: Invalid tag name. %s:%u\n",__FILE__,__LINE__);
    MPI_Abort(comm_,1);
  }

  it->second.deviceElapsed+=device_.timeBetween(it->second.startTag,stopTag);
  it->second.count++;
}

void hostTic(const std::string tag,int ifSync){
  if(ifSync) MPI_Barrier(comm_);
  m_[tag].startTime=MPI_Wtime();
}

void hostTic(const std::string tag){
  if(ifSync()) MPI_Barrier(comm_);
  m_[tag].startTime=MPI_Wtime();
}

void hostToc(const std::string tag){
  double stopTime=MPI_Wtime();

  auto it=m_.find(tag);
  if(it==m_.end()){
    printf("Error in deviceToc: Invalid tag name. %s:%u\n",__FILE__,__LINE__);
    MPI_Abort(comm_,1);
  }

  it->second.hostElapsed+=(stopTime - it->second.startTime);
  it->second.count++;
}

void tic(const std::string tag,int ifSync){
  if(ifSync) MPI_Barrier(comm_);
  m_[tag].startTime=MPI_Wtime();
  m_[tag].startTag =device_.tagStream();
}

void tic(const std::string tag){
  if(ifSync()) MPI_Barrier(comm_);
  m_[tag].startTime=MPI_Wtime();
  m_[tag].startTag =device_.tagStream();
}

void toc(const std::string tag){
  auto stopTime=MPI_Wtime();
  auto stopTag =device_.tagStream();

  auto it=m_.find(tag);
  if(it==m_.end()){
    printf("Error in deviceToc: Invalid tag name. %s:%u\n",__FILE__,__LINE__);
    MPI_Abort(comm_,1);
  }

  it->second.hostElapsed  +=(stopTime - it->second.startTime);
  it->second.deviceElapsed+=device_.timeBetween(it->second.startTag,stopTag);
  it->second.count++;
}

double hostElapsed(const std::string tag){
  auto it=m_.find(tag);
  if(it==m_.end()){ return NEKRS_TIMER_INVALID_KEY; }
  return it->second.hostElapsed;
}

double deviceElapsed(const std::string tag){
  auto it=m_.find(tag);
  if(it==m_.end()){ return NEKRS_TIMER_INVALID_KEY; }
  return it->second.deviceElapsed;
}

int count(const std::string tag){
  auto it=m_.find(tag);
  if(it==m_.end()){ return NEKRS_TIMER_INVALID_KEY; }
  return it->second.count;
}

double query(const std::string tag,const std::string metric){
  int size;
  MPI_Comm_size(comm_,&size);

  auto it=m_.find(tag);
  if(it==m_.end()){ return NEKRS_TIMER_INVALID_KEY; }
  auto hostElapsed=it->second.hostElapsed;
  auto deviceElapsed=it->second.deviceElapsed;
  auto count=it->second.count;

  double retVal;

  std::string upperMetric=metric;
  std::transform(upperMetric.begin(),upperMetric.end(),upperMetric.begin(),::toupper);

  if(upperMetric.compare("HOST:MIN"  )==0) {
    MPI_Allreduce(&hostElapsed,&retVal,1,MPI_DOUBLE,MPI_MIN,comm_);
    return retVal;
  }
  if(upperMetric.compare("HOST:MAX"  )==0) {
    MPI_Allreduce(&hostElapsed,&retVal,1,MPI_DOUBLE,MPI_MAX,comm_);
    return retVal;
  }
  if(upperMetric.compare("HOST:SUM"  )==0) {
    MPI_Allreduce(&hostElapsed,&retVal,1,MPI_DOUBLE,MPI_SUM,comm_);
    return retVal;
  }
  if(upperMetric.compare("HOST:AVG"  )==0) {
    MPI_Allreduce(&hostElapsed,&retVal,1,MPI_DOUBLE,MPI_SUM,comm_);
    return retVal/(size*count);
  }
  if(upperMetric.compare("DEVICE:MIN")==0) {
    MPI_Allreduce(&deviceElapsed,&retVal,1,MPI_DOUBLE,MPI_MIN,comm_);
    return retVal;
  }
  if(upperMetric.compare("DEVICE:MAX")==0) {
    MPI_Allreduce(&deviceElapsed,&retVal,1,MPI_DOUBLE,MPI_MAX,comm_);
    return retVal;
  }
  if(upperMetric.compare("DEVICE:SUM")==0) {
    MPI_Allreduce(&deviceElapsed,&retVal,1,MPI_DOUBLE,MPI_SUM,comm_);
    return retVal;
  }
  if(upperMetric.compare("DEVICE:AVG")==0) {
    MPI_Allreduce(&deviceElapsed,&retVal,1,MPI_DOUBLE,MPI_SUM,comm_);
    return retVal/(size*count);
  }
  return NEKRS_TIMER_INVALID_METRIC;
}

void printStat()
{
  int rank;
  MPI_Comm_rank(comm_, &rank);

  double dEtime[6];
  dEtime[0] = timer::query("makef", "DEVICE:MAX");
  dEtime[1] = timer::query("velocitySolve", "DEVICE:MAX");
  dEtime[2] = timer::query("pressureSolve", "DEVICE:MAX");
  dEtime[3] = timer::query("makeq", "DEVICE:MAX"); 
  dEtime[4] = timer::query("scalarSolve", "DEVICE:MAX");

  double hEtime[6];
  hEtime[0] = timer::query("makef", "HOST:MAX");
  hEtime[1] = timer::query("velocitySolve", "HOST:MAX");
  hEtime[2] = timer::query("pressureSolve", "HOST:MAX");
  hEtime[3] = timer::query("makeq", "HOST:MAX"); 
  hEtime[4] = timer::query("scalarSolve", "HOST:MAX");

  if (rank == 0) {
    std::cout.setf ( std::ios::scientific );

    std::cout << "runtime statistics\n\n";
    std::cout << "                      device(s)      host(s)\n"
              << "  makef               " << dEtime[0] << "   " << hEtime[0] << "\n"  
              << "  velocitySolve       " << dEtime[1] << "   " << hEtime[1] << "\n"  
              << "  pressureSolve       " << dEtime[2] << "   " << hEtime[2] << "\n" 
              << "  makeq               " << dEtime[3] << "   " << hEtime[3] << "\n"
              << "  scalarSolve         " << dEtime[4] << "   " << hEtime[4] << "\n" 
              << std::endl; 

    std::cout.unsetf ( std::ios::scientific ); 
  }
}

} // namespace

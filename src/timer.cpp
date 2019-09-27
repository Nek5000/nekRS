#include "timer.hpp"

#include <algorithm>
#include <string>

int timer::init(MPI_Comm comm,occa::device &device,int ifSync){
  MPI_Comm_dup(comm,&comm_);
  device_=device;
  ifSync_=ifSync;
}

int timer::reset(){
  ifSync_=0;
  m_.clear();
}

int timer::reset(const std::string tag){
  auto it=m_.find(tag);
  if(it!=m_.end()) m_.erase(it);
  return 0;
}

int timer::finalize(){
  reset();
}

int timer::deviceTic(const std::string tag){
  tagData data={.count=0, .hostElapsed=0, .deviceElapsed=0};

  auto it=m_.find(tag);
  if(it==m_.end()){ m_.insert(std::pair<std::string,tagData>(tag,data)); it=m_.find(tag); }

  if(ifSync()) MPI_Barrier(comm_);
  it->second.startTag =device_.tagStream();
}

int timer::deviceToc(const std::string tag){
  auto stopTag=device_.tagStream();

  auto it=m_.find(tag);
  if(it==m_.end()){ return NEKRS_TIMER_INVALID_KEY; }

  it->second.deviceElapsed=device_.timeBetween(it->second.startTag,stopTag);
  it->second.count++;
  return 0;
}

int timer::hostTic(const std::string tag){
  tagData data={.count=0, .hostElapsed=0, .deviceElapsed=0};

  auto it=m_.find(tag);
  if(it==m_.end()){ m_.insert(std::pair<std::string,tagData>(tag,data)); it=m_.find(tag); }

  if(ifSync()) MPI_Barrier(comm_);
  it->second.startTime=MPI_Wtime();
}

int timer::hostToc(const std::string tag){
  auto stopTime=MPI_Wtime();

  auto it=m_.find(tag);
  if(it==m_.end()){ return NEKRS_TIMER_INVALID_KEY; }

  it->second.hostElapsed+=(stopTime-it->second.startTime);
  it->second.count++;
  return 0;
}

int timer::tic(const std::string tag){
  tagData data={.count=0, .hostElapsed=0, .deviceElapsed=0};

  auto it=m_.find(tag);
  if(it==m_.end()){ m_.insert(std::pair<std::string,tagData>(tag,data)); it=m_.find(tag); }

  if(ifSync()) MPI_Barrier(comm_);
  it->second.startTime=MPI_Wtime();
  it->second.startTag =device_.tagStream();
}

int timer::toc(const std::string tag){
  auto stopTime=MPI_Wtime();
  auto stopTag =device_.tagStream();

  auto it=m_.find(tag);
  if(it==m_.end()){ return NEKRS_TIMER_INVALID_KEY; }

  it->second.hostElapsed +=(stopTime-it->second.startTime);
  it->second.deviceElapsed=device_.timeBetween(it->second.startTag,stopTag);
  it->second.count++;
  return 0;
}

double timer::hostElapsed(const std::string tag){
  auto it=m_.find(tag);
  if(it==m_.end()){ return NEKRS_TIMER_INVALID_KEY; }
  return it->second.hostElapsed;
}

double timer::deviceElapsed(const std::string tag){
  auto it=m_.find(tag);
  if(it==m_.end()){ return NEKRS_TIMER_INVALID_KEY; }
  return it->second.deviceElapsed;
}

int timer::count(const std::string tag){
  auto it=m_.find(tag);
  if(it==m_.end()){ return NEKRS_TIMER_INVALID_KEY; }
  return it->second.count;
}

timer::tagStats timer::query(const std::string tag){
  auto it=m_.find(tag);
  if(it==m_.end()){ tagStats t={0,0,0,0,0,0,0}; return t; }

  auto hostElapsed=it->second.hostElapsed;
  auto deviceElapsed=it->second.deviceElapsed;
  tagStats t;

  MPI_Allreduce(&hostElapsed  ,&t.hostMin  ,1,MPI_DOUBLE,MPI_MIN,comm_);
  MPI_Allreduce(&hostElapsed  ,&t.hostMax  ,1,MPI_DOUBLE,MPI_MAX,comm_);
  MPI_Allreduce(&hostElapsed  ,&t.hostSum  ,1,MPI_DOUBLE,MPI_SUM,comm_);
  MPI_Allreduce(&deviceElapsed,&t.deviceMin,1,MPI_DOUBLE,MPI_MIN,comm_);
  MPI_Allreduce(&deviceElapsed,&t.deviceMax,1,MPI_DOUBLE,MPI_MAX,comm_);
  MPI_Allreduce(&deviceElapsed,&t.deviceSum,1,MPI_DOUBLE,MPI_SUM,comm_);
  t.count=it->second.count;

  return t;
}

double timer::query(const std::string tag,const std::string metric){
  std::string upperMetric=metric;
  std::transform(upperMetric.begin(),upperMetric.end(),upperMetric.begin(),::toupper);

  tagStats t=query(tag);
  if(upperMetric.compare("HOST:MIN")==0)   return t.hostMin;
  if(upperMetric.compare("HOST:MAX")==0)   return t.hostMax;
  if(upperMetric.compare("HOST:SUM")==0)   return t.hostSum;
  if(upperMetric.compare("HOST:AVG")==0)   return t.hostSum/t.count;
  if(upperMetric.compare("DEVICE:MIN")==0) return t.deviceMin;
  if(upperMetric.compare("DEVICE:MAX")==0) return t.deviceMax;
  if(upperMetric.compare("DEVICE:SUM")==0) return t.deviceSum;
  if(upperMetric.compare("DEVICE:AVG")==0) return t.deviceSum/t.count;
  return NEKRS_TIMER_INVALID_METRIC;
}

#include "timer.hpp"

#include <algorithm>
#include <string>

void timer::init(MPI_Comm comm,occa::device &device,int ifSync){
  device_=device;
  ifSync_=ifSync;
  if(MPI_Comm_dup(comm,&comm_)!=MPI_SUCCESS){
    printf("Error in MPI_Comm_dup: %s:%u\n",__FILE__,__LINE__);
    MPI_Abort(comm,1);
  }
}

void timer::reset(){
  ifSync_=0;
  m_.clear();
}

void timer::reset(const std::string tag){
  std::map<std::string,tagData>::iterator it=m_.find(tag);
  if(it!=m_.end()) m_.erase(it);
}

void timer::finalize(){
  reset();
}

void timer::deviceTic(const std::string tag){
  tagData data={.count=0, .hostElapsed=0, .deviceElapsed=0};
  m_[tag]=data;

  if(ifSync()) MPI_Barrier(comm_);
  m_[tag].startTag =device_.tagStream();
}

void timer::deviceToc(const std::string tag){
  occa::streamTag stopTag=device_.tagStream();

  std::map<std::string,tagData>::iterator it=m_.find(tag);
  if(it==m_.end()){
    printf("Error in deviceToc: Invalid tag name. %s:%u\n",__FILE__,__LINE__);
    MPI_Abort(comm_,1);
  }

  it->second.deviceElapsed+=device_.timeBetween(it->second.startTag,stopTag);
  it->second.count++;
}

void timer::hostTic(const std::string tag){
  tagData data={.count=0, .hostElapsed=0, .deviceElapsed=0};
  m_[tag]=data;

  if(ifSync()) MPI_Barrier(comm_);
  m_[tag].startTime=MPI_Wtime();
}

void timer::hostToc(const std::string tag){
  double stopTime=MPI_Wtime();

  auto it=m_.find(tag);
  if(it==m_.end()){
    printf("Error in deviceToc: Invalid tag name. %s:%u\n",__FILE__,__LINE__);
    MPI_Abort(comm_,1);
  }

  it->second.hostElapsed+=(stopTime-it->second.startTime);
  it->second.count++;
}

void timer::tic(const std::string tag){
  tagData data={.count=0, .hostElapsed=0, .deviceElapsed=0};
  m_[tag]=data;

  if(ifSync()) MPI_Barrier(comm_);
  m_[tag].startTime=MPI_Wtime();
  m_[tag].startTag =device_.tagStream();
}

void timer::toc(const std::string tag){
  auto stopTime=MPI_Wtime();
  auto stopTag =device_.tagStream();

  auto it=m_.find(tag);
  if(it==m_.end()){
    printf("Error in deviceToc: Invalid tag name. %s:%u\n",__FILE__,__LINE__);
    MPI_Abort(comm_,1);
  }

  it->second.hostElapsed  +=(stopTime-it->second.startTime);
  it->second.deviceElapsed+=device_.timeBetween(it->second.startTag,stopTag);
  it->second.count++;
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

  int size;
  MPI_Comm_size(comm_,&size); t.deviceSum/=size; t.hostSum/=size;

  return t;
}

double timer::query(const std::string tag,const std::string metric){
  std::string upperMetric=metric;
  std::transform(upperMetric.begin(),upperMetric.end(),upperMetric.begin(),::toupper);

  tagStats t=query(tag);
  if(upperMetric.compare("HOST:MIN"  )==0) return t.hostMin;
  if(upperMetric.compare("HOST:MAX"  )==0) return t.hostMax;
  if(upperMetric.compare("HOST:SUM"  )==0) return t.hostSum;
  if(upperMetric.compare("HOST:AVG"  )==0) return t.hostSum/t.count;
  if(upperMetric.compare("DEVICE:MIN")==0) return t.deviceMin;
  if(upperMetric.compare("DEVICE:MAX")==0) return t.deviceMax;
  if(upperMetric.compare("DEVICE:SUM")==0) return t.deviceSum;
  if(upperMetric.compare("DEVICE:AVG")==0) return t.deviceSum/t.count;
  return NEKRS_TIMER_INVALID_METRIC;
}

std::map<std::string,timer::tagData>::const_iterator timer::getTags(){
  return m_.begin();
}

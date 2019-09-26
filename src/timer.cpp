#include "timer.hpp"

int nekrs::timer::timerInit(MPI_Comm comm,occa::device &device,int ifSync){
  MPI_Comm_dup(comm,&comm_);
  device_=device;
  ifSync_=ifSync;
}

int nekrs::timer::timerReset(){
  ifSync_=0;
  m_.clear();
}

int nekrs::timer::timerFinalize(){
  timerReset();
}

int nekrs::timer::deviceTic(std::string tag){
  startTag = device_.tagStream();

  struct tagData data;
  auto it=m_.find(tag);
  if(it==m_.end()){ // new key
    data.count=0; data.elapsed=0; data.host=0;
  } else {
    data=it->second;
    if(data.host==1) { printf("nekrs::timer: deviceTic used with a host tag.\n"); return 1; }
    m_.erase(it);
  }

  data.startTag=startTag;
  m_.insert(std::pair<std::string,struct tagData>(tag,data));
}

int nekrs::timer::hostTic(std::string tag){
  if(ifSync()) MPI_Barrier(comm_);
  startTime=MPI_Wtime();

  struct tagData data;
  auto it=m_.find(tag);
  if(it==m_.end()){ data.count=0; data.elapsed=0; data.host=1; }
  else {
    data=it->second;
    if(data.host==0) { printf("nekrs::timer: hostTic used with a device tag.\n"); return 1; }
    m_.erase(it);
  }

  data.startTime=startTime;
  m_.insert(std::pair<std::string,struct tagData>(tag,data));
}

double nekrs::timer::hostToc(std::string tag){
  if(ifSync()) MPI_Barrier(comm_);
  stopTime=MPI_Wtime();

  struct tagData data;
  auto it=m_.find(tag);
  if(it==m_.end()){ printf("nekrs::timer: invalid key.\n"); return -1; }
  else {
    data=it->second;
    if(data.host==0) { printf("nekrs::timer: hostToc used with a device tag.\n"); return -1; }
    m_.erase(it);
  }

  data.elapsed+= (stopTime - data.startTime);
  data.count++;
  m_.insert(std::pair<std::string,struct tagData>(tag,data));
  return data.elapsed;
}

double nekrs::timer::deviceToc(std::string tag){
  stopTag = device_.tagStream();

  struct tagData data;
  auto it=m_.find(tag);
  if(it==m_.end()){ printf("nekrs::timer: invalid key.\n"); return 1; }
  else {
    data=it->second;
    if(data.host==1) { printf("nekrs::timer: deviceToc used with a host tag.\n"); return 1; }
    m_.erase(it);
  }

  data.elapsed=device_.timeBetween(startTag,stopTag);
  data.count++;
  m_.insert(std::pair<std::string,struct tagData>(tag,data));
  return data.elapsed;
}

double nekrs::timer::elapsed(std::string tag){
  auto it=m_.find(tag);
  if(it==m_.end()){ printf("nekrs::timer invalid key.\n"); return 1;}
  return it->second.elapsed;
}

int nekrs::timer::count(std::string tag){
  auto it=m_.find(tag);
  if(it==m_.end()){ printf("nekrs::timer invalid key.\n"); return 1;}
  return it->second.count;
}

int nekrs::timer::reset(std::string tag){
  auto it=m_.find(tag);
  if(it!=m_.end()) m_.erase(it);
  return 0;
}

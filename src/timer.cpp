#include "timer.hpp"

int timer::timerInit(int ifsync){
  ifsync_=ifsync;
}

int timer::timerInit(occa::device &device,int ifsync){
  device_=device;
  ifdevice_=1;
  timerInit(ifsync);
}

int timer::timerReset(){
  ifsync_=0;
  m_.clear();
}

int timer::timerFinalize(){
  timerReset();
}

static int timer::updateTic(std::string &tag,double time){
  auto it=m_.find(tag);
  if(it==m_.end()){
    // new key
    struct tagData value;
    value.count=0;
    value.ticTime=time;
    value.elapsed=0.0;
    m_.insert(std::pair<std::string,struct tagData>(tag,value));
  } else {
    if(it->second.ticTime>=0) printf("nekrs::timer is not in a valid statement.\n");
    it->second.ticTime=time;
  }
}

static int timer::updateToc(std::string &tag,double time){
  auto it=m_.find(tag);
  if(it==m_.end()){
    printf("nekrs::timer::htoc called before nekrs::timer::htic or "
      "nekrs::timer is not in a valid state.\n");
    return 1;
  } else {
    if(it->second.ticTime<0) printf("nekrs::timer is not in a valid state.\n");
    it->second.elapsed+=(time-it->second.ticTime);
    it->second.count++;
    it->second.ticTime=-1;
  }
}

int timer::tic(std::string &tag){
  if(isSync()) mpi::Barrier();
  if(isDevice()) device_.finish(); 
  return updateTic(tag,mpi::Wtime());
}

int timer::toc(std::string &tag){
  if(isSync()) mpi::Barrier();
  if(isDevice()) device_.finish(); 
  return updateToc(tag,mpi::Wtime());
}

/*
int timer::hostTic(std::string &tag){
  if(isSync()) mpi::Barrier();
  double time=mpi::Wtime();
  return updateTic(tag,time);
}

int timer::hostToc(std::string &tag){
  if(isSync()) mpi::Barrier();
  double time=mpi::Wtime();
  return updateToc(tag,time);
}

int timer::deviceTic(std::string &tag){
  if(isSync()) mpi::Barrier();
  if(isDevice()) device.finish(); 
  double time=mpi::Wtime();
  updateTic(tag,time);
}
int timer::deviceToc(std::string &tag){
  if(isSync()) mpi::Barrier();
  if(isDevice()) device_.finish(); 
  double time=mpi::Wtime();
  updateToc(tag,time);
}
*/

double timer::elapsed(std::string &tag){
  auto it=m_.find(tag);
  if(it==m_.end()){ printf("nekrs::timer invalid key.\n"); return 1;}
  return it->second.elapsed;
}

int timer::count(std::string &tag){
  auto it=m_.find(tag);
  if(it==m_.end()){ printf("nekrs::timer invalid key.\n"); return 1;}
  return it->second.count;
}

int timer::reset(std::string &tag){
  auto it=m_.find(tag);
  if(it!=m_.end()) m_.erase(it);
  return 0;
}

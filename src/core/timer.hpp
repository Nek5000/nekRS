#if !defined(nekrs_timer_hpp_)
#define nekrs_timer_hpp_

#include <string>

#include "nrssys.hpp"

namespace timer{
struct timer_t
{
private:
  typedef struct tagData_
  {
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
public:
  timer_t(){}
  void init(MPI_Comm comm,occa::device device,int ifsync = 0);
  void reset();
  void reset(const std::string tag);
  void finalize();
  
  void tic(const std::string tag);
  void tic(const std::string tag,int ifSync);
  void toc(const std::string tag);
  void hostTic(const std::string tag);
  void hostTic(const std::string tag,int ifSync);
  void hostToc(const std::string tag);
  void deviceTic(const std::string tag);
  void deviceTic(const std::string tag,int ifSync);
  void deviceToc(const std::string tag);
  
  void set(const std::string tag, double time);
  
  double hostElapsed(const std::string tag);
  double deviceElapsed(const std::string tag);
  int count(const std::string tag);
  double query(const std::string tag,std::string metric);
  void printRunStat();
};
}

#endif

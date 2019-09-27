#if !defined(nekrs_timer_hpp_)
#define nekrs_timer_hpp_

#include <string>
#include <map>

#include "occa.hpp"
#include "mpi.h"

namespace timer{
  const int NEKRS_TIMER_INVALID_KEY   =-1;
  const int NEKRS_TIMER_INVALID_METRIC=-2;

  typedef struct tagData_{
    int count;
    double hostElapsed;
    double deviceElapsed;
    double startTime;
    occa::streamTag startTag;
  } tagData;

  typedef struct tagStats_{
    double hostMin;
    double hostMax;
    double hostSum;
    double deviceMin;
    double deviceMax;
    double deviceSum;
    int count;
  } tagStats;

  /* private variables */
  static int ifSync_;
  static occa::device device_;
  static MPI_Comm comm_;
  static std::map<std::string,tagData> m_;

  inline int ifSync(){ return ifSync_; }

  int init(MPI_Comm comm,occa::device &device,int ifsync=1);
  int reset();
  int reset(const std::string tag);
  int finalize();

  int tic(const std::string tag);
  int toc(const std::string tag);
  int hostTic(const std::string tag);
  int hostToc(const std::string tag);
  int deviceTic(const std::string tag);
  int deviceToc(const std::string tag);

  double hostElapsed(const std::string tag);
  double deviceElapsed(const std::string tag);
  int count(const std::string tag);
  tagStats query(const std::string tag);
  double query(const std::string tag,std::string metric);
}

#endif

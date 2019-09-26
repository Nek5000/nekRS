#if !defined(nekrs_timer_hpp_)
#define nekrs_timer_hpp_

#include <string>
#include <map>
#include <stack>

#include "occa.hpp"
#include "mpi.h"

namespace nekrs::timer{
  struct tagData{
    int count;
    double elapsed;
    double startTime;
    occa::streamTag startTag;
    int host;
  };

  /* private variables */
  static int ifSync_;
  static occa::device device_;
  static MPI_Comm comm_;
  static std::map<std::string,struct tagData> m_;
  static occa::streamTag startTag ,stopTag;
  static double          startTime,stopTime;

  inline int ifSync(){ return ifSync_;}

  int timerInit(MPI_Comm comm,int ifsync=1);
  int timerInit(occa::device &device);
  int timerReset();
  int timerFinalize();

  int hostTic(std::string tag);
  int hostToc(std::string tag);
  int deviceTic(std::string tag);
  int deviceToc(std::string tag);

  // Replaced query by elapsed and count;
  double elapsed(std::string tag);
  int count(std::string tag);
  int reset(std::string tag);
}

#endif

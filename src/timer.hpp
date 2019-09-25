#if !defined(nekrs_timer_hpp_)
#define nekrs_timer_hpp_

#include <string>
#include <map>

#include "occa.hpp"
#include "mpi_wrapper.hpp"

namespace timer{
  struct tagData{
    int count;
    double elapsed;
    double ticTime;
  };

  /* private variables */
  static int ifsync_,ifdevice_;
  static occa::device device_;
  static std::map<std::string,struct tagData> m_;
  /* private functions */
  static int updateTic(std::string &tag,double time);
  static int updateToc(std::string &tag,double time);

  inline int isSync(){ return ifsync_;}
  inline int isDevice(){ return ifdevice_;}

  int timerInit(int ifsync=1);
  int timerInit(occa::device &device,int ifsync=1);
  int timerReset();
  int timerFinalize();

  int tic(std::string &tag);
  int toc(std::string &tag);

  /*
  int hostTic(std::string &tag);
  int hostToc(std::string &tag);

  int deviceTic(std::string &tag);
  int deviceToc(std::string &tag);
  */

  // Replaced query by elapsed and count;
  double elapsed(std::string &tag);
  int count(std::string &tag);
  int reset(std::string &tag);
}

#endif

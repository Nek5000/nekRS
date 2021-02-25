#if !defined(nekrs_timer_hpp_)
#define nekrs_timer_hpp_

#include <string>

#include "nrssys.hpp"

namespace timer{
void tic(const std::string tag);
void tic(const std::string tag,int ifSync);
void toc(const std::string tag);
void hostTic(const std::string tag);
void hostTic(const std::string tag,int ifSync);
void hostToc(const std::string tag);

struct timer_t
{
timer_t(MPI_Comm comm,occa::device device,int ifsync = 0);
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

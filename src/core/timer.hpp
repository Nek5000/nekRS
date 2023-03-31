#if !defined(nekrs_timer_hpp_)
#define nekrs_timer_hpp_

#include <string>
#include "nrssys.hpp"


namespace timer{
struct timer_t
{
timer_t(MPI_Comm comm,occa::device device,int ifsync, int enable_sync);
void init(MPI_Comm comm,occa::device device,int ifsync, int enable_sync);
void reset();
void reset(const std::string tag);
void finalize();
void enableSync();
void disableSync();
void enable();
void disable();

void tic(const std::string tag);
void tic(const std::string tag,int ifSync);
void toc(const std::string tag);
void hostTic(const std::string tag);
void hostTic(const std::string tag,int ifSync);
void hostToc(const std::string tag);
void deviceTic(const std::string tag);
void deviceTic(const std::string tag,int ifSync);
void deviceToc(const std::string tag);

void set(const std::string tag, double time, long long int count = 1);

double hostElapsed(const std::string tag);
double deviceElapsed(const std::string tag);
long long int count(const std::string tag);
double query(const std::string tag,std::string metric);
void printRunStat(int step);
void printStatEntry(std::string name, std::string tag, std::string type, double tNorm);
void printStatEntry(std::string name, double time, double tNorm);
void printStatEntry(std::string name, double tTag, long long int nCalls, double tNorm);

// print every entry in the map
void printAll();

// obtain all tags registered with the timer
std::vector<std::string> tags();
};
}

#endif

#if !defined(nekrs_nekrsascent_hpp_)
#define nekrs_nekrsascent_hpp_

/*
   insitu visualizatoin via Ascent 
*/

#include <vector>
#include <tuple>

#include "nrs.hpp"

#ifdef ENABLE_ASCENT
#include <ascent.hpp>
#include <mpi.h>

namespace nekrsAscent
{
typedef std::vector< std::tuple<std::string, occa::memory, dlong> > fields;
void setup(mesh_t *mesh_, const dlong fieldOffset_, const fields& flds);
void setup(nrs_t *nrs_);
void run(const double time, const int tstep);
void finalize();
occa::memory o_getAscentFields();
static ascent::Ascent mAscent;
}

void initializeAscent();
void printStat();

#else

namespace nekrsAscent
{
typedef std::vector< std::tuple<std::string, occa::memory, dlong> > fields;
void setup(mesh_t *mesh_, const dlong fieldOffset_, const fields& flds);
void setup(nrs_t *nrs_);
void run(const double time, const int tstep);
void finalize();
occa::memory o_getAscentFields();
}
#endif // ascent
#endif // hpp

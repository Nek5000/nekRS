#if !defined(nekrs_nekrsascent_hpp_)
#define nekrs_nekrsascent_hpp_

/*
   insitu visualizatoin via Ascent 
*/

#include <vector>
#include <tuple>

#include <ascent.hpp>
#include <mpi.h>

#include "nrs.hpp"


namespace nekrsAscent
{
typedef std::vector< std::tuple<std::string, occa::memory, dlong> > fields;
// TODO: add another setup for all active fields
void setup(mesh_t *mesh_, const dlong fieldOffset_, const fields& flds);
void run(const double time, const int tstep);
void printStat(const int tstep);
static ascent::Ascent mAscent;
occa::memory o_avg(); // TODO, add interface to directly access the data?
void finalize();
}

#endif 

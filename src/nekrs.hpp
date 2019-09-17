#if !defined(nekrs_nekrs_hpp_)
#define nekrs_nekrs_hpp_

#define NEKRS_VERSION "19"
#define NEKRS_SUBVERSION "0"

#define EXIT(a)  { MPI_Finalize(); exit(a); } 

#include "libParanumal.hpp"
#include "ins.h"

// std::to_string might be not accurate enough 
static string to_string_f(double a) {
  stringstream s;
  s << std::scientific << a;
  return s.str();
}

ins_t *setup(mesh_t *mesh, setupAide &options);
void runPlan4(ins_t *ins);
void restartRead(ins_t *ins, setupAide &options);
void report(ins_t *ins, dfloat time, int tstep);

libParanumal::setupAide parRead(std::string &setupFile, MPI_Comm comm); 

#endif

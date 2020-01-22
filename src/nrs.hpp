#if !defined(nekrs_nekrs_hpp_)
#define nekrs_nekrs_hpp_

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <getopt.h>
#include <mpi.h>

#define NEKRS_VERSION "20"
#define NEKRS_SUBVERSION "0"

#define EXIT(a)  { MPI_Finalize(); exit(a); } 

#include "libParanumal.hpp"
#include "ins.h"
#include "timer.hpp"

// std::to_string might be not accurate enough 
static string to_string_f(double a) {
  stringstream s;
  s << std::scientific << a;
  return s.str();
}

#endif

#if !defined(nekrs_nekrs_hpp_)
#define nekrs_nekrs_hpp_

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <getopt.h>

#include "mpi_wrapper.hpp"

#define NEKRS_VERSION "0"
#define NEKRS_SUBVERSION "2"

#define NEKLDIMT 2

#define EXIT(a)  { MPI_Finalize(); exit(a); }

#include "libParanumal.hpp"
#include "ins.h"

// std::to_string might be not accurate enough
static string to_string_f(double a) {
  stringstream s;
  s << std::scientific << a;
  return s.str();
}

extern int ciMode;

#endif

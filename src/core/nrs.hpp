#if !defined(nekrs_nekrs_hpp_)
#define nekrs_nekrs_hpp_

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <getopt.h>

#include "mesh.h"
#include "mesh3D.h"
#include "ins.h"
#include "timer.hpp"

occa::device occaDeviceConfig(setupAide &options, MPI_Comm comm);

// std::to_string might be not accurate enough
static string to_string_f(double a)
{
  stringstream s;
  s << std::scientific << a;
  return s.str();
}

static std::vector<std::string> serializeString(const std::string sin)
{
  std::vector<std::string> slist;
  string s(sin);
  s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
  std::stringstream ss;
  ss.str(s);
  while( ss.good() ) {
    std::string substr;
    std::getline(ss, substr, ',');
    slist.push_back(substr);
  }
  return slist;
}

#endif

#include <stdlib.h>
#include "nrs.hpp"

void printHeader()
{
  std::cout << R"(                 __    ____  _____)" << std::endl
            << R"(   ____   ___   / /__ / __ \/ ___/)" << std::endl
            << R"(  / __ \ / _ \ / //_// /_/ /\__ \ )" << std::endl
            << R"( / / / //  __// ,<  / _, _/___/ / )" << std::endl
            << R"(/_/ /_/ \___//_/|_|/_/ |_|/____/  )"
            << "v" << NEKRS_VERSION << "." << NEKRS_SUBVERSION << "." << NEKRS_PATCHVERSION 
	    << " (" GITCOMMITHASH << ")" << std::endl
            << std::endl
            << "COPYRIGHT (c) 2019-2023 UCHICAGO ARGONNE, LLC" << std::endl
            << std::endl;
}

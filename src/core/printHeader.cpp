#include "nekrsSys.hpp"

void printHeader()
{

std::cout << R"(███    ██ ███████ ██   ██ ██████  ███████)" << std::endl
          << R"(████   ██ ██      ██  ██  ██   ██ ██     )" << std::endl
          << R"(██ ██  ██ █████   █████   ██████  ███████)" << std::endl
          << R"(██  ██ ██ ██      ██  ██  ██   ██      ██)" << std::endl
          << R"(██   ████ ███████ ██   ██ ██   ██ ███████)" << std::endl
            << "(c) 2019-2024 UCHICAGO ARGONNE, LLC" << std::endl
            << "v" << NEKRS_VERSION << "." << NEKRS_SUBVERSION << "." << NEKRS_PATCHVERSION 
	    << " (sha:" << GITCOMMITHASH << ")" << std::endl
            << std::endl;
}

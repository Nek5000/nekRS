#if !defined(nekrs_environment_hpp_)
#define nekrs_environment_hpp_

#include <string>

#include "mpi.h"

#define stringify_(x) #x
#define stringify(x) stringify_(x)

namespace os{
  /* private const variables */
  const static std::string separator_(stringify(NEKRS_PATH_SEPARATOR));

  int exist(std::string fname);
  int readable(std::string fname);

  std::string joinPath(std::string root,std::string relative);
  std::string getWorkingDir();
  void makeDir(std::string dirName);
}

namespace env{
  /* private const variables */
  const static std::string installDir_(stringify(NEKRS_INSTALL_DIR));
  const static std::string delim_="::=";

  /* private variables */
  static MPI_Comm comm_;
  static std::map<std::string,std::string> m_;

  /* private functions */
  static int readConfig();

  int init(MPI_Comm comm);
  int set(const std::string envVar,const std::string value);
  std::string get(std::string envVar);

  std::string installDir();
  std::string binDir();
  std::string shareDir();
  std::string nekDir();
  std::string nekPpList();
  std::string libPDir();
  std::string libPDefines();
  std::string udfDir();
  std::string nekInterfaceDir();
  std::string occaDir();
  std::string cacheDir();
  void makeCacheDir();
  std::string occaCacheDir();
  void makeOccaCacheDir();

  std::string cxxCompiler();
  std::string cxxFlags();
  std::string cCompiler();
  std::string cFlags();
  std::string fortranCompiler();
  std::string fortranFlags();
}

#endif

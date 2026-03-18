#include <dlfcn.h>
#include <stdlib.h>
#include "unifdef.h"

#include "udf.hpp"
#include "fileUtils.hpp"
#include "fileBcast.hpp"
#include "platform.hpp"
#include "bdryBase.hpp"
#include "sha1.hpp"

#include "solver.hpp"

#include "udfMake.hpp"

UDF udf = {NULL, NULL, NULL, NULL};

static int dirichletConditions = 0;
static int neumannConditions = 0;
static int RobinConditions = 0;

static void *libudfHandle = nullptr;
static std::string udfFile;

static void verifyOudf()
{
  for (auto &[key, value] : platform->app->bc->bIdToTypeId()) {
    const auto field = key.first;
    const auto typeId = value;

    if (typeId == bdryBase::bcType_udfDirichlet || typeId == bdryBase::bcType_interpolation) {
      oudfFindDirichlet(field);
    }

    if (typeId == bdryBase::bcType_udfNeumann || typeId == bdryBase::bcType_zeroDirichletX_udfNeumann ||
        typeId == bdryBase::bcType_zeroDirichletY_udfNeumann ||
        typeId == bdryBase::bcType_zeroDirichletZ_udfNeumann ||
        typeId == bdryBase::bcType_zeroDirichletN_udfNeumann) {
      oudfFindNeumann(field);
    }

    if (typeId == bdryBase::bcType_udfRobin) {
      oudfFindRobin(field);
    }
  }
}

void oudfFindDirichlet(const std::string &field)
{
  nekrsCheck(!dirichletConditions,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s\n",
             "Cannot find required okl function udfDirichlet!");
}

void oudfFindNeumann(const std::string &field)
{
  nekrsCheck(!neumannConditions,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s\n",
             "Cannot find required okl function udfNeumann!");
}

void oudfFindRobin(const std::string &field)
{
  nekrsCheck(!RobinConditions,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s\n",
             "Cannot find required okl function udfRobin!");
}

void adjustOudf(bool buildRequired, const std::string &postOklSource, const std::string &filePath)
{
  std::stringstream buffer;
  {
    std::ifstream f;
    f.open(postOklSource);
    buffer << f.rdbuf();
    f.close();
  }

  std::ofstream f;
  f.open(filePath, std::ios_base::app);

  if (buildRequired) {
    f << "#ifdef __okl__\n";
  }

  dirichletConditions = std::regex_search(buffer.str(), std::regex(R"(\s*void\s+udfDirichlet)"));
  if (!dirichletConditions && buildRequired) {
    f << "void udfDirichlet(bcData *bc){}\n";
  }

  neumannConditions = std::regex_search(buffer.str(), std::regex(R"(\s*void\s+udfNeumann)"));
  if (!neumannConditions && buildRequired) {
    f << "void udfNeumann(bcData *bc){}\n";
  }

  RobinConditions = std::regex_search(buffer.str(), std::regex(R"(\s*void\s+udfRobin)"));
  if (!RobinConditions && buildRequired) {
    f << "void udfRobin(bcData *bc){}\n";
  }

  if (buildRequired) {
    f << "#endif\n";
  }

  f.close();
}

void udfBuild(setupAide &options)
{
  options.getArgs("UDF FILE", udfFile);

  udfFile = fs::absolute(udfFile);
  if (platform->comm.mpiRank() == 0) {
    nekrsCheck(!fs::exists(udfFile), MPI_COMM_SELF, EXIT_FAILURE, "Cannot find %s!\n", udfFile.c_str());
  }

  const int verbose = platform->verbose() ? 1 : 0;
  const std::string installDir(getenv("NEKRS_HOME"));
  const std::string cache_dir(getenv("NEKRS_CACHE_DIR"));
  const std::string udfLib = cache_dir + "/udf/libudf.so";
  const std::string udfFileCache = cache_dir + "/udf/udf.cpp";
  const std::string udfHashFile = cache_dir + "/udf/udf.hash";
  const std::string oudfFileCache = cache_dir + "/udf/udf.okl";
  const std::string case_dir(fs::current_path());
  const std::string casename = options.getArgs("CASENAME");

  const std::string cmakeBuildDir = cache_dir + "/udf";
  const std::string postOklSource = cmakeBuildDir + "/CMakeFiles/OKL.dir/okl.cpp.i";

  const std::string libnekrsFile = (sizeof(dfloat) == sizeof(float)) ? installDir + "/lib/libnekrs-fp32.so"
                                                                     : installDir + "/lib/libnekrs.so";
  const std::string libnekrsHashFile = cache_dir + "/udf/libnekrs.hash";

  std::string oudfFile;
  options.getArgs("UDF OKL FILE", oudfFile);
  oudfFile = fs::absolute(oudfFile);

  MPI_Comm comm = (platform->cacheLocal) ? platform->comm.mpiCommLocal() : platform->comm.mpiComm();
  int buildRank;
  MPI_Comm_rank(comm, &buildRank);

  int buildRequired = 0;
  if (platform->comm.mpiRank() == 0) {

    auto getHash = [&](const std::string &fname) {
      std::ifstream f(fname);
      if (!f.is_open()) {
        return std::string("");
      }
      std::stringstream buffer;
      buffer << f.rdbuf();
      f.close();

      return buffer.str();
    };

    // changes in udf include files + env-vars are currently not detected
    // note, we want to avoid calling system()
    if (options.compareArgs("BUILD ONLY", "TRUE")) {
      buildRequired = 1;
    } else if (!fs::exists(udfLib) || !fs::exists(oudfFileCache)) {
      buildRequired = 1;
    } else if (SHA1::from_file(udfFile) != getHash(udfHashFile)) {
      buildRequired = 1;
    } else if (SHA1::from_file(libnekrsFile) != getHash(libnekrsHashFile)) {
      buildRequired = 1;
    }

    if (fs::exists(std::string(case_dir + "/ci.inc"))) {
      if (isFileNewer(std::string(case_dir + "/ci.inc").c_str(), udfFileCache.c_str())) {
        buildRequired = 1;
      }
    }

    if (fs::exists(oudfFile)) {
      if (isFileNewer(oudfFile.c_str(), oudfFileCache.c_str())) {
        buildRequired = 1;
      }
    }

    // check for a typical include file
    if (fs::exists(std::string(case_dir + "/" + casename + ".okl"))) {
      if (isFileNewer(std::string(case_dir + "/" + casename + ".okl").c_str(), oudfFileCache.c_str())) {
        buildRequired = 1;
      }
    }
  }
  MPI_Bcast(&buildRequired, 1, MPI_INT, 0, comm);

  int oudfFileExists;
  if (platform->comm.mpiRank() == 0) {
    oudfFileExists = fs::exists(oudfFile);
  }
  MPI_Bcast(&oudfFileExists, 1, MPI_INT, 0, comm);
  if (!oudfFileExists) {
    options.removeArgs("UDF OKL FILE");
  }

  if (platform->cacheBcast || platform->cacheLocal) {
    if (oudfFileExists) {
      fileBcast(oudfFile, platform->tmpDir, comm, platform->verbose());
      oudfFile = platform->tmpDir / fs::path(oudfFile).filename();
      options.setArgs("UDF OKL FILE", std::string(oudfFile));
    }
  }

  {
    const auto err = (buildRank == 0 && buildRequired)
                         ? udfMake(options, platform->app->id(), platform->comm.mpiRank())
                         : 0;
    auto log = cmakeBuildDir + "/cmake.log";
    nekrsCheck(err, platform->comm.mpiComm(), EXIT_FAILURE, "see %s for more details\n", log.c_str());
  }

  if (buildRank == 0) {
    if (fs::exists(cache_dir + "/udf/okl.cpp")) {
      fs::rename(cache_dir + "/udf/okl.cpp", oudfFileCache);
    }

    adjustOudf(buildRequired, postOklSource, oudfFileCache); // call every time for verifyOudf
    verifyOudf();

    fileSync(oudfFileCache.c_str());
  }

  if (platform->cacheBcast || platform->cacheLocal) {
    const auto dst = fs::path(platform->tmpDir) / "udf";
    fileBcast(fs::path(udfLib), dst, comm, platform->verbose());
    fileBcast(fs::path(oudfFileCache), dst, comm, platform->verbose());
  }

  // some BC kernels will include this file
  options.setArgs("OKL FILE CACHE", oudfFileCache);
  if (platform->cacheBcast || platform->cacheLocal) {
    options.setArgs("OKL FILE CACHE", std::string(platform->tmpDir + "/udf/udf.okl"));
  }

  MPI_Bcast(&dirichletConditions, 1, MPI_INT, 0, comm);
  MPI_Bcast(&neumannConditions, 1, MPI_INT, 0, comm);
  MPI_Bcast(&RobinConditions, 1, MPI_INT, 0, comm);
}

void *udfLoadFunction(const char *fname, int errchk)
{
  if (!libudfHandle) {
    std::string cache_dir(getenv("NEKRS_CACHE_DIR"));
    if (platform->cacheBcast) {
      cache_dir = fs::path(platform->tmpDir);
    }

    const auto udfLib = std::string(fs::path(cache_dir) / "udf/libudf.so");

    if (platform->comm.mpiRank() == 0 && platform->verbose()) {
      std::cout << "loading " << udfLib << std::endl;
    }

    libudfHandle = dlopen(udfLib.c_str(), RTLD_NOW | RTLD_GLOBAL);
    nekrsCheck(!libudfHandle, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", dlerror());
  }

  void *fptr = dlsym(libudfHandle, fname);
  nekrsCheck(!fptr && errchk, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", dlerror());

  dlerror();

  return fptr;
}

void udfUnload()
{
  dlclose(libudfHandle);
}

void udfLoad()
{
  *(void **)(&udf.setup0) = udfLoadFunction("UDF_Setup0", 0);
  *(void **)(&udf.setup) = udfLoadFunction("UDF_Setup", 1);
  *(void **)(&udf.loadKernels) = udfLoadFunction("UDF_LoadKernels", 0);
  *(void **)(&udf.autoloadKernels) = udfLoadFunction("UDF_AutoLoadKernels", 0);
  *(void **)(&udf.autoloadPlugins) = udfLoadFunction("UDF_AutoLoadPlugins", 1);
  *(void **)(&udf.executeStep) = udfLoadFunction("UDF_ExecuteStep", 1);
  *(void **)(&udf.finalize) = udfLoadFunction("UDF_Finalize", 0);
}

void udfEcho()
{
  const std::string cache_dir(getenv("NEKRS_CACHE_DIR"));
  const std::string oudfFileCache = cache_dir + "/udf/udf.okl";

  const auto tmpFile = udfFile + ".unifdef";
  unifdef("__okl__", udfFile.c_str(), tmpFile.c_str());

  std::ifstream fudf(tmpFile);
  std::string text;
  while (std::getline(fudf, text)) {
    std::cout << "<<< " << text << "\n";
  }
  std::cout << std::endl;
  fs::remove(tmpFile);

  std::ifstream foudf(oudfFileCache);
  while (std::getline(foudf, text)) {
    std::cout << "<<< " << text << "\n";
  }
  std::cout << std::endl;
}

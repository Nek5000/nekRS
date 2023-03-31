#include <iostream>
#include <unistd.h>
#include <dlfcn.h>
#include <stdlib.h>
#include <regex>

#include "udf.hpp"
#include "fileUtils.hpp"
#include "platform.hpp"
#include "bcMap.hpp"
#include "unifdef.h"

UDF udf = {NULL, NULL, NULL, NULL};

static int velocityDirichletConditions = 0;
static int meshVelocityDirichletConditions = 0;
static int velocityNeumannConditions = 0;
static int pressureDirichletConditions = 0;
static int scalarDirichletConditions = 0;
static int scalarNeumannConditions = 0;

void oudfFindDirichlet(std::string &field)
{
  nrsCheck(field.find("velocity") != std::string::npos && !velocityDirichletConditions,
           platform->comm.mpiComm, EXIT_FAILURE,
           "Cannot find velocityDirichletConditions!\n", "");

  nrsCheck(field.find("scalar") != std::string::npos && !scalarDirichletConditions,
           platform->comm.mpiComm, EXIT_FAILURE,
           "Cannot find scalarDirichletConditions!\n", "");

  if (field == "pressure" && !pressureDirichletConditions) {
    if (platform->comm.mpiRank == 0)
      std::cout << "WARNING: Cannot find pressureDirichletConditions!\n";
  }
  nrsCheck(field.find("mesh") != std::string::npos && !meshVelocityDirichletConditions &&
           !bcMap::useDerivedMeshBoundaryConditions(), platform->comm.mpiComm, EXIT_FAILURE,
           "Cannot find meshVelocityDirichletConditions!\n", "");
}

void oudfFindNeumann(std::string &field)
{
  nrsCheck(field.find("velocity") != std::string::npos && !velocityNeumannConditions,
           platform->comm.mpiComm, EXIT_FAILURE,
           "Cannot find velocityNeumannConditions!\n", "");
  nrsCheck(field.find("scalar") != std::string::npos && !scalarNeumannConditions,
           platform->comm.mpiComm, EXIT_FAILURE,
           "Cannot find scalarNeumannConditions!\n", "");
}

void adjustOudf(const std::string &filePath)
{
  std::fstream df;
  df.open(filePath, std::fstream::out | std::fstream::in);

  std::stringstream buffer;
  buffer << df.rdbuf();

  bool found = std::regex_search(buffer.str(), std::regex(R"(\s*void\s+velocityDirichletConditions)"));
  velocityDirichletConditions = found;
  if (!found)
    df << "void velocityDirichletConditions(bcData *bc){}\n";

  found = std::regex_search(buffer.str(), std::regex(R"(\s*void\s+meshVelocityDirichletConditions)"));
  meshVelocityDirichletConditions = found;
  if (!found)
    df << "void meshVelocityDirichletConditions(bcData *bc){\n"
          "  velocityDirichletConditions(bc);\n"
          "}\n";

  found = std::regex_search(buffer.str(), std::regex(R"(\s*void\s+velocityNeumannConditions)"));
  velocityNeumannConditions = found;
  if (!found)
    df << "void velocityNeumannConditions(bcData *bc){}\n";

  found = std::regex_search(buffer.str(), std::regex(R"(\s*void\s+pressureDirichletConditions)"));
  pressureDirichletConditions = found;
  if (!found)
    df << "void pressureDirichletConditions(bcData *bc){}\n";

  found = std::regex_search(buffer.str(), std::regex(R"(\s*void\s+scalarNeumannConditions)"));
  scalarNeumannConditions = found;
  if (!found)
    df << "void scalarNeumannConditions(bcData *bc){}\n";

  found = std::regex_search(buffer.str(), std::regex(R"(\s*void\s+scalarDirichletConditions)"));
  scalarDirichletConditions = found;

  if (!found)
    df << "void scalarDirichletConditions(bcData *bc){}\n";

  df.close();
}

void udfAutoKernels(const std::string& udfFileCache, const std::string& oudfFileCache)
{
  const std::string includeFile = fs::path(udfFileCache).parent_path() / fs::path("udfAutoLoadKernel.hpp");

  std::ifstream stream(oudfFileCache);
  std::regex exp(R"(\s*@kernel\s*void\s*([\S]*)\s*\()");
  std::smatch res;
  std::ofstream f(includeFile);

  std::string line;
  std::vector<std::string> kernelNameList;
  while (std::getline(stream, line)) {
    if (std::regex_search(line, res, exp)) {
      std::string kernelName = res[1];
      kernelNameList.push_back(kernelName);
      f << "static occa::kernel " << kernelName << ";\n";
    }
  }

  f << "void UDF_AutoLoadKernels(occa::properties& kernelInfo)" << std::endl << "{" << std::endl;

  for (auto entry : kernelNameList) {
    f << "  " << entry << " = "
      << "oudfBuildKernel(kernelInfo, \"" << entry << "\");" << std::endl;
  }

  f << "}" << std::endl;

  f.close();
  fileSync(includeFile.c_str());
  stream.close();
}

void udfBuild(const char *_udfFile, setupAide &options)
{
  std::string udfFile = fs::absolute(_udfFile);
  if (platform->comm.mpiRank == 0)
    nrsCheck(!fs::exists(udfFile), MPI_COMM_SELF, EXIT_FAILURE, "Cannot find %s!\n", udfFile.c_str()); 

  const int verbose = options.compareArgs("VERBOSE", "TRUE") ? 1 : 0;
  const std::string installDir(getenv("NEKRS_HOME"));
  const std::string udf_dir = installDir + "/udf";
  const std::string cache_dir(getenv("NEKRS_CACHE_DIR"));
  const std::string udfLib = cache_dir + "/udf/libUDF.so";
  const std::string udfFileCache = cache_dir + "/udf/udf.cpp";
  const std::string oudfFileCache = cache_dir + "/udf/udf.okl";
  const std::string case_dir(fs::current_path());
  const std::string casename = options.getArgs("CASENAME");

  std::string oudfFile;
  options.getArgs("UDF OKL FILE", oudfFile);
  oudfFile = fs::absolute(oudfFile);


  MPI_Comm comm = (platform->cacheLocal) ? platform->comm.mpiCommLocal : platform->comm.mpiComm;
  int buildRank;
  MPI_Comm_rank(comm, &buildRank);

  int buildRequired = 0;
  if(platform->comm.mpiRank == 0) {

    // there is no mechansism to check for dependency changes
    // without invoking system, comparing file timestamps for now 
    // note, this will not work for unknown include files and env-var changes! 
    if (platform->options.compareArgs("BUILD ONLY", "TRUE")) {
      buildRequired = 1;
    } else if (!fs::exists(udfLib)) {
      buildRequired = 1;
    } else if (isFileNewer(udfFile.c_str(), udfFileCache.c_str())) {
      buildRequired = 1;
    };

    if (fs::exists(oudfFile)) {
      if (isFileNewer(oudfFile.c_str(), oudfFileCache.c_str())) 
        buildRequired = 1;
    }

    if (fs::exists(std::string(case_dir + "/ci.inc"))) {
      if (isFileNewer(std::string(case_dir + "/ci.inc").c_str(), udfFileCache.c_str()))
        buildRequired = 1;
    }

    if (fs::exists(std::string(case_dir + "/" + casename + ".okl"))) {
      if (isFileNewer(std::string(case_dir + "/" + casename + ".okl").c_str(), oudfFileCache.c_str())) 
        buildRequired = 1;
    }

  }
  MPI_Bcast(&buildRequired, 1, MPI_INT, 0, comm);

  int oudfFileExists;
  if(platform->comm.mpiRank == 0) oudfFileExists = fs::exists(oudfFile);
  MPI_Bcast(&oudfFileExists, 1, MPI_INT, 0, comm);
  if(!oudfFileExists) options.removeArgs("UDF OKL FILE");

  if (platform->cacheBcast || platform->cacheLocal) {
    if (oudfFileExists) {
      fileBcast(oudfFile, platform->tmpDir, comm, platform->verbose);
      oudfFile = platform->tmpDir / fs::path(oudfFile).filename();
      options.setArgs("UDF OKL FILE", std::string(oudfFile));
    }
  }

  if (platform->cacheBcast) {
    fileBcast(udfFile, platform->tmpDir, comm, platform->verbose);
    udfFile = std::string(platform->tmpDir / fs::path(udfFile).filename());
    options.setArgs("UDF FILE", udfFile);
  }

  int err = 0;
  err += [&]() {
    if (buildRank == 0) {
      double tStart = MPI_Wtime();

      char cmd[4096];
      mkdir(std::string(cache_dir + "/udf").c_str(), S_IRWXU);

      const std::string pipeToNull =
          (platform->comm.mpiRank == 0) ? std::string("") : std::string("> /dev/null 2>&1");

      if (buildRequired) {
        if (platform->comm.mpiRank == 0)
          printf("building udf ... \n");
        fflush(stdout);

        copyFile(udfFile.c_str(), udfFileCache.c_str());

        // generate udfFileCache
        {
          std::ofstream f(udfFileCache, std::ios::trunc);
          f << "#include \"" << udfFile << "\"" << std::endl;

          // autoload plugins
          std::map<std::string, std::string> pluginTable =
          {
            {"nekrs_tavg_hpp_"        , "tavg::buildKernel"},
            {"nekrs_RANSktau_hpp_"    , "RANSktau::buildKernel"},
            {"nekrs_lowMach_hpp_"     , "lowMach::buildKernel"},
            {"nekrs_velRecycling_hpp_", "velRecycling::buildKernel"},
            {"nekrs_lpm_hpp_"         , "lpm_t::registerKernels"}
          };

          f << "void UDF_AutoLoadPlugins(occa::properties& kernelInfo)" << std::endl
            << "{" << std::endl;

            for (auto const& plugin : pluginTable) {
              f << "#ifdef " << plugin.first << std::endl
                << "  " << plugin.second << "(kernelInfo);" << std::endl
                << "#endif" << std::endl;
            }

          f << "}" << std::endl; 
          f.close();
        }

        copyFile(std::string(udf_dir + std::string("/CMakeLists.txt")).c_str(),
                 std::string(cache_dir + std::string("/udf/CMakeLists.txt")).c_str());

        std::string cmakeFlags("-Wno-dev");
        if (verbose)
          cmakeFlags += " --trace-expand";
        std::string cmakeBuildDir = cache_dir + "/udf";

        { // generate dummy to make cmake happy that the file exists
          const std::string includeFile = std::string(cache_dir + std::string("/udf/udfAutoLoadKernel.hpp"));
          std::ofstream f(includeFile);
          f << "// dummy";
          f.close();
          fileSync(includeFile.c_str());
        }

        extract_ifdef("__okl__", udfFile.c_str(), std::string(cache_dir + "/udf/okl.cpp").c_str());
        bool oklSectionFound = fs::file_size(cache_dir + "/udf/okl.cpp");

        if(!oklSectionFound && oudfFileExists) {
           copyFile(oudfFile.c_str(), std::string(cache_dir + "/udf/okl.cpp").c_str());
           if(platform->comm.mpiRank == 0)
             printf("Cannot find okl section in udf (oudf will be deprecated in next version!)\n");
        } else if (!oklSectionFound) {
          if(platform->comm.mpiRank == 0)
            printf("Cannot find oudf or okl section in udf\n");
          return EXIT_FAILURE;
        }

        sprintf(cmd,
                "rm -f %s/*.so && cmake %s -S %s -B %s "
                "-DNEKRS_INSTALL_DIR=\"%s\" -DCASE_DIR=\"%s\" -DCMAKE_CXX_COMPILER=\"$NEKRS_CXX\" "
                "-DCMAKE_CXX_FLAGS=\"$NEKRS_CXXFLAGS\" %s >cmake.log 2>&1",
                cmakeBuildDir.c_str(),
                cmakeFlags.c_str(),
                cmakeBuildDir.c_str(),
                cmakeBuildDir.c_str(),
                installDir.c_str(),
                case_dir.c_str(),
                pipeToNull.c_str());
        const int retVal = system(cmd);
        if (verbose && platform->comm.mpiRank == 0)
          printf("%s (cmake retVal: %d)\n", cmd, retVal);
        if (retVal)
          return EXIT_FAILURE;

        { // generate pre-processed okl
          sprintf(cmd, "cd %s && make -j1 okl.i %s", cmakeBuildDir.c_str(), pipeToNull.c_str());
          const std::string postOklSource = cmakeBuildDir + "/CMakeFiles/OKL.dir/okl.cpp.i";
          const int retVal = system(cmd);
          if (verbose && platform->comm.mpiRank == 0)
            printf("%s (preprocessing retVal: %d)\n", cmd, retVal);
          if (retVal)
            return EXIT_FAILURE;
          
          std::smatch match;
          std::ifstream fi(postOklSource);
          std::ofstream fo(cache_dir + "/udf/okl.cpp", std::ios::trunc);
          std::stringstream buffer;
          buffer << fi.rdbuf();

          // filter line info as not supported by occa's parser
          fo << std::regex_replace(buffer.str(), std::regex(R"(#\s*\d.*\n)"), "");

          fi.close();
          fo.close();
        }

        if(oklSectionFound)
          udfAutoKernels(udfFileCache, cache_dir + "/udf/okl.cpp");

        { // build
          sprintf(cmd, "cd %s/udf && make -j1 %s", cache_dir.c_str(), pipeToNull.c_str());
          const int retVal = system(cmd);
          if (verbose && platform->comm.mpiRank == 0) {
            printf("%s (make retVal: %d)\n", cmd, retVal);
          }
          if (retVal)
            return EXIT_FAILURE;
          fileSync(udfLib.c_str());
        }

        if (platform->comm.mpiRank == 0)
          printf("done (%gs)\n", MPI_Wtime() - tStart);
        fflush(stdout);

      } // buildRequired

      if(fs::exists(cache_dir + "/udf/okl.cpp")) 
        fs::rename(cache_dir + "/udf/okl.cpp", oudfFileCache);

      adjustOudf(oudfFileCache); // called every time to check if required device BC functions exist

      fileSync(oudfFileCache.c_str());

    } // rank 0

    if (platform->cacheBcast || platform->cacheLocal)
      fileBcast(fs::path(udfFileCache).parent_path(), platform->tmpDir, comm, platform->verbose);

    return 0;
  }();

  options.setArgs("OKL FILE CACHE", oudfFileCache);
  if (platform->cacheBcast || platform->cacheLocal)
    options.setArgs("OKL FILE CACHE", std::string(platform->tmpDir + "/udf/udf.okl"));

  MPI_Bcast(&velocityDirichletConditions, 1, MPI_INT, 0, comm);
  MPI_Bcast(&meshVelocityDirichletConditions, 1, MPI_INT, 0, comm);
  MPI_Bcast(&velocityNeumannConditions, 1, MPI_INT, 0, comm);
  MPI_Bcast(&pressureDirichletConditions, 1, MPI_INT, 0, comm);
  MPI_Bcast(&scalarNeumannConditions, 1, MPI_INT, 0, comm);
  MPI_Bcast(&scalarDirichletConditions, 1, MPI_INT, 0, comm);

  MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_SUM, platform->comm.mpiComm);

  nrsCheck(err, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "see above and cmake.log for more details");
}

void *udfLoadFunction(const char *fname, int errchk)
{
  static void *h = nullptr;
  if(!h) {
    std::string cache_dir(getenv("NEKRS_CACHE_DIR"));
    if (platform->cacheBcast)
      cache_dir = fs::path(platform->tmpDir);

    const auto udfLib = std::string(fs::path(cache_dir) / "udf/libUDF.so");

    if(platform->comm.mpiRank == 0 && platform->verbose)
      std::cout << "loading " << udfLib << std::endl;

    h = dlopen(udfLib.c_str(), RTLD_NOW | RTLD_GLOBAL);
    nrsCheck(!h, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", dlerror());
  }

  void *fptr = dlsym(h, fname);
  nrsCheck(!fptr && errchk, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", dlerror());

  dlerror();

  return fptr;
}

void udfLoad()
{
  *(void **)(&udf.setup0) = udfLoadFunction("UDF_Setup0", 0);
  *(void **)(&udf.setup) = udfLoadFunction("UDF_Setup", 0);
  *(void **)(&udf.loadKernels) = udfLoadFunction("UDF_LoadKernels", 0);
  *(void **)(&udf.autoloadKernels) = udfLoadFunction("UDF_AutoLoadKernels", 0);
  *(void **)(&udf.autoloadPlugins) = udfLoadFunction("UDF_AutoLoadPlugins", 1);
  *(void **)(&udf.executeStep) = udfLoadFunction("UDF_ExecuteStep", 0);
}

occa::kernel oudfBuildKernel(occa::properties kernelInfo, const char *function)
{
  std::string installDir;
  installDir.assign(getenv("NEKRS_HOME"));
  const std::string bcDataFile = installDir + "/include/bdry/bcData.h";
  kernelInfo["includes"] += bcDataFile.c_str();

  // provide some common kernel args
  int N;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  const int Nq = N + 1;
  const int Np = Nq * Nq * Nq;
  kernelInfo["defines/p_Nq"] = Nq;
  kernelInfo["defines/p_Np"] = Np;

  std::string oudfCache;
  platform->options.getArgs("OKL FILE CACHE", oudfCache);

  return platform->device.buildKernel(oudfCache.c_str(), function, kernelInfo);
}

#include <iostream>
#include <unistd.h>
#include <dlfcn.h>
#include <stdlib.h>
#include <regex>
#include <set>
#include <chrono>
#include <fcntl.h>

#include "unifdef.h"
#include "sha1.hpp"

#include "setupAide.hpp"

static void udfAutoKernels(const std::string &udfFileCache,
                           const std::string &postOklSource,
                           const std::string &oudfFileCache)
{
  const std::string includeFile = fs::path(udfFileCache).parent_path() / fs::path("udfAutoLoadKernel.hpp");
  std::ofstream f(includeFile);

  auto buffer = [&]() {
    std::stringstream fbuffer;
    std::ifstream postOklSourceStream(postOklSource);
    fbuffer << postOklSourceStream.rdbuf();
    postOklSourceStream.close();
    return fbuffer.str();
  }();

  std::set<std::string> kernelNameList;
  std::regex rexp(R"(\s*@kernel\s*void\s*([\S]*)\s*\()");
  std::smatch res;

  {
    std::regex_token_iterator<std::string::iterator> end; // default constructor = end-of-sequence:
    std::regex_token_iterator<std::string::iterator> token(buffer.begin(), buffer.end(), rexp, 1);
    while (token != end) {
      const auto kernelName = *token++;
      kernelNameList.insert(kernelName);
    }
  }

  for (auto entry : kernelNameList) {
    f << "static occa::kernel " << entry << ";\n";
  }

  f << "void UDF_AutoLoadKernels(occa::properties& kernelInfo)" << std::endl << "{" << std::endl;

  for (auto entry : kernelNameList) {
    f << "  " << entry << " = "
      << "oudfBuildKernel(kernelInfo, \"" << entry << "\");" << std::endl;
  }

  f << "}" << std::endl;

  f.close();
  fileSync(includeFile.c_str());
}

int udfMake(setupAide &options, const std::string &solverName, int rank)
{
  std::string udfFile;
  options.getArgs("UDF FILE", udfFile);
  udfFile = fs::absolute(udfFile);

  const int verbose = platform->verbose() ? 1 : 0;
  const std::string installDir(getenv("NEKRS_HOME"));
  const std::string udf_dir = installDir + "/udf";
  const std::string cache_dir(getenv("NEKRS_CACHE_DIR"));
  const std::string udfLib = cache_dir + "/udf/libudf.so";
  const std::string udfFileCache = cache_dir + "/udf/udf.cpp";
  const std::string udfHashFile = cache_dir + "/udf/udf.hash";
  const std::string case_dir(fs::current_path());
  const std::string casename = options.getArgs("CASENAME");
  const std::string cmakeBuildDir = cache_dir + "/udf";
  const std::string postOklSource = cmakeBuildDir + "/CMakeFiles/OKL.dir/okl.cpp.i";
  const std::string libnekrsFile = (sizeof(dfloat) == sizeof(float)) ? installDir + "/lib/libnekrs-fp32.so"
                                                                     : installDir + "/lib/libnekrs.so";
  const std::string libnekrsHashFile = cache_dir + "/udf/libnekrs.hash";

  std::string oudfFile;
  int oudfFileExists = 0;
  options.getArgs("UDF OKL FILE", oudfFile);
  if (!oudfFile.empty()) {
    oudfFileExists = 1;
    oudfFile = fs::absolute(oudfFile);
  }

  auto tStart = std::chrono::high_resolution_clock::now();

  const int cmdSize = 4096;
  char cmd[cmdSize];

  fs::create_directories(std::string(cache_dir + "/udf"));

  const std::string pipeToNull = (rank == 0) ? std::string("") : std::string("> /dev/null 2>&1");

  if (rank == 0) {
    printf("building udf ... \n");
  }
  fflush(stdout);

  {
    std::ofstream f(udfHashFile, std::ios::trunc);
    f << SHA1::from_file(udfFile);
    f.close();
  }

  {
    std::ofstream f(libnekrsHashFile, std::ios::trunc);
    f << SHA1::from_file(libnekrsFile);
    f.close();
  }

  auto solverIncludes = [&]() {
    std::string txt;
    if (solverName == "nrs") {
      txt = "#include \"nrs.hpp\"";
      txt += "\n";
      txt += "const auto nrs = dynamic_cast<nrs_t*>(platform->app);";
      txt += "\n";
    }

    return txt;
  };

  // generate udfFileCache
  {
    std::ofstream f(udfFileCache, std::ios::trunc);
    f << "#include \"udf.hpp\"" << std::endl
      << "#include \"udfAutoLoadKernel.hpp\"" << std::endl

      << solverIncludes()

      << "#include \"udfHelper.hpp\"" << std::endl
      << "#include \"ci.hpp\"" << std::endl
      << "#include \"" << udfFile << "\"" << std::endl;

    // autoload plugins
    std::map<std::string, std::string> pluginTable = {{"nekrs_tavg_hpp_", "tavg::registerKernels"},
                                                      {"nekrs_RANSktau_hpp_", "RANSktau::buildKernel"},
                                                      {"nekrs_lowMach_hpp_", "lowMach::buildKernel"},
                                                      {"nekrs_recycling_hpp_", "planarCopy::buildKernel"},
                                                      {"nekrs_lpm_hpp_", "lpm_t::registerKernels"},
                                                      {"nekrs_Radiation_hpp_", "Radiation::buildKernel"}};

    f << "void UDF_AutoLoadPlugins(occa::properties& kernelInfo)" << std::endl << "{" << std::endl;

    for (auto const &plugin : pluginTable) {
      f << "#ifdef " << plugin.first << std::endl
        << "  " << plugin.second << "(kernelInfo);" << std::endl
        << "#endif" << std::endl;
    }

    f << "}" << std::endl;
    f.close();
  }

  fs::copy(std::string(udf_dir + "/CMakeLists.txt"),
           std::string(cache_dir + "/udf/CMakeLists.txt"),
           fs::copy_options::overwrite_existing);

  std::string cmakeFlags("-Wno-dev");
  if (verbose) {
    cmakeFlags += " --trace-expand";
  }

  { // generate dummy to make cmake happy that the file exists
    const std::string includeFile = std::string(cache_dir + std::string("/udf/udfAutoLoadKernel.hpp"));
    std::ofstream f(includeFile);
    f << "// dummy";
    f.close();
    fileSync(includeFile.c_str());
  }

  extract_ifdef("__okl__", udfFile.c_str(), std::string(cache_dir + "/udf/okl.cpp").c_str());
  bool oklSectionFound = fs::file_size(cache_dir + "/udf/okl.cpp");

  if (!oklSectionFound && oudfFileExists) {
    fs::copy(oudfFile, std::string(cache_dir + "/udf/okl.cpp"), fs::copy_options::overwrite_existing);
    if (rank == 0) {
      printf("Cannot find okl section in udf (oudf will be deprecated in next version!)\n");
    }
  } else if (!oklSectionFound) {
    if (rank == 0) {
      printf("Cannot find oudf or okl section in udf\n");
    }
    return EXIT_FAILURE;
  }

  const std::string useFloat = (sizeof(dfloat) == sizeof(float)) ? "ON" : "OFF";
  const std::string cmakeVerbose = (verbose) ? "ON" : "OFF";

  snprintf(cmd,
           cmdSize,
           "cmake %s -S %s -B %s "
           "-DNEKRS_USE_DFLOAT_FLOAT=%s "
           "-DNEKRS_INSTALL_DIR=\"%s\" -DCASE_DIR=\"%s\" -DCMAKE_CXX_COMPILER=\"$NEKRS_CXX\" "
           "-DCMAKE_CXX_FLAGS=\"$NEKRS_CXXFLAGS\" -DCMAKE_VERBOSE_MAKEFILE=%s >cmake.log 2>&1",
           cmakeFlags.c_str(),
           cmakeBuildDir.c_str(),
           cmakeBuildDir.c_str(),
           useFloat.c_str(),
           installDir.c_str(),
           case_dir.c_str(),
           cmakeVerbose.c_str());

  fs::remove(cmakeBuildDir + "/libudf.so");

  const int retVal = system(cmd);
  if (verbose && rank == 0) {
    printf("%s (cmake retVal: %d)\n", cmd, retVal);
  }
  if (retVal) {
    return EXIT_FAILURE;
  }

  auto stdoutFlag = (verbose) ? std::string("") : ">>cmake.log 2>&1";

  { // generate pre-processed okl
    snprintf(cmd, cmdSize, "cd %s && make -j1 okl.i %s", cmakeBuildDir.c_str(), stdoutFlag.c_str());
    const int retVal = system(cmd);
    if (verbose && rank == 0) {
      printf("%s (preprocessing retVal: %d)\n", cmd, retVal);
    }
    if (retVal) {
      return EXIT_FAILURE;
    }
  }

  if (oklSectionFound) {
    udfAutoKernels(udfFileCache, postOklSource, cache_dir + "/udf/okl.cpp");
  }

  { // build
    snprintf(cmd, cmdSize, "cd %s/udf && make -j1", cache_dir.c_str());
    const int retVal = system(cmd);
    if (verbose && rank == 0) {
      printf("%s (make retVal: %d)\n", cmd, retVal);
    }
    if (retVal) {
      return EXIT_FAILURE;
    }
    fileSync(udfLib.c_str());
  }

  if (rank == 0) {
    auto tEnd = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart);
    std::cout << "done (" << elapsed.count() << "s)" << std::endl;
  }
  fflush(stdout);

  return 0;
}

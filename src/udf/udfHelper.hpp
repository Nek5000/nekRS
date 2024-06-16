#if !defined(nekrs_udfhelp_hpp_)
#define nekrs_udfhelp_hpp_

static int failCnt = 0;
static std::vector <std::tuple<std::string, bool> > ciTestList;

static void CiPassTest(const std::string& txt)
{
  if (platform->comm.mpiRank == 0)
    printf("CI test <%s> passed\n", txt.c_str());
    if(!failCnt) platform->exitValue = 0;
}

static void CiFailTest(const std::string& txt, const std::string suffix = "")
{
  if (platform->comm.mpiRank == 0)
    printf("CI test <%s> failed %s\n", txt.c_str(), suffix.c_str());
    failCnt++;
    if(failCnt) platform->exitValue = 1;
}

static void CiPassTest() { CiPassTest(std::string("")); }
static void CiFailTest() { CiFailTest(std::string("")); }

static void CiEvalTest(const std::string& txt, bool passed, std::string suffix = "")
{
  ciTestList.push_back({txt, passed});

  if (passed) 
    CiPassTest(txt);
  else
    CiFailTest(txt, suffix);
}

static void CiAddTest(const std::string& txt, bool passed)
{
  ciTestList.push_back({txt, passed});
}

static void CiFinalize()
{
  if (platform->comm.mpiRank == 0) {
    std::cout << "CI test summary:\n"; 
    for(auto &entry : ciTestList) {
      std::cout << std::get<0>(entry) << ":" << std::get<0>(entry) << std::endl; 
    }

    std::cout << std::endl;
  }
}

// aliases
static const auto& _linAlg = platform->linAlg;
static const auto& _par = platform->par;
static const auto& _options = platform->options;
static const auto& _comm = platform->comm;

#endif

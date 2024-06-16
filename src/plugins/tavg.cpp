#include "nekInterfaceAdapter.hpp"
#include "tavg.hpp"
#include "platform.hpp"
#include "linAlg.hpp"
#include "fld.hpp"

// private members
namespace
{
ogs_t *ogs;

dlong fieldOffset;

occa::memory o_Uavg, o_Urms;
occa::memory o_Urm2;
occa::memory o_Pavg, o_Prms;
occa::memory o_Savg, o_Srms;

std::vector< std::vector<deviceMemory<dfloat>> > userFieldList;
occa::memory o_AVG;

occa::kernel EXKernel;
occa::kernel EXYKernel;
occa::kernel EXYZKernel;
occa::kernel E4Kernel;

bool buildKernelCalled = false;
bool setupCalled = false;

int counter = 0;

double atime;
double timel;

int outfldCounter = 0;
} // namespace

static void EX(dlong N, dfloat a, dfloat b, int nflds, occa::memory o_x, occa::memory o_EX)
{
  EXKernel(N, fieldOffset, nflds, a, b, o_x, o_EX);
}

static void EX(dlong N, dlong fieldOffset, dfloat a, dfloat b, int nflds, occa::memory o_x, occa::memory o_EX)
{
  EXKernel(N, fieldOffset, nflds, a, b, o_x, o_EX);
}

static void
EXY(dlong N, dfloat a, dfloat b, int nflds, occa::memory o_x, occa::memory o_y, occa::memory o_EXY)
{
  EXYKernel(N, fieldOffset, nflds, a, b, o_x, o_y, o_EXY);
}

static void EXYZ(dlong N,
                 dfloat a,
                 dfloat b,
                 int nflds,
                 occa::memory o_x,
                 occa::memory o_y,
                 occa::memory o_z,
                 occa::memory &o_EXYZ)
{
  EXYZKernel(N, fieldOffset, nflds, a, b, o_x, o_y, o_z, o_EXYZ);
}

static void E4(dlong N,
               dfloat a,
               dfloat b,
               int nflds,
               occa::memory o_1,
               occa::memory o_2,
               occa::memory o_3,
               occa::memory o_4,
               occa::memory &o_E4)
{
  E4Kernel(N, fieldOffset, nflds, a, b, o_1, o_2, o_3, o_4, o_E4);
}

void tavg::buildKernel(occa::properties kernelInfo)
{
  auto buildKernel = [&kernelInfo](const std::string& kernelName)
  {
    const auto path = getenv("NEKRS_KERNEL_DIR") + std::string("/plugins/");
    const auto fileName = path + "E.okl";
    const auto reqName = "tavg::";
    if (platform->options.compareArgs("REGISTER ONLY", "TRUE")) {
      platform->kernelRequests.add(reqName, fileName, kernelInfo);
      return occa::kernel();
    } else {
      buildKernelCalled = 1;
      return platform->kernelRequests.load(reqName, kernelName);
    }
  };

  EXKernel = buildKernel("EX");
  EXYKernel = buildKernel("EXY");
  EXYZKernel = buildKernel("EXYZ");
  E4Kernel = buildKernel("E4");
}

void tavg::reset()
{
  counter = 0;
  atime = 0;
}

void tavg::run(double time)
{
  nekrsCheck(!setupCalled || !buildKernelCalled,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s\n",
             "called prior to tavg::setup()!");

  if (!counter) {
    atime = 0;
    timel = time;
  }
  counter++;

  const double dtime = time - timel;
  atime += dtime;

  if (atime == 0 || dtime == 0) {
    return;
  }

  const dfloat b = dtime / atime;
  const dfloat a = 1 - b;

  if (userFieldList.size()) {
    int cnt = 0;
    for (auto &entry : userFieldList) {
      auto o_avg = o_AVG.slice(cnt * fieldOffset, fieldOffset);
      const auto N = fieldOffset;

      if (entry.size() == 1) {
        EX(N, a, b, 1, entry.at(0), o_avg);
      } else if (entry.size() == 2) {
        EXY(N, a, b, 1, entry.at(0), entry.at(1), o_avg);
      } else if (entry.size() == 3) {
        EXYZ(N, a, b, 1, entry.at(0), entry.at(1), entry.at(2), o_avg);
      } else if (entry.size() == 4) {
        E4(N, a, b, 1, entry.at(0), entry.at(1), entry.at(2), entry.at(3), o_avg);
      }
      cnt++;
    }
  }

  timel = time;
}

void tavg::setup(dlong _fieldOffset, const std::vector< std::vector<deviceMemory<dfloat>> >& flds)
{
  static bool isInitialized = false;
  if (isInitialized) {
    return;
  }
  isInitialized = true;

  nekrsCheck(!buildKernelCalled, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "called prior tavg::buildKernel()!");

  userFieldList = flds;

  for (auto &entry : userFieldList) {
    nekrsCheck(entry.size() < 1 || entry.size() > 4,
               platform->comm.mpiComm,
               EXIT_FAILURE,
               "%s\n",
               "invalid number of vectors in one of the user list entries!");

    for (auto &entry_i : entry) {
      nekrsCheck(entry_i.length() < _fieldOffset,
                 platform->comm.mpiComm,
                 EXIT_FAILURE,
                 "%s\n",
                 "vector size in one of the user list entries smaller than fieldOffset");
    }
  }

  fieldOffset = _fieldOffset;
  o_AVG = platform->device.malloc<dfloat>(userFieldList.size() * fieldOffset);

  setupCalled = true;
}

void tavg::outfld(int _outXYZ, int FP64)
{
  nekrsCheck(!setupCalled || !buildKernelCalled,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s\n",
             "called prior to tavg::setup()!");

  int outXYZ = _outXYZ;
  if (!outfldCounter) {
    outXYZ = 1;
  }

  if (userFieldList.size()) {
    std::vector<occa::memory> o_s;
    for(int i = 0; i < userFieldList.size(); i++) {
      o_s.push_back(o_AVG + i*fieldOffset);
    }   
    fld::write("tavg", atime, outfldCounter, o_s, outXYZ, FP64);
  }

  atime = 0; // reset
  outfldCounter++;
}

void tavg::outfld()
{
  tavg::outfld(/* outXYZ */ 0, /* FP64 */ 1);
}

deviceMemory<dfloat> tavg::o_avg()
{
  nekrsCheck(!setupCalled || !buildKernelCalled,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s\n",
             "called prior to tavg::setup()!");
  deviceMemory<dfloat> d_AVG(o_AVG);
  return d_AVG;
}

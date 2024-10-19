#include "platform.hpp"
#include "tavg.hpp"
#include "nekInterfaceAdapter.hpp"
#include "iofldFactory.hpp"

// private members
namespace
{
dlong fieldOffset;

std::vector< std::vector<deviceMemory<dfloat>> > userFieldList;
occa::memory o_AVG;

occa::kernel E1Kernel;
occa::kernel E2Kernel;
occa::kernel E3Kernel;
occa::kernel E4Kernel;

std::unique_ptr<iofld> fldWriter;

bool buildKernelCalled = false;
bool setupCalled = false;

int counter = 0;

double atime;
double timel;

int outfldCounter = 0;
} // namespace

static void E1(dlong N, dfloat a, dfloat b, int nflds, occa::memory o_x, occa::memory o_EX)
{
  E1Kernel(N, fieldOffset, nflds, a, b, o_x, o_EX);
}

static void E1(dlong N, dlong fieldOffset, dfloat a, dfloat b, int nflds, occa::memory o_x, occa::memory o_EX)
{
  E1Kernel(N, fieldOffset, nflds, a, b, o_x, o_EX);
}

static void
E2(dlong N, dfloat a, dfloat b, int nflds, occa::memory o_x, occa::memory o_y, occa::memory o_EXY)
{
  E2Kernel(N, fieldOffset, nflds, a, b, o_x, o_y, o_EXY);
}

static void E3(dlong N,
               dfloat a,
               dfloat b,
               int nflds,
               occa::memory o_x,
               occa::memory o_y,
               occa::memory o_z,
               occa::memory &o_EXYZ)
{
  E3Kernel(N, fieldOffset, nflds, a, b, o_x, o_y, o_z, o_EXYZ);
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

  E1Kernel = buildKernel("E1");
  E2Kernel = buildKernel("E2");
  E3Kernel = buildKernel("E3");
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
        E1(N, a, b, 1, entry.at(0), o_avg);
      } else if (entry.size() == 2) {
        E2(N, a, b, 1, entry.at(0), entry.at(1), o_avg);
      } else if (entry.size() == 3) {
        E3(N, a, b, 1, entry.at(0), entry.at(1), entry.at(2), o_avg);
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
  o_AVG = platform->device.malloc<double>(userFieldList.size() * fieldOffset);

  setupCalled = true;
}

void tavg::outfld(mesh_t *mesh)
{
  nekrsCheck(!setupCalled || !buildKernelCalled,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s\n",
             "called prior to tavg::setup()!");

  if (userFieldList.size() == 0) return;

  const bool outXYZ = mesh && outfldCounter == 0; 

  fldWriter = iofldFactory::create();

  if (!fldWriter->isInitialized()) {
    fldWriter->open(mesh, iofld::mode::write, "tavg");

    fldWriter->writeAttribute("precision", "64");
 
    fldWriter->addVariable("time", atime);
 
    for(int i = 0; i < userFieldList.size(); i++) {
      fldWriter->addVariable("scalar" + scalarDigitStr(i), std::vector<occa::memory>{o_AVG.slice(i * fieldOffset, mesh->Nlocal)});
    }  
  }

  fldWriter->writeAttribute("outputmesh", (outXYZ) ? "true" : "false");
  fldWriter->process(); 

  atime = 0; // reset
  outfldCounter++;
}

deviceMemory<double> tavg::o_avg()
{
  nekrsCheck(!setupCalled || !buildKernelCalled,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s\n",
             "called prior to tavg::setup()!");
  deviceMemory<double> d_AVG(o_AVG);
  return d_AVG;
}

void tavg::free()
{
  userFieldList.clear();
  o_AVG.free();
  fldWriter.reset();
}

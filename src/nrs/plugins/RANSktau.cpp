#include "nrs.hpp"
#include "platform.hpp"
#include "nekInterfaceAdapter.hpp"
#include "RANSktau.hpp"
#include "linAlg.hpp"

// private members
namespace
{
static nrs_t *nrs;

int kFieldIndex;

dfloat rho;
dfloat mueLam;

static occa::memory o_mut;

static occa::memory o_k;
static occa::memory o_tau;

static occa::memory o_implicitKtau;

static occa::kernel computeKernel;
static occa::kernel mueKernel;
static occa::kernel limitKernel;

static occa::kernel SijMag2OiOjSkKernel;

static bool buildKernelCalled = false;
static bool setupCalled = false;

static dfloat coeff[] = {
    0.6,       // sigma_k
    0.5,       // sigma_tau
    1.0,       // alpinf_str
    0.0708,    // beta0
    0.41,      // kappa
    0.09,      // betainf_str
    0.0,       // sigd_min
    1.0 / 8.0, // sigd_max
    400.0,     // fb_c1st
    400.0,     // fb_c2st
    85.0,      // fb_c1
    100.0,     // fb_c2
    0.52,      // alp_inf
    1e-8,      // TINY
    0          // Pope correction
};

occa::memory implicitK(double time, int scalarIdx)
{
  if (scalarIdx == kFieldIndex) {
    return o_implicitKtau.slice(0 * nrs->fieldOffset, nrs->fieldOffset);
  }
  if (scalarIdx == kFieldIndex + 1) {
    return o_implicitKtau.slice(1 * nrs->fieldOffset, nrs->fieldOffset);
  }
  return o_NULL;
}

} // namespace

void RANSktau::buildKernel(occa::properties _kernelInfo)
{
  occa::properties kernelInfo;
  if (!kernelInfo.get<std::string>("defines/p_sigma_k").size()) {
    kernelInfo["defines/p_sigma_k"] = coeff[0];
  }
  if (!kernelInfo.get<std::string>("defines/p_sigma_tau").size()) {
    kernelInfo["defines/p_sigma_tau"] = coeff[1];
  }
  if (!kernelInfo.get<std::string>("defines/p_alpinf_str").size()) {
    kernelInfo["defines/p_alpinf_str"] = coeff[2];
  }
  if (!kernelInfo.get<std::string>("defines/p_beta0").size()) {
    kernelInfo["defines/p_beta0"] = coeff[3];
  }
  if (!kernelInfo.get<std::string>("defines/p_kappa").size()) {
    kernelInfo["defines/p_kappa"] = coeff[4];
  }
  if (!kernelInfo.get<std::string>("defines/p_betainf_str").size()) {
    kernelInfo["defines/p_betainf_str"] = coeff[5];
  }
  if (!kernelInfo.get<std::string>("defines/p_ibetainf_str3").size()) {
    kernelInfo["defines/p_ibetainf_str3"] = 1 / pow(coeff[5], 3);
  }
  if (!kernelInfo.get<std::string>("defines/p_sigd_min").size()) {
    kernelInfo["defines/p_sigd_min"] = coeff[6];
  }
  if (!kernelInfo.get<std::string>("defines/p_sigd_max").size()) {
    kernelInfo["defines/p_sigd_max"] = coeff[7];
  }
  if (!kernelInfo.get<std::string>("defines/p_fb_c1st").size()) {
    kernelInfo["defines/p_fb_c1st"] = coeff[8];
  }
  if (!kernelInfo.get<std::string>("defines/p_fb_c2st").size()) {
    kernelInfo["defines/p_fb_c2st"] = coeff[9];
  }
  if (!kernelInfo.get<std::string>("defines/p_fb_c1").size()) {
    kernelInfo["defines/p_fb_c1"] = coeff[10];
  }
  if (!kernelInfo.get<std::string>("defines/p_fb_c2").size()) {
    kernelInfo["defines/p_fb_c2"] = coeff[11];
  }
  if (!kernelInfo.get<std::string>("defines/p_alp_inf").size()) {
    kernelInfo["defines/p_alp_inf"] = coeff[12];
  }
  if (!kernelInfo.get<std::string>("defines/p_tiny").size()) {
    kernelInfo["defines/p_tiny"] = coeff[13];
  }
  if (!kernelInfo.get<std::string>("defines/p_pope").size()) {
    kernelInfo["defines/p_pope"] = coeff[14];
  }

  const int verbose = platform->options.compareArgs("VERBOSE", "TRUE") ? 1 : 0;

  if (platform->comm.mpiRank == 0 && verbose) {
    std::cout << "\nRANSktau settings\n";
    std::cout << kernelInfo << std::endl;
  }

  kernelInfo += _kernelInfo;

  auto buildKernel = [&kernelInfo](const std::string &kernelName) {
    const auto path = getenv("NEKRS_KERNEL_DIR") + std::string("/nrs/plugins/");
    const auto fileName = path + "RANSktau.okl";
    const auto reqName = "RANSktau::";
    if (platform->options.compareArgs("REGISTER ONLY", "TRUE")) {
      platform->kernelRequests.add(reqName, fileName, kernelInfo);
      return occa::kernel();
    } else {
      buildKernelCalled = 1;
      return platform->kernelRequests.load(reqName, kernelName);
    }
  };

  computeKernel = buildKernel("RANSktauComputeHex3D");
  mueKernel = buildKernel("mue");
  limitKernel = buildKernel("limit");
  SijMag2OiOjSkKernel = buildKernel("SijMag2OiOjSk");

  int Nscalar;
  platform->options.getArgs("NUMBER OF SCALARS", Nscalar);

  nekrsCheck(Nscalar < 2, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "Nscalar needs to be >= 2!");
  platform->options.setArgs("VELOCITY STRESSFORMULATION", "TRUE");
}

void RANSktau::updateProperties()
{
  nekrsCheck(!setupCalled || !buildKernelCalled,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s\n",
             "called prior to tavg::setup()!");

  auto mesh = nrs->mesh;
  auto cds = nrs->cds;

  occa::memory o_mue = nrs->o_mue;
  occa::memory o_diff = cds->o_diff + cds->fieldOffsetScan[kFieldIndex];

  limitKernel(mesh->Nelements * mesh->Np, o_k, o_tau);
  mueKernel(mesh->Nelements * mesh->Np, nrs->fieldOffset, rho, mueLam, o_k, o_tau, o_mut, o_mue, o_diff);
}

const deviceMemory<dfloat> RANSktau::o_mue_t()
{
  deviceMemory<dfloat> out(o_mut);
  return out;
}

void RANSktau::updateSourceTerms()
{
  nekrsCheck(!setupCalled || !buildKernelCalled,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s\n",
             "called prior to tavg::setup()!");

  auto mesh = nrs->mesh;
  cds_t *cds = nrs->cds;

  occa::memory o_OiOjSk = platform->deviceMemoryPool.reserve<dfloat>(nrs->fieldOffset);
  occa::memory o_SijMag2 = platform->deviceMemoryPool.reserve<dfloat>(nrs->fieldOffset);

  occa::memory o_FS = cds->o_NLT + cds->fieldOffsetScan[kFieldIndex];

  auto o_SijOij = nrs->strainRotationRate();

  SijMag2OiOjSkKernel(mesh->Nelements * mesh->Np, nrs->fieldOffset, 1, o_SijOij, o_OiOjSk, o_SijMag2);

  computeKernel(mesh->Nelements,
                nrs->fieldOffset, // assumes offset is always the same
                rho,
                mueLam,
                mesh->o_vgeo,
                mesh->o_D,
                o_k,
                o_tau,
                o_SijMag2,
                o_OiOjSk,
                o_implicitKtau,
                o_FS);
}

void RANSktau::setup(int ifld)
{
  static bool isInitialized = false;
  if (isInitialized) {
    return;
  }
  isInitialized = true;

  nrs = dynamic_cast<nrs_t *>(platform->solver);
  kFieldIndex = ifld; // tauFieldIndex is assumed to be kFieldIndex+1

  platform->options.getArgs("VISCOSITY", mueLam);
  platform->options.getArgs("DENSITY", rho);

  for (int i = 0; i < 2; i++) {
    auto cds = nrs->cds;
    auto mesh = (kFieldIndex + i) ? cds->meshV : cds->mesh[0]; // only first scalar can be a CHT mesh
    auto o_rho = cds->o_prop.slice(cds->fieldOffsetSum + cds->fieldOffsetScan[kFieldIndex + i]);
    platform->linAlg->fill(mesh->Nlocal, rho, o_rho);

    const std::string sid = scalarDigitStr(kFieldIndex + i);
    nekrsCheck(!platform->options.getArgs("SCALAR" + sid + " DIFFUSIVITY").empty() ||
                   !platform->options.getArgs("SCALAR" + sid + " DENSITY").empty(),
               platform->comm.mpiComm,
               EXIT_FAILURE,
               "%s\n",
               "illegal property specificition for k/tau in par!");
  }

  auto cds = nrs->cds;

  nekrsCheck(cds->NSfields < kFieldIndex + 1,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "number of scalar fields too low!");

  o_k = cds->o_S + cds->fieldOffsetScan[kFieldIndex];
  o_tau = cds->o_S + cds->fieldOffsetScan[kFieldIndex + 1];

  o_mut = platform->device.malloc<dfloat>(cds->fieldOffset[kFieldIndex]);

  o_implicitKtau = platform->device.malloc<dfloat>(2 * nrs->fieldOffset);
  cds->userImplicitLinearTerm = implicitK;

  setupCalled = true;
}

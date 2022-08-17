#include "nrs.hpp"
#include "platform.hpp"
#include "nekInterfaceAdapter.hpp"
#include "RANSktau.hpp"
#include "linAlg.hpp"

// private members
namespace
{
static nrs_t* nrs;

int kFieldIndex;

dfloat rho;
dfloat mueLam;

static occa::memory o_mut;

static occa::memory o_k;
static occa::memory o_tau;

static occa::kernel computeKernel;
static occa::kernel SijOijKernel;
static occa::kernel SijOijMag2Kernel;
static occa::kernel limitKernel;
static occa::kernel mueKernel;

static bool setupCalled = 0;

static dfloat coeff[] = {
  0.6,                      // sigma_k
  0.5,                      // sigma_tau
  1.0,                      // alpinf_str
  0.0708,                   // beta0
  0.41,                     // kappa
  0.09,                     // betainf_str
  0.0,                      // sigd_min
  1.0 / 8.0,                // sigd_max
  400.0,                    // fb_c1st
  400.0,                    // fb_c2st
  85.0,                     // fb_c1
  100.0,                    // fb_c2
  0.52,                     // alp_inf
  1e-8,                     // TINY
  0                         // Pope correction
};
}

void RANSktau::buildKernel(occa::properties _kernelInfo)
{

  occa::properties kernelInfo;
  if(!kernelInfo.get<std::string>("defines/p_sigma_k").size())
    kernelInfo["defines/p_sigma_k"]       = coeff[0];
  if(!kernelInfo.get<std::string>("defines/p_sigma_tau").size())
    kernelInfo["defines/p_sigma_tau"]     = coeff[1];
  if(!kernelInfo.get<std::string>("defines/p_alpinf_str").size())
    kernelInfo["defines/p_alpinf_str"]    = coeff[2];
  if(!kernelInfo.get<std::string>("defines/p_beta0").size())
    kernelInfo["defines/p_beta0"]         = coeff[3];
  if(!kernelInfo.get<std::string>("defines/p_kappa").size())
    kernelInfo["defines/p_kappa"]         = coeff[4];
  if(!kernelInfo.get<std::string>("defines/p_betainf_str").size())
    kernelInfo["defines/p_betainf_str"]   = coeff[5];
  if(!kernelInfo.get<std::string>("defines/p_ibetainf_str3").size())
    kernelInfo["defines/p_ibetainf_str3"] = 1 / pow(coeff[5],3);
  if(!kernelInfo.get<std::string>("defines/p_sigd_min").size())
    kernelInfo["defines/p_sigd_min"]      = coeff[6];
  if(!kernelInfo.get<std::string>("defines/p_sigd_max").size())
    kernelInfo["defines/p_sigd_max"]      = coeff[7];
  if(!kernelInfo.get<std::string>("defines/p_fb_c1st").size())
    kernelInfo["defines/p_fb_c1st"]       = coeff[8];
  if(!kernelInfo.get<std::string>("defines/p_fb_c2st").size())
    kernelInfo["defines/p_fb_c2st"]       = coeff[9];
  if(!kernelInfo.get<std::string>("defines/p_fb_c1").size())
    kernelInfo["defines/p_fb_c1"]         = coeff[10];
  if(!kernelInfo.get<std::string>("defines/p_fb_c2").size())
    kernelInfo["defines/p_fb_c2"]         = coeff[11];
  if(!kernelInfo.get<std::string>("defines/p_alp_inf").size())
    kernelInfo["defines/p_alp_inf"]       = coeff[12];
  if(!kernelInfo.get<std::string>("defines/p_tiny").size())
    kernelInfo["defines/p_tiny"]          = coeff[13];
  if(!kernelInfo.get<std::string>("defines/p_pope").size())
    kernelInfo["defines/p_pope"]          = coeff[14];

  const int verbose = platform->options.compareArgs("VERBOSE","TRUE") ? 1:0;

  if(platform->comm.mpiRank == 0 && verbose) {
    std::cout << "\nRANSktau settings\n"; 
    std::cout << kernelInfo << std::endl;
  }

  kernelInfo += _kernelInfo;

  int rank = platform->comm.mpiRank;
  const std::string install_dir(getenv("NEKRS_INSTALL_DIR"));
  const std::string path = install_dir + "/okl/plugins/";
  std::string fileName, kernelName;
  const std::string extension = ".okl";
  {
      kernelName = "RANSktauComputeHex3D";
      fileName = path + kernelName + extension;
      computeKernel    = platform->device.buildKernel(fileName, kernelInfo, true);

      kernelName = "SijOijHex3D";
      fileName = install_dir + "/okl/nrs/" + kernelName + extension;
      SijOijKernel     = platform->device.buildKernel(fileName, kernelInfo, true);

      kernelName = "SijOijMag2";
      fileName = install_dir + "/okl/nrs/" + kernelName + extension;
      SijOijMag2Kernel = platform->device.buildKernel(fileName, kernelInfo, true);

      kernelName = "limit";
      fileName = path + kernelName + extension;
      limitKernel      = platform->device.buildKernel(fileName, kernelInfo, true);

      kernelName = "mue";
      fileName = path + kernelName + extension;
      mueKernel        = platform->device.buildKernel(fileName, kernelInfo, true);
  }

  int Nscalar;
  platform->options.getArgs("NUMBER OF SCALARS", Nscalar);

  if(Nscalar < 2) {
    if(platform->comm.mpiRank == 0) std::cout << "RANSktau: Nscalar needs to be >= 2!\n";
    ABORT(1);
  }
  platform->options.setArgs("STRESSFORMULATION", "TRUE");
}

void RANSktau::updateProperties()
{
  mesh_t* mesh = nrs->meshV;
  cds_t* cds = nrs->cds;

  occa::memory o_mue  = nrs->o_mue;
  occa::memory o_diff = cds->o_diff + cds->fieldOffsetScan[kFieldIndex] * sizeof(dfloat);

  limitKernel(mesh->Nelements * mesh->Np, o_k, o_tau);
  mueKernel(mesh->Nelements * mesh->Np,
            nrs->fieldOffset,
            rho,
            mueLam,
            o_k,
            o_tau,
            o_mut,
            o_mue,
            o_diff);
}

occa::memory RANSktau::o_mue_t()
{
  return o_mut;
}

void RANSktau::updateSourceTerms()
{
  mesh_t* mesh = nrs->meshV;
  cds_t* cds = nrs->cds;
  

  occa::memory o_OiOjSk  = platform->o_mempool.slice0;
  occa::memory o_SijMag2 = platform->o_mempool.slice1;
  occa::memory o_SijOij  = platform->o_mempool.slice2;

  occa::memory o_FS      = cds->o_FS     + cds->fieldOffsetScan[kFieldIndex] * sizeof(dfloat);
  occa::memory o_BFDiag  = cds->o_BFDiag + cds->fieldOffsetScan[kFieldIndex] * sizeof(dfloat);

  const int NSOfields = 9;
  SijOijKernel(mesh->Nelements,
               nrs->fieldOffset,
               1,
               mesh->o_vgeo,
               mesh->o_D,
               nrs->o_U,
               o_SijOij);

  ogsGatherScatterMany(o_SijOij,
                       NSOfields,
                       nrs->fieldOffset,
                       ogsDfloat,
                       ogsAdd,
                       mesh->ogs);

  platform->linAlg->axmyMany(
    mesh->Nlocal,
    NSOfields,
    nrs->fieldOffset,
    0,
    1.0,
    nrs->meshV->o_invLMM,
    o_SijOij);

  SijOijMag2Kernel(mesh->Nelements * mesh->Np,
                   nrs->fieldOffset,
                   1,
                   o_SijOij,
                   o_OiOjSk,
                   o_SijMag2);

  limitKernel(mesh->Nelements * mesh->Np, o_k, o_tau);

  computeKernel(mesh->Nelements,
                nrs->cds->fieldOffset[kFieldIndex],
                rho,
                mueLam,
                mesh->o_vgeo,
                mesh->o_D,
                o_k,
                o_tau,
                o_SijMag2,
                o_OiOjSk,
                o_BFDiag,
                o_FS);
}

void RANSktau::setup(nrs_t* nrsIn, dfloat mueIn, dfloat rhoIn, int ifld)
{
  setup(nrsIn, mueIn, rhoIn, ifld, NULL);
}

void RANSktau::setup(nrs_t* nrsIn, dfloat mueIn, dfloat rhoIn,
                     int ifld, const dfloat* coeffIn)
{
  if(setupCalled) return;
  
  nrs    = nrsIn;
  mueLam = mueIn;
  rho    = rhoIn;
  kFieldIndex = ifld;

  cds_t* cds = nrs->cds;
  mesh_t* mesh = nrs->meshV;

  if(coeffIn) memcpy(coeff, coeffIn, sizeof(coeff));

  o_k   = cds->o_S + cds->fieldOffsetScan[kFieldIndex] * sizeof(dfloat);
  o_tau = cds->o_S + cds->fieldOffsetScan[kFieldIndex+1] * sizeof(dfloat);

  o_mut = platform->device.malloc(cds->fieldOffset[kFieldIndex] ,  sizeof(dfloat));

  if(!cds->o_BFDiag.ptr()) {
    cds->o_BFDiag = platform->device.malloc(cds->fieldOffsetSum,  sizeof(dfloat));
    platform->linAlg->fill(cds->fieldOffsetSum, 0.0, cds->o_BFDiag);
  }

  setupCalled = 1;
}

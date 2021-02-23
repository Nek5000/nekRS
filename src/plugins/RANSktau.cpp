#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "RANSktau.hpp"

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
  1e-8                      // TINY
};
}

void RANSktau::buildKernel(nrs_t* nrs)
{
  mesh_t* mesh = nrs->mesh;

  occa::properties kernelInfo = *(nrs->kernelInfo);
  kernelInfo["defines/p_sigma_k"]       = coeff[0];
  kernelInfo["defines/p_sigma_tau"]     = coeff[1];
  kernelInfo["defines/p_alpinf_str"]    = coeff[2];
  kernelInfo["defines/p_beta0"]         = coeff[3];
  kernelInfo["defines/p_kappa"]         = coeff[4];
  kernelInfo["defines/p_betainf_str"]   = coeff[5];
  kernelInfo["defines/p_ibetainf_str3"] = 1 / pow(coeff[5],3);
  kernelInfo["defines/p_sigd_min"]      = coeff[6];
  kernelInfo["defines/p_sigd_max"]      = coeff[7];
  kernelInfo["defines/p_fb_c1st"]       = coeff[8];
  kernelInfo["defines/p_fb_c2st"]       = coeff[9];
  kernelInfo["defines/p_fb_c1"]         = coeff[10];
  kernelInfo["defines/p_fb_c2"]         = coeff[11];
  kernelInfo["defines/p_alp_inf"]       = coeff[12];
  kernelInfo["defines/p_tiny"]          = coeff[13];

  string fileName;
  int rank = mesh->rank;
  fileName.assign(getenv("NEKRS_INSTALL_DIR"));
  fileName += "/okl/plugins/RANSktau.okl";
  for (int r = 0; r < 2; r++) {
    if ((r == 0 && rank == 0) || (r == 1 && rank > 0)) {
      computeKernel    = mesh->device.buildKernel(fileName.c_str(), "computeHex3D", kernelInfo);
      SijOijKernel     = mesh->device.buildKernel(fileName.c_str(), "SijOijHex3D", kernelInfo);
      SijOijMag2Kernel = mesh->device.buildKernel(fileName.c_str(), "SijOijMag2", kernelInfo);
      limitKernel      = mesh->device.buildKernel(fileName.c_str(), "limit", kernelInfo);
      mueKernel        = mesh->device.buildKernel(fileName.c_str(), "mue", kernelInfo);
    }
    MPI_Barrier(mesh->comm);
  }

  if(nrs->Nscalar < 2) {
    if(mesh->rank == 0) cout << "RANSktau: Nscalar needs to be >= 2!\n";
    ABORT(1);
  }
  nrs->options.setArgs("STRESSFORMULATION", "TRUE");
}

void RANSktau::updateProperties()
{
  mesh_t* mesh = nrs->mesh;
  cds_t* cds = nrs->cds;

  occa::memory o_mue  = nrs->o_mue;
  occa::memory o_diff = cds->o_diff + kFieldIndex * cds->fieldOffset * sizeof(dfloat);

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
  mesh_t* mesh = nrs->mesh;
  cds_t* cds = nrs->cds;

  occa::memory o_OiOjSk  = nrs->o_wrk0;
  occa::memory o_SijMag2 = nrs->o_wrk1;
  occa::memory o_SijOij  = nrs->o_wrk2;

  occa::memory o_FS      = cds->o_FS     + kFieldIndex * cds->fieldOffset * sizeof(dfloat);
  occa::memory o_BFDiag  = cds->o_BFDiag + kFieldIndex * cds->fieldOffset * sizeof(dfloat);

  const int NSOfields = 9;
  SijOijKernel(mesh->Nelements,
               nrs->fieldOffset,
               mesh->o_vgeo,
               mesh->o_Dmatrices,
               nrs->o_U,
               o_SijOij);

  ogsGatherScatterMany(o_SijOij,
                       NSOfields,
                       nrs->fieldOffset,
                       ogsDfloat,
                       ogsAdd,
                       mesh->ogs);

  nrs->invMassMatrixKernel(
    mesh->Nelements,
    nrs->fieldOffset,
    NSOfields,
    mesh->o_vgeo,
    nrs->mesh->o_invLMM,
    o_SijOij);

  SijOijMag2Kernel(mesh->Nelements * mesh->Np,
                   nrs->fieldOffset,
                   o_SijOij,
                   o_OiOjSk,
                   o_SijMag2);

  limitKernel(mesh->Nelements * mesh->Np, o_k, o_tau);

  computeKernel(mesh->Nelements,
                nrs->cds->fieldOffset,
                rho,
                mueLam,
                mesh->o_vgeo,
                mesh->o_Dmatrices,
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
  mesh_t* mesh = nrs->mesh;

  if(coeffIn) memcpy(coeff, coeffIn, sizeof(coeff));

  o_k   = cds->o_S + kFieldIndex * cds->fieldOffset * sizeof(dfloat);
  o_tau = cds->o_S + (kFieldIndex + 1) * cds->fieldOffset * sizeof(dfloat);

  o_mut = mesh->device.malloc(cds->fieldOffset * sizeof(dfloat));

  if(!cds->o_BFDiag.ptr()) {
    cds->o_BFDiag = mesh->device.malloc(cds->NSfields * cds->fieldOffset * sizeof(dfloat));
    nrs->fillKernel(cds->NSfields * cds->fieldOffset, 0.0, cds->o_BFDiag);
  }

  setupCalled = 1;
}

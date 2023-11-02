#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "tavg.hpp"
#include "platform.hpp"
#include "linAlg.hpp"

// private members
namespace {
static ogs_t *ogs;
static nrs_t *nrs;

static occa::memory o_Uavg, o_Urms;
static occa::memory o_Urm2;
static occa::memory o_Pavg, o_Prms;
static occa::memory o_Savg, o_Srms;

tavg::fields userFieldList;
static occa::memory o_AVG;

static occa::kernel EXKernel;
static occa::kernel EXYKernel;
static occa::kernel EXYZKernel;
static occa::kernel E4Kernel;

static bool buildKernelCalled = false;
static bool setupCalled = false;

static int counter = 0;

static double atime;
static double timel;

static int outfldCounter = 0;

} // namespace


static void EX(dlong N, dfloat a, dfloat b, int nflds, occa::memory o_x, occa::memory o_EX)
{
  EXKernel(N, nrs->fieldOffset, nflds, a, b, o_x, o_EX);
}

static void EX(dlong N, dlong fieldOffset, dfloat a, dfloat b, int nflds, occa::memory o_x, occa::memory o_EX)
{
  EXKernel(N, fieldOffset, nflds, a, b, o_x, o_EX);
}

static void EXY(dlong N, dfloat a, dfloat b, int nflds, occa::memory o_x, occa::memory o_y, occa::memory o_EXY)
{
  EXYKernel(N, nrs->fieldOffset, nflds, a, b, o_x, o_y, o_EXY);
}

static void EXYZ(dlong N, dfloat a, dfloat b, int nflds, occa::memory o_x, occa::memory o_y, occa::memory o_z, 
                 occa::memory& o_EXYZ)
{
  EXYZKernel(N, nrs->fieldOffset, nflds, a, b, o_x, o_y, o_z, o_EXYZ);
}

static void E4(dlong N, dfloat a, dfloat b, int nflds, occa::memory o_1, occa::memory o_2, occa::memory o_3, 
               occa::memory o_4, occa::memory& o_E4)
{
  E4Kernel(N, nrs->fieldOffset, nflds, a, b, o_1, o_2, o_3, o_4, o_E4);
}


void tavg::buildKernel(occa::properties kernelInfo)
{
  static bool isInitialized = false;
  if (isInitialized) return;
  isInitialized = true;

  const std::string path = getenv("NEKRS_KERNEL_DIR") + std::string("/plugins/");
  std::string kernelName, fileName;
  const std::string extension = ".okl";
  {
    kernelName = "EX";
    fileName = path + kernelName + extension;
    EXKernel = platform->device.buildKernel(fileName, kernelInfo, true);

    kernelName = "EXY";
    fileName = path + kernelName + extension;
    EXYKernel = platform->device.buildKernel(fileName, kernelInfo, true);

    kernelName = "EXYZ";
    fileName = path + kernelName + extension;
    EXYZKernel = platform->device.buildKernel(fileName, kernelInfo, true);

    kernelName = "E4";
    fileName = path + kernelName + extension;
    E4Kernel = platform->device.buildKernel(fileName, kernelInfo, true);

  }
  buildKernelCalled = true;
}

void tavg::reset()
{
  counter = 0;
  atime = 0;
}

void tavg::run(double time)
{
  nrsCheck(!setupCalled || !buildKernelCalled, MPI_COMM_SELF, EXIT_FAILURE,
           "%s\n", "called prior to tavg::setup()!");

  if (!nrs->timeStepConverged)
    return;

  if (!counter) {
    atime = 0;
    timel = time;
  }
  counter++;

  const double dtime = time - timel;
  atime += dtime;

  if (atime == 0 || dtime == 0)
    return;

  const dfloat b = dtime / atime;
  const dfloat a = 1 - b;

  if(userFieldList.size()) {
    int cnt = 0;
    for(auto& entry : userFieldList) {
      auto o_avg = o_AVG.slice(cnt*nrs->fieldOffset, nrs->fieldOffset);
      const auto N = nrs->fieldOffset;

      if(entry.size() == 1) 
      {
        EX(N, a, b, 1, entry.at(0), o_avg);
      }
      else if(entry.size() == 2) 
      {
        EXY(N, a, b, 1, entry.at(0), entry.at(1), o_avg);
      } 
      else if(entry.size() == 3)
      {
        EXYZ(N, a, b, 1, entry.at(0), entry.at(1), entry.at(2), o_avg);
      }     
      else if(entry.size() == 4)
      {
        E4(N, a, b, 1, entry.at(0), entry.at(1), entry.at(2), entry.at(3), o_avg);
      }  
      cnt++;
    }
  } else {
    const auto N = nrs->fieldOffset;

    // velocity
    EX(N, a, b, nrs->NVfields, nrs->o_U, o_Uavg);
    EXY(N, a, b, nrs->NVfields, nrs->o_U, nrs->o_U, o_Urms);
 
    const dlong offset = nrs->fieldOffset;
    occa::memory o_vx = nrs->o_U + 0 * offset;
    occa::memory o_vy = nrs->o_U + 1 * offset;
    occa::memory o_vz = nrs->o_U + 2 * offset;
 
    EXY(N, a, b, 1, o_vx, o_vy, o_Urm2 + 0 * offset);
    EXY(N, a, b, 1, o_vy, o_vz, o_Urm2 + 1 * offset);
    EXY(N, a, b, 1, o_vz, o_vx, o_Urm2 + 2 * offset);
 
    // pressure
    EX(N, a, b, 1, nrs->o_P, o_Pavg);
    EXY(N, a, b, 1, nrs->o_P, nrs->o_P, o_Prms);
 
    // scalars
    if (nrs->Nscalar) {
      cds_t *cds = nrs->cds;
      EX(N, a, b, cds->NSfields, cds->o_S, o_Savg);
      EXY(N, a, b, cds->NSfields, cds->o_S, cds->o_S, o_Srms);
    }
  }

  timel = time;
}

void tavg::setup(nrs_t *nrs_, const fields& flds) 
{
  userFieldList = flds;

  for(auto& entry : userFieldList) {
    nrsCheck(entry.size() < 1 || entry.size() > 4, 
             platform->comm.mpiComm, EXIT_FAILURE,
             "%s\n", "invalid number of vectors in one of the user list entries!");

    for(auto& entry_i : entry) {
      nrsCheck(entry_i.length() < nrs_->fieldOffset, 
               platform->comm.mpiComm, EXIT_FAILURE,
               "%s\n", "vector size in one of the user list entries smaller than nrs_t::fieldOffset");
    }
  }
 
  setup(nrs_);
}

void tavg::setup(nrs_t *nrs_)
{
  static bool isInitialized = false;
  if (isInitialized) return;
  isInitialized = true;


  nrsCheck(!buildKernelCalled, MPI_COMM_SELF, EXIT_FAILURE,
           "%s\n", "called prior tavg::buildKernel()!");

  nrs = nrs_;

  if(userFieldList.size() == 0) {
    o_Uavg = platform->device.malloc<dfloat>(nrs->fieldOffset * nrs->NVfields);
    o_Urms = platform->device.malloc<dfloat>(nrs->fieldOffset * nrs->NVfields);
    platform->linAlg->fill(nrs->fieldOffset * nrs->NVfields, 0.0, o_Uavg);
    platform->linAlg->fill(nrs->fieldOffset * nrs->NVfields, 0.0, o_Urms);
 
    o_Urm2 = platform->device.malloc<dfloat>(nrs->fieldOffset * nrs->NVfields);
    platform->linAlg->fill(nrs->fieldOffset * nrs->NVfields, 0.0, o_Urm2);
 
    o_Pavg = platform->device.malloc<dfloat>(nrs->fieldOffset);
    o_Prms = platform->device.malloc<dfloat>(nrs->fieldOffset);
    platform->linAlg->fill(nrs->fieldOffset, 0.0, o_Pavg);
    platform->linAlg->fill(nrs->fieldOffset, 0.0, o_Prms);
 
    if (nrs->Nscalar) {
      cds_t *cds = nrs->cds;
      o_Savg = platform->device.malloc<dfloat>(cds->fieldOffsetSum);
      o_Srms = platform->device.malloc<dfloat>(cds->fieldOffsetSum);
      platform->linAlg->fill(cds->fieldOffsetSum, 0.0, o_Savg);
      platform->linAlg->fill(nrs->fieldOffset, 0.0, o_Srms);
    }
  } else {
    o_AVG = platform->device.malloc<dfloat>(userFieldList.size()*nrs->fieldOffset); 
  }

  setupCalled = true;
}

void tavg::outfld(int _outXYZ, int FP64)
{
  nrsCheck(!setupCalled || !buildKernelCalled, MPI_COMM_SELF, EXIT_FAILURE,
           "%s\n", "called prior to tavg::setup()!");

  if (!nrs->timeStepConverged)
    return;

  int outXYZ = _outXYZ;
  if (!outfldCounter)
    outXYZ = 1;

  occa::memory o_Tavg, o_Trms;

  if(userFieldList.size()) {
    writeFld("ust", atime, outfldCounter, outXYZ, FP64, o_NULL, o_NULL, o_AVG, userFieldList.size());
  } else {
    const int Nscalar = nrs->Nscalar;
    if (nrs->Nscalar) {
      o_Tavg = o_Savg;
      o_Trms = o_Srms;
    }
 
    writeFld("avg", atime, outfldCounter, outXYZ, FP64, o_Uavg, o_Pavg, o_Tavg, Nscalar);
 
    writeFld("rms", atime, outfldCounter, outXYZ, FP64, o_Urms, o_Prms, o_Trms, Nscalar);
 
    writeFld("rm2", atime, outfldCounter, outXYZ, FP64, o_Urm2, o_NULL, o_NULL, 0);
  }

  atime = 0;
  outfldCounter++;
}

void tavg::outfld()
{ 
  tavg::outfld(/* outXYZ */ 0, /* FP64 */ 1); 
}

occa::memory tavg::o_avg()
{
  nrsCheck(!setupCalled || !buildKernelCalled, MPI_COMM_SELF, EXIT_FAILURE,
           "%s\n", "called prior to tavg::setup()!");

  return o_AVG;
}

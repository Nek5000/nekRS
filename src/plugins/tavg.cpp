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
static occa::memory o_userFieldAvg;

static occa::kernel EXKernel;
static occa::kernel EXXKernel;
static occa::kernel EXYKernel;
static occa::kernel EXYZKernel;

static bool buildKernelCalled = 0;
static bool setupCalled = 0;

static int counter = 0;

static dfloat atime;
static dfloat timel;

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

static void EXX(dlong N, dfloat a, dfloat b, int nflds, occa::memory o_x, occa::memory o_EXX)
{
  EXXKernel(N, nrs->fieldOffset, nflds, a, b, o_x, o_EXX);
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


void tavg::buildKernel(occa::properties kernelInfo)
{

  std::string path;
  int rank = platform->comm.mpiRank;
  path.assign(getenv("NEKRS_INSTALL_DIR"));
  path += "/kernels/plugins/";
  std::string kernelName, fileName;
  const std::string extension = ".okl";
  {
    kernelName = "EX";
    fileName = path + kernelName + extension;
    EXKernel = platform->device.buildKernel(fileName, kernelInfo, true);

    kernelName = "EXX";
    fileName = path + kernelName + extension;
    EXXKernel = platform->device.buildKernel(fileName, kernelInfo, true);

    kernelName = "EXY";
    fileName = path + kernelName + extension;
    EXYKernel = platform->device.buildKernel(fileName, kernelInfo, true);

    kernelName = "EXYZ";
    fileName = path + kernelName + extension;
    EXYZKernel = platform->device.buildKernel(fileName, kernelInfo, true);
  }
  buildKernelCalled = 1;
}

void tavg::reset()
{
  counter = 0;
  atime = 0;
}

void tavg::run(dfloat time)
{
  if (!setupCalled || !buildKernelCalled) {
    if(platform->comm.mpiRank == 0)
      std::cout << "tavg::run() was called prior to tavg::setup()!\n";
    ABORT(1);
  }

  if (!nrs->timeStepConverged)
    return;

  if (!counter) {
    atime = 0;
    timel = time;
  }
  counter++;

  const dfloat dtime = time - timel;
  atime += dtime;

  if (atime == 0 || dtime == 0)
    return;

  const dfloat b = dtime / atime;
  const dfloat a = 1 - b;

  mesh_t *mesh = nrs->meshV;
  const dlong N = mesh->Nelements * mesh->Np;

  if(userFieldList.size()) {
    int cnt = 0;
    for(auto& entry : userFieldList) {
      auto o_avg = o_userFieldAvg.slice(cnt*nrs->fieldOffset*sizeof(dfloat), nrs->fieldOffset*sizeof(dfloat));
      if(entry.size() == 1) 
      {
        EX(N, a, b, 1, entry.at(0), o_avg);
      } else 
      if(entry.size() == 2) 
      {
        if(entry.at(0).ptr() == entry.at(1).ptr())
          EXX(N, a, b, 1, entry.at(0), o_avg);
        else
          EXY(N, a, b, 1, entry.at(0), entry.at(1), o_avg);
      } else
      if(entry.size() == 3)
      {
        EXYZ(N, a, b, 1, entry.at(0), entry.at(1), entry.at(2), o_avg);
      }      
      cnt++;
    }
  } else {
    // velocity
    EX(N, a, b, nrs->NVfields, nrs->o_U, o_Uavg);
    EXX(N, a, b, nrs->NVfields, nrs->o_U, o_Urms);
 
    const dlong offsetByte = nrs->fieldOffset * sizeof(dfloat);
    occa::memory o_vx = nrs->o_U + 0 * offsetByte;
    occa::memory o_vy = nrs->o_U + 1 * offsetByte;
    occa::memory o_vz = nrs->o_U + 2 * offsetByte;
 
    EXY(N, a, b, 1, o_vx, o_vy, o_Urm2 + 0 * offsetByte);
    EXY(N, a, b, 1, o_vy, o_vz, o_Urm2 + 1 * offsetByte);
    EXY(N, a, b, 1, o_vz, o_vx, o_Urm2 + 2 * offsetByte);
 
    // pressure
    EX(N, a, b, 1, nrs->o_P, o_Pavg);
    EXX(N, a, b, 1, nrs->o_P, o_Prms);
 
    // scalars
    if (nrs->Nscalar) {
      cds_t *cds = nrs->cds;
      const dlong N = cds->mesh[0]->Nelements * cds->mesh[0]->Np;
      EX(N, a, b, cds->NSfields, cds->o_S, o_Savg);
      EXX(N, a, b, cds->NSfields, cds->o_S, o_Srms);
    }
  }

  timel = time;
}

void tavg::setup(nrs_t *nrs_, const fields& flds) 
{
  userFieldList = flds;

  for(auto& entry : userFieldList) {
    if(entry.size() > 3) {
      if(platform->comm.mpiRank == 0)
        std::cout << "tavg supports fields up to 3 entries!\n";
      ABORT(1);
    }
 }
 
  setup(nrs_);
}

void tavg::setup(nrs_t *nrs_)
{
  if (setupCalled) {
    if(platform->comm.mpiRank == 0)
      std::cout << "invalid second call to tavg::setup()!\n";
    ABORT(1);
  }

  if (!buildKernelCalled) {
    if(platform->comm.mpiRank == 0)
      std::cout << "tavg::setup() was called prior tavg::buildKernel()!\n";
    ABORT(1);
  }

  nrs = nrs_;
  mesh_t *mesh = nrs->meshV;

  if(userFieldList.size() == 0) {
    o_Uavg = platform->device.malloc(nrs->fieldOffset * nrs->NVfields, sizeof(dfloat));
    o_Urms = platform->device.malloc(nrs->fieldOffset * nrs->NVfields, sizeof(dfloat));
    platform->linAlg->fill(nrs->fieldOffset * nrs->NVfields, 0.0, o_Uavg);
    platform->linAlg->fill(nrs->fieldOffset * nrs->NVfields, 0.0, o_Urms);
 
    o_Urm2 = platform->device.malloc(nrs->fieldOffset * nrs->NVfields, sizeof(dfloat));
    platform->linAlg->fill(nrs->fieldOffset * nrs->NVfields, 0.0, o_Urm2);
 
    o_Pavg = platform->device.malloc(nrs->fieldOffset, sizeof(dfloat));
    o_Prms = platform->device.malloc(nrs->fieldOffset, sizeof(dfloat));
    platform->linAlg->fill(nrs->fieldOffset, 0.0, o_Pavg);
    platform->linAlg->fill(nrs->fieldOffset, 0.0, o_Prms);
 
    if (nrs->Nscalar) {
      cds_t *cds = nrs->cds;
      o_Savg = platform->device.malloc(cds->fieldOffsetSum, sizeof(dfloat));
      o_Srms = platform->device.malloc(cds->fieldOffsetSum, sizeof(dfloat));
      platform->linAlg->fill(cds->fieldOffsetSum, 0.0, o_Savg);
      platform->linAlg->fill(nrs->fieldOffset, 0.0, o_Srms);
    }
  } else {
    o_userFieldAvg = platform->device.malloc(userFieldList.size()*nrs->fieldOffset*sizeof(dfloat)); 
  }

  setupCalled = 1;
}

void tavg::outfld(int _outXYZ, int FP64)
{
  if (!setupCalled || !buildKernelCalled) {
    if(platform->comm.mpiRank == 0)
      std::cout << "tavg::outfld() was called prior to tavg::setup()!\n";
    ABORT(1);
  }

  if (!nrs->timeStepConverged)
    return;

  cds_t *cds = nrs->cds;
  mesh_t *mesh = nrs->meshV;

  int outXYZ = _outXYZ;
  if (!outfldCounter)
    outXYZ = 1;

  occa::memory o_Tavg, o_Trms;

  if(userFieldList.size()) {
    writeFld("ust", atime, outfldCounter, outXYZ, FP64, &o_NULL, &o_NULL, &o_userFieldAvg, userFieldList.size());
  } else {
    const int Nscalar = nrs->Nscalar;
    if (nrs->Nscalar) {
      o_Tavg = o_Savg;
      o_Trms = o_Srms;
    }
 
    writeFld("avg", atime, outfldCounter, outXYZ, FP64, &o_Uavg, &o_Pavg, &o_Tavg, Nscalar);
 
    writeFld("rms", atime, outfldCounter, outXYZ, FP64, &o_Urms, &o_Prms, &o_Trms, Nscalar);
 
    writeFld("rm2", atime, outfldCounter, outXYZ, FP64, &o_Urm2, &o_NULL, &o_NULL, 0);
  }

  atime = 0;
  outfldCounter++;
}

void tavg::outfld()
{ 
  tavg::outfld(/* outXYZ */ 0, /* FP64 */ 1); 
}

occa::memory tavg::userFieldAvg()
{
  if (!setupCalled || !buildKernelCalled) {
    if(platform->comm.mpiRank == 0)
      std::cout << "tavg::userFieldAvg() was called prior to tavg::setup()!\n";
    ABORT(1);
  }

  return o_userFieldAvg;
}

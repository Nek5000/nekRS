#include "nrs.hpp"
#include "platform.hpp"
#include "nekInterfaceAdapter.hpp"
#include "velRecycling.hpp"

// private members
namespace
{
oogs_t *ogs;
nrs_t *nrs;

occa::memory o_wrk;

occa::kernel setBCVectorValueKernel;
occa::kernel maskCopyKernel;

pointInterpolation_t *interp;
occa::memory o_Uint;
occa::memory o_maskIds;

int bID;
occa::memory o_bID;
dfloat wbar;
dfloat area;

static bool buildKernelCalled = false;
static bool setupCalled = false;

} // namespace

static void _setup(occa::memory &o_wrk_, const int bID_, const dfloat wbar_)
{
  static bool isInitialized = false;
  if (isInitialized) {
    return;
  }
  isInitialized = true;

  nrs = dynamic_cast<nrs_t*>(platform->solver);
  o_wrk = o_wrk_;
  bID = bID_;
  wbar = wbar_;

  auto mesh = nrs->mesh;

  nekrsCheck(o_wrk.length() < nrs->NVfields * nrs->fieldOffset,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "o_wrk too small!\n");

  {
    std::vector<int> tmp{bID};
    o_bID = platform->device.malloc<int>(tmp.size());
    o_bID.copyFrom(tmp.data());
  }

  auto o_tmp = platform->device.malloc<dfloat>(mesh->Nlocal);
  platform->linAlg->fill(mesh->Nlocal, 1.0, o_tmp);
  area = mesh->surfaceAreaMultiplyIntegrate(o_bID.length(), o_bID, o_tmp).at(0);

  setupCalled = true;
}

void velRecycling::buildKernel(occa::properties kernelInfo)
{
  auto buildKernel = [&kernelInfo](const std::string& kernelName)
  {
    const auto path = getenv("NEKRS_KERNEL_DIR") + std::string("/nrs/plugins/");
    const auto fileName = path + "velRecycling.okl";
    const auto reqName = "velRecycling::";
    if (platform->options.compareArgs("REGISTER ONLY", "TRUE")) {
      platform->kernelRequests.add(reqName, fileName, kernelInfo);
      return occa::kernel();
    } else {
      buildKernelCalled = 1;
      return platform->kernelRequests.load(reqName, kernelName);
    }
  };

  setBCVectorValueKernel = buildKernel("setBCVectorValue");
  maskCopyKernel = buildKernel("velRecyclingMaskCopy");
}

void velRecycling::copy()
{
  nekrsCheck(!setupCalled || !buildKernelCalled,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s\n",
             "called prior to tavg::setup()!");

  auto mesh = nrs->mesh;

  if (interp) {
    const dlong offset = o_Uint.length() / nrs->NVfields;
    interp->eval(nrs->NVfields, nrs->fieldOffset, nrs->o_U, offset, o_Uint);
    maskCopyKernel(interp->numPoints(), offset, nrs->fieldOffset, o_maskIds, o_Uint, o_wrk);
  } else {
    o_wrk.copyFrom(nrs->o_U, nrs->NVfields * nrs->fieldOffset);
    const dfloat zero = 0.0;
    setBCVectorValueKernel(mesh->Nelements, zero, bID, nrs->fieldOffset, o_wrk, mesh->o_vmapM, mesh->o_EToB);
    oogs::startFinish(o_wrk, nrs->NVfields, nrs->fieldOffset, ogsDfloat, ogsAdd, ogs);
  }

  auto flux = mesh->surfaceAreaNormalMultiplyIntegrate(nrs->fieldOffset, o_bID.length(), o_bID, o_wrk);

  const dfloat scale = -wbar * area / flux[0];
  nekrsCheck(std::isnan(scale) || std::isinf(scale), 
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s\n",
             "unreasonable scaling factor!");

  platform->linAlg->scale(nrs->NVfields * nrs->fieldOffset, scale, o_wrk);
}

void velRecycling::setup(occa::memory o_wrk_,
                         const hlong eOffset,
                         const int bID_,
                         const dfloat wbar_)
{
  _setup(o_wrk_, bID_, wbar_);

  auto mesh = nrs->mesh;

  const dlong Ntotal = mesh->Np * mesh->Nelements;

  // establish a unique numbering
  // relies on a special global element numbering (extruded mesh)
  auto ids = (hlong *)calloc(Ntotal, sizeof(hlong));
  for (int e = 0; e < mesh->Nelements; e++) {
    auto eg = nek::localElementIdToGlobal(e);

    for (int n = 0; n < mesh->Np; n++) {
      ids[e * mesh->Np + n] = eg * mesh->Np + (n + 1);
    }

    for (int n = 0; n < mesh->Nfp * mesh->Nfaces; n++) {
      const int f = n / mesh->Nfp;
      const int idM = mesh->vmapM[e * mesh->Nfp * mesh->Nfaces + n];
      if (mesh->EToB[f + e * mesh->Nfaces] == bID) {
        ids[idM] += eOffset * mesh->Np;
      }
    }
  }

  ogs = oogs::setup(Ntotal,
                    ids,
                    nrs->NVfields,
                    nrs->fieldOffset,
                    ogsDfloat,
                    platform->comm.mpiComm,
                    0,
                    platform->device.occaDevice(),
                    NULL,
                    OOGS_AUTO);

  free(ids);
}

void velRecycling::setup(occa::memory o_wrk_,
                         const dfloat xOffset,
                         const dfloat yOffset,
                         const dfloat zOffset,
                         const int bID_,
                         const dfloat wbar_)
{
  _setup(o_wrk_, bID_, wbar_);

  auto mesh = nrs->mesh;

  int cnt = 0;
  for (int e = 0; e < mesh->Nelements; e++) {
    for (int n = 0; n < mesh->Nfp * mesh->Nfaces; n++) {
      const int f = n / mesh->Nfp;
      if (mesh->EToB[f + e * mesh->Nfaces] == bID) {
        cnt++;
      }
    }
  }

  const auto nPoints = cnt;

  o_Uint = platform->device.malloc<dfloat>(nrs->NVfields * alignStride<dfloat>(nPoints));
  o_maskIds = platform->device.malloc<dlong>(nPoints);
  std::vector<dlong> maskIds(nPoints);

  std::vector<dfloat> xBid(nPoints);
  std::vector<dfloat> yBid(nPoints);
  std::vector<dfloat> zBid(nPoints);

  auto [x, y, z] = mesh->xyzHost();

  cnt = 0;
  for (int e = 0; e < mesh->Nelements; e++) {
    for (int n = 0; n < mesh->Nfp * mesh->Nfaces; n++) {
      const int f = n / mesh->Nfp;
      const int idM = mesh->vmapM[e * mesh->Nfp * mesh->Nfaces + n];
      if (mesh->EToB[f + e * mesh->Nfaces] == bID) {
        maskIds[cnt] = idM;
        xBid[cnt] = x[idM] + xOffset;
        yBid[cnt] = y[idM] + yOffset;
        zBid[cnt] = z[idM] + zOffset;
        cnt++;
      }
    }
  }
  o_maskIds.copyFrom(maskIds.data());

  interp = new pointInterpolation_t(mesh, platform->comm.mpiComm);

  auto o_xBid = platform->device.malloc<dfloat>(nPoints, xBid.data());
  auto o_yBid = platform->device.malloc<dfloat>(nPoints, yBid.data());
  auto o_zBid = platform->device.malloc<dfloat>(nPoints, zBid.data());

  interp->setPoints(o_xBid, o_yBid, o_zBid);
  const auto verbosity = pointInterpolation_t::VerbosityLevel::Basic;
  interp->find(verbosity);
}

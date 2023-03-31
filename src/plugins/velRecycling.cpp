/*
   copy velocity data of a given slab (slabIdSrc) to another slab
   (slideIdDst) also known as recycling

   Note: This implementation relies on a special global element
         numbering which is only true for extruded meshes in z from nek!
 */

#include "nrs.hpp"
#include "platform.hpp"
#include "nekInterfaceAdapter.hpp"
#include "velRecycling.hpp"

// private members
namespace {
static oogs_t *ogs;
static nrs_t *nrs;

static occa::memory o_wrk;

static dfloat *flux, *area;
static occa::memory o_flux, o_area;

static dfloat *tmp1, *tmp2;
static occa::memory o_tmp1, o_tmp2;

static occa::kernel setBCVectorValueKernel;
static occa::kernel getBCFluxKernel;
static occa::kernel sumReductionKernel;

static bool buildKernelCalled = 0;
static bool setupCalled = 0;

static int bID;
static dfloat wbar;

static int Nblock;
} // namespace

void velRecycling::buildKernel(occa::properties kernelInfo)
{
  const std::string path = getenv("NEKRS_KERNEL_DIR") + std::string("/plugins/");

  std::string fileName, kernelName;
  const std::string extension = ".okl";
  {
    kernelName = "setBCVectorValue";
    fileName = path + kernelName + extension;
    setBCVectorValueKernel = platform->device.buildKernel(fileName, kernelInfo, true);

    kernelName = "getBCFlux";
    fileName = path + kernelName + extension;
    getBCFluxKernel = platform->device.buildKernel(fileName, kernelInfo, true);

    kernelName = "sumReduction";
    fileName = path + kernelName + extension;
    sumReductionKernel = platform->device.buildKernel(fileName, kernelInfo, true);
  }
}

void velRecycling::copy()
{
  mesh_t *mesh = nrs->meshV;

  const dfloat zero = 0.0;

  // copy recycling plane in interior to inlet
  o_wrk.copyFrom(nrs->o_U, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));
  setBCVectorValueKernel(mesh->Nelements, zero, bID, nrs->fieldOffset, o_wrk, mesh->o_vmapM, mesh->o_EToB);

  oogs::startFinish(o_wrk, nrs->NVfields, nrs->fieldOffset, ogsDfloat, ogsAdd, ogs);

  // rescale
  getBCFluxKernel(mesh->Nelements,
                  bID,
                  nrs->fieldOffset,
                  o_wrk,
                  mesh->o_vmapM,
                  mesh->o_EToB,
                  mesh->o_sgeo,
                  o_area,
                  o_flux);

  const int NfpTotal = mesh->Nelements * mesh->Nfaces * mesh->Nfp;
  sumReductionKernel(NfpTotal, o_area, o_flux, o_tmp1, o_tmp2);

  o_tmp1.copyTo(tmp1);
  o_tmp2.copyTo(tmp2);
  dfloat sbuf[2] = {0, 0};
  for (int n = 0; n < Nblock; n++) {
    sbuf[0] += tmp1[n];
    sbuf[1] += tmp2[n];
  }
  MPI_Allreduce(MPI_IN_PLACE, sbuf, 2, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);

  const dfloat scale = -wbar * sbuf[0] / sbuf[1];
  // printf("rescaling inflow: %f\n", scale);
  platform->linAlg->scale(nrs->NVfields * nrs->fieldOffset, scale, o_wrk);
}

void velRecycling::setup(nrs_t *nrs_,
                         occa::memory o_wrk_,
                         const hlong eOffset,
                         const int bID_,
                         const dfloat wbar_)
{
  nrs = nrs_;
  o_wrk = o_wrk_;
  bID = bID_;
  wbar = wbar_;

  mesh_t *mesh = nrs->meshV;

  const dlong Ntotal = mesh->Np * mesh->Nelements;
  hlong *ids = (hlong *)calloc(Ntotal, sizeof(hlong));

  for (int e = 0; e < mesh->Nelements; e++) {
    // establish a unique numbering
    const hlong eg = nek::lglel(e); // 0-based

    for (int n = 0; n < mesh->Np; n++)
      ids[e * mesh->Np + n] = eg * mesh->Np + (n + 1);

    for (int n = 0; n < mesh->Nfp * mesh->Nfaces; n++) {
      const int f = n / mesh->Nfp;
      const int idM = mesh->vmapM[e * mesh->Nfp * mesh->Nfaces + n];
      if (mesh->EToB[f + e * mesh->Nfaces] == bID)
        ids[idM] += eOffset * mesh->Np;
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

  const int NfpTotal = mesh->Nelements * mesh->Nfaces * mesh->Nfp;

  Nblock = (NfpTotal + BLOCKSIZE - 1) / BLOCKSIZE;
  tmp1 = (dfloat *)calloc(Nblock, sizeof(dfloat));
  tmp2 = (dfloat *)calloc(Nblock, sizeof(dfloat));

  o_tmp1 = platform->device.malloc(Nblock * sizeof(dfloat), tmp1);
  o_tmp2 = platform->device.malloc(Nblock * sizeof(dfloat), tmp2);

  flux = (dfloat *)calloc(NfpTotal, sizeof(dfloat));
  area = (dfloat *)calloc(NfpTotal, sizeof(dfloat));
  o_flux = platform->device.malloc(NfpTotal * sizeof(dfloat), flux);
  o_area = platform->device.malloc(NfpTotal * sizeof(dfloat), area);
}

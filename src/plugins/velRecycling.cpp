/*
   copy velocity data of a given slab (slabIdSrc) to another slab
   (slideIdDst) also known as recycling

   Note: This implementation relies on a special global element
         numbering which is only true for extruded meshes in z from nek!
 */

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "velRecycling.hpp"

// private members
namespace
{
static ogs_t* ogs;
static nrs_t* nrs;

static occa::memory o_wrk;

static dfloat* flux, * area;
static occa::memory o_flux, o_area;

static dfloat* tmp1, * tmp2;
static occa::memory o_tmp1, o_tmp2;

static occa::kernel setBCVectorValueKernel;
static occa::kernel getBCFluxKernel;
static occa::kernel sumReductionKernel;
static occa::kernel scalarMultiplyKernel;

static bool buildKernelCalled = 0;
static bool setupCalled = 0;

static int bID;
static dfloat wbar;

static int Nblock;
}

void velRecycling::buildKernel(nrs_t* nrs)
{
  mesh_t* mesh = nrs->mesh;

  string fileName;
  int rank = mesh->rank;
  fileName.assign(getenv("NEKRS_INSTALL_DIR"));
  fileName += "/okl/plugins/velRecycling.okl";
  occa::properties& kernelInfo = *nrs->kernelInfo;
  for (int r = 0; r < 2; r++) {
    if ((r == 0 && rank == 0) || (r == 1 && rank > 0)) {
      setBCVectorValueKernel =  mesh->device.buildKernel(fileName.c_str(),
                                                         "setBCVectorValue",
                                                         kernelInfo);
      getBCFluxKernel        =  mesh->device.buildKernel(fileName.c_str(), "getBCFlux", kernelInfo);
      sumReductionKernel     =  mesh->device.buildKernel(fileName.c_str(),
                                                         "sumReduction",
                                                         kernelInfo);
      scalarMultiplyKernel   =  mesh->device.buildKernel(fileName.c_str(),
                                                         "scalarMultiply",
                                                         kernelInfo);
    }
    MPI_Barrier(mesh->comm);
  }
}

void velRecycling::copy()
{
  mesh_t* mesh = nrs->mesh;
  const dfloat zero = 0.0;

  // copy recycling plane in interior to inlet
  o_wrk.copyFrom(nrs->o_U, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));
  setBCVectorValueKernel(mesh->Nelements, zero, bID, nrs->fieldOffset,
                         o_wrk, mesh->o_vmapM, mesh->o_EToB);

  ogsGatherScatterMany(o_wrk, nrs->NVfields, nrs->fieldOffset,
                       ogsDfloat, ogsAdd, ogs);

/*
   for(int k=0;k<nrs->dim;++k)
    ogsGatherScatter(o_wrk+k*nrs->fieldOffset*sizeof(dfloat),
                     ogsDfloat, ogsAdd, ogs);
 */

  // rescale
  getBCFluxKernel(mesh->Nelements, bID, nrs->fieldOffset, o_wrk,
                  mesh->o_vmapM, mesh->o_EToB, mesh->o_sgeo, o_area, o_flux);

  const int NfpTotal = mesh->Nelements * mesh->Nfaces * mesh->Nfp;
  sumReductionKernel(NfpTotal, o_area, o_flux, o_tmp1, o_tmp2);

  o_tmp1.copyTo(tmp1);
  o_tmp2.copyTo(tmp2);
  dfloat sbuf[2] = {0,0};
  for(int n = 0; n < Nblock; n++) {
    sbuf[0] += tmp1[n];
    sbuf[1] += tmp2[n];
  }
  MPI_Allreduce(MPI_IN_PLACE, sbuf, 2, MPI_DFLOAT, MPI_SUM, mesh->comm);

  const dfloat scale = -wbar * sbuf[0] / sbuf[1];
  //printf("rescaling inflow: %f\n", scale);
  scalarMultiplyKernel(nrs->NVfields * nrs->fieldOffset, scale, o_wrk);
}

void velRecycling::setup(nrs_t* nrs_, occa::memory o_wrk_, const hlong eOffset, const int bID_,
                         const dfloat wbar_)
{
  nrs = nrs_;
  o_wrk = o_wrk_;
  bID = bID_;
  wbar = wbar_;

  mesh_t* mesh = nrs->mesh;

  const dlong Ntotal = mesh->Np * mesh->Nelements;
  hlong* ids = (hlong*) calloc(Ntotal, sizeof(hlong));

  for (int e = 0; e < mesh->Nelements; e++) {
    // establish a unique numbering
    const hlong eg = nek_lglel(e); // 0-based

    for (int n = 0; n < mesh->Np; n++)
      ids[e * mesh->Np + n] = eg * mesh->Np + (n + 1);

    for (int n = 0; n < mesh->Nfp * mesh->Nfaces; n++) {
      const int f = n / mesh->Nfp;
      const int idM = mesh->vmapM[e * mesh->Nfp * mesh->Nfaces + n];
      if (mesh->EToB[f + e * mesh->Nfaces] == bID)
        ids[idM] += eOffset * mesh->Np;
    }
  }

  ogs = ogsSetup(Ntotal, ids, mesh->comm, 1, mesh->device);
  free(ids);

  const int NfpTotal = mesh->Nelements * mesh->Nfaces * mesh->Nfp;

  Nblock = (NfpTotal + BLOCKSIZE - 1) / BLOCKSIZE;
  tmp1   = (dfloat*) calloc(Nblock, sizeof(dfloat));
  tmp2   = (dfloat*) calloc(Nblock, sizeof(dfloat));

  o_tmp1 = mesh->device.malloc(Nblock * sizeof(dfloat), tmp1);
  o_tmp2 = mesh->device.malloc(Nblock * sizeof(dfloat), tmp2);

  flux   = (dfloat*)calloc(NfpTotal, sizeof(dfloat));
  area   = (dfloat*)calloc(NfpTotal, sizeof(dfloat));
  o_flux = mesh->device.malloc(NfpTotal * sizeof(dfloat), flux);
  o_area = mesh->device.malloc(NfpTotal * sizeof(dfloat), area);
}

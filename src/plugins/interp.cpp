
#include <mpi.h>
#include "nrs.hpp"
#include "platform.hpp"
#include <vector>

#include "interp.hpp"

interp_t::interp_t(nrs_t* nrs, double newton_tol)
  : nrs(nrs_), newton_tol_(newton_tol)
{

  if (newton_tol < 5e-13) {
    newton_tol = 5e-13;
  }
  int npt_max = 128;
  int bb_tol = 0.01;

  mesh_t* mesh = nrs->meshV;

  dlong nmsh = mesh->N;
  dlong nelm = mesh->Nelements;

  // element geometry
  dfloat* elx[3] = {mesh->x, mesh->y, mesh->z};

  // element dimensions
  dlong n1[3] = {mesh->N+1, mesh->N+1, mesh->N+1};

  dlong m1[3] = {2*n1[0], 2*n1[1], 2*n1[2]};

  // used for # of cells in hash tables
  dlong hash_size = nelm*n1[0]*n1[1]*n1[2];

  MPI_Comm comm = platform_t::getInstance()->comm.mpiComm;

  findpts_handle = ogsFindptsSetup(3, comm, elx, n1, nelm, m1, bb_tol,
                                   hash_size, hash_size, npt_max, newton_tol,
                                   (occa::device*)&platform_t::getInstance()->device);
}


interp_t::~interp_t()
{
  ogsFindptsFree(handle->findpts);
}

void interp_t::findPoints(dfloat* x[], dlong xStride[],
                          dlong*  code,  dlong  codeStride,
                          dlong*  proc,  dlong  procStride,
                          dlong*  el,    dlong    elStride,
                          dfloat* r,     dlong     rStride,
                          dfloat* dist2, dlong dist2Stride,
                          dlong n,
                          bool printWarnings)
{
  // findpts takes strides in terms of bytes, but findPoints takes strides in terms of elements
  dlong xStrideBytes[3] = {xStride[0]*sizeof(dfloat),
                           xStride[1]*sizeof(dfloat),
                           xStride[2]*sizeof(dfloat)};

  ogsFindpts(code,   codeStride*sizeof(dlong),
             proc,   procStride*sizeof(dlong),
             el,       elStride*sizeof(dlong),
             r,         rStride*sizeof(dfloat),
             dist2, dist2Stride*sizeof(dfloat),
             x,     xStrideBytes,
             n, findpts);


  if(printWarnings) {
    dlong nFail = 0;
    for (int in = 0; in < n; ++in) {
      if (code[in] == 1) {
        if (dist2[in] > 10*newton_tol) {
          nFail += 1;
          if (nFail < 5){
            std::cerr << " WARNING: point on boundary or outside the mesh xy[z]d^2: "
                      << x[0][in] << "," << x[1][in] << ", " << x[2][in] << ", " << dist2[in] << std::endl;
          }
        }
      } else if (code[in] == 2) {
        nFail += 1;
        if (nFail < 5){
          std::cerr << " WARNING: point not within mesh xy[z]: "
                    << x[0][in] << "," << x[1][in] << ", " << x[2][in] << std::endl;
        }
      }
    }
    hlong counts[2] = {n, nFail};
    MPI_Reduce(MPI_IN_PLACE, counts, 2, MPI_HLONG, MPI_SUM, 0, platform_t::getInstance()->comm.mpiComm);
    if (platform_t::getInstance()->comm.mpiRank == 0 && counts[1] > 0) {
      std::cout << "interp::findPoints - Total number of points = " << counts[0]
                << ", failed = " << counts[1] << std::endl;
    }
  }
}

void interp_t::evalPoints(dfloat* fields, dlong nFields,
                          dlong*   code,  dlong  codeStride,
                          dlong*   proc,  dlong  procStride,
                          dlong*   el,    dlong    elStride,
                          dfloat*  r,     dlong     rStride,
                          dfloat** out,   dlong   outStride,
                          dlong n)
{
  dlong fieldOffset = nrs->fieldOffset;
  for (int i = 0; i < nFields; ++i) {
    ogsFindptsEval(out[i],  outStride[i]*sizeof(dfloat),
                   code,   codeStride   *sizeof(dlong),
                   proc,   procStride   *sizeof(dlong),
                   el,       elStride   *sizeof(dlong),
                   r,         rStride   *sizeof(dfloat),
                   n, fields+i*fieldOffset, findpts);
  }
}

void interp_t::evalPoints(occa::memory fields, dlong nFields,
                          dlong*   code,  dlong  codeStride,
                          dlong*   proc,  dlong  procStride,
                          dlong*   el,    dlong    elStride,
                          dfloat*  r,     dlong     rStride,
                          dfloat** out,   dlong   outStride,
                          dlong n)
{
  occa::memory o_fields = fields.cast(occa::dtype::get<dfloat>());
  for (int i = 0; i < nFields; ++i) {
    ogsFindptsEval(out[i],  outStride[i]*sizeof(dfloat),
                   code,   codeStride   *sizeof(dlong),
                   proc,   procStride   *sizeof(dlong),
                   el,       elStride   *sizeof(dlong),
                   r,         rStride   *sizeof(dfloat),
                   n, o_fields+i*fieldOffset, findpts);
  }
}


void interp_t::evalLocalPoints(dfloat fields, dlong nFields,
                               dlong*   el,    dlong    elStride,
                               dfloat*  r,     dlong     rStride,
                               dfloat** out,   dlong   outStride,
                               dlong n)
{
  if (n == 0 || nFields == 0) {
     return;
  }
  dlong fieldOffset = nrs->fieldOffset;
  for (int i = 0; i < nFields; ++i) {
    ogsFindptsLocalEval(out[i],  outStride[i]*sizeof(dfloat),
                        el,       elStride   *sizeof(dlong),
                        r,         rStride   *sizeof(dfloat),
                        n, fields+i*fieldOffset, findpts);
  }
}

void interp_t::evalLocalPoints(occa::memory fields, dlong nFields,
                               dlong*   el,    dlong    elStride,
                               dfloat*  r,     dlong     rStride,
                               dfloat** out,   dlong   outStride,
                               dlong n)
{
  if (n == 0 || nFields == 0) {
     return;
  }

  dlong  elStrideBytes =  elStride*sizeof(dlong);
  dlong   rStrideBytes =   rStride*sizeof(dfloat);
  dlong outStrideBytes = outStride*sizeof(dfloat);

  occa::device device = *findpts->device;
  dlong allocSize = (nfield*outStrideBytes+rStrideBytes+elStrideBytes)*pn;
  occa::memory workspace;
  occa::memory mempool = platform_t::getInstance()->o_mempool.o_ptr;
  if(alloc_size < mempool.size()) {
    workspace = mempool.cast(occa::dtype::byte);
  } else {
    workspace = device.malloc(allocSize, occa::dtype::byte);
  }
  occa::memory o_out = workspace; workspace += n*nField*sizeof(dfloat);
  occa::memory o_r   = workspace; workspace += n* rStrideBytes;
  occa::memory o_el  = workspace; workspace += n*elStrideBytes;
  o_r .copyFrom(r,   rStrideBytes*n);
  o_el.copyFrom(el, elStrideBytes*n);

  dlong fieldOffset = nrs->fieldOffset;
  for (int i = 0; i < nField; ++i) {
    occa::memory o_out_i = o_out.slice(i*sizeof(dfloat)*n, sizeof(dfloat)*n);
    ogsFindptsLocalEval(o_out_i, sizeof(dfloat),
                        o_el,     elStrideBytes,
                        o_r,       rStrideBytes,
                        n, o_fields+i*fieldOffset, findpts);
  }
  if(sizeof(dlong) == outStrideBytes){
    // combine d->h copies where able
    dlong i = 0;
    while (i < nField) {
      dlong start = out[i];
      dlong j = 0;
      while (i+j<nField && out[i+j] == start+j*outStrideBytes*n) ++j;
      o_out.copyTo(start, j*outStrideBytes*n);
      i += j;
    }
  } else {
    dfloat* outTemp = new dfloat[nField*sizeof(dfloat)*n];
    o_out.copyTo(outTemp, nField*sizeof(dfloat)*n);
    for(dlong i = 0; i < nFields; ++i) {
      for(dlong j = 0; j < n; ++j) {
        out[i][j*outStride] = outTemp[i*n + j];
      }
    }
  }
}

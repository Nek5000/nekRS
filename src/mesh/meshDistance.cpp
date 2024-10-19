#include "mesh.h"
#include "ogs.hpp"
#include "platform.hpp"

namespace
{

hlong hsum(mesh_t *mesh, const dlong N, occa::memory &o_a, MPI_Comm _comm)
{
  const auto blocksize = BLOCKSIZE;
  int Nblock = (N + blocksize - 1) / blocksize;
  const size_t Nbytes = Nblock * sizeof(hlong);

  static occa::memory o_scratch;
  static occa::memory h_scratch;
  static hlong *scratch;

  if (o_scratch.byte_size() < Nbytes) {
    o_scratch = platform->device.malloc(Nbytes);
    h_scratch = platform->device.mallocHost(Nbytes);
    scratch = (hlong *)h_scratch.ptr();
  }

  if (N > 1) {
    mesh->hlongSumKernel(Nblock, N, 0, o_a, o_scratch);
    o_scratch.copyTo(scratch, Nbytes);
  } else {
    o_a.copyTo(scratch, N);
  }

  hlong sum = 0;
  for (dlong n = 0; n < Nblock; ++n) {
    sum += scratch[n];
  }

  if (_comm != MPI_COMM_SELF) {
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_HLONG, MPI_SUM, _comm);
  }

  return sum;
}

// GPU-enabled mechanism of computing distances not yet supported
occa::memory
cheapDist(mesh_t *mesh, int nbID, const occa::memory &o_bID, dlong offsetFld, bool minDist, int maxIter)
{
  bool verbose = platform->options.compareArgs("VERBOSE", "TRUE");
  const auto [minCoord, maxCoord] = [&]() {
    const auto xMin = platform->linAlg->min(mesh->Nlocal, mesh->o_x, platform->comm.mpiComm);
    const auto yMin = platform->linAlg->min(mesh->Nlocal, mesh->o_y, platform->comm.mpiComm);
    const auto zMin = platform->linAlg->min(mesh->Nlocal, mesh->o_z, platform->comm.mpiComm);

    const auto xMax = platform->linAlg->max(mesh->Nlocal, mesh->o_x, platform->comm.mpiComm);
    const auto yMax = platform->linAlg->max(mesh->Nlocal, mesh->o_y, platform->comm.mpiComm);
    const auto zMax = platform->linAlg->max(mesh->Nlocal, mesh->o_z, platform->comm.mpiComm);

    return std::make_tuple(std::min({xMin, yMin, zMin}), std::max({xMax, yMax, zMax}));
  }();

  const auto Nfields = minDist ? 1 : nbID;
  auto o_dist = platform->device.malloc<dfloat>(Nfields * offsetFld);
  platform->linAlg->fill(Nfields * offsetFld, 0.0, o_dist);

  // as defined in nek5000 cheap_dist
  const auto largeNum = 10.0 * (maxCoord - minCoord);
  for (int fld = 0; fld < Nfields; ++fld) {
    auto o_dist_fld = o_dist + fld * offsetFld;
    platform->linAlg->fill(mesh->Nlocal, largeNum, o_dist_fld);
  }

  const auto zero = 0.0;

  if (minDist) {
    // zero-out distance across all boundaries specified in o_bid
    mesh->setBIDKernel(mesh->Nelements, 1, 0, nbID, zero, o_bID, 0, mesh->o_vmapM, mesh->o_EToB, o_dist);
  } else {
    // zero-out different distance for each boundary specified in o_bid
    mesh->setBIDKernel(mesh->Nelements,
                       nbID,
                       offsetFld,
                       nbID,
                       zero,
                       o_bID,
                       1,
                       mesh->o_vmapM,
                       mesh->o_EToB,
                       o_dist);
  }

  auto o_changed = platform->deviceMemoryPool.reserve<hlong>(mesh->Nlocal);

  for (int iter = 0; iter < maxIter; ++iter) {
    mesh->distanceKernel(mesh->Nelements,
                         Nfields,
                         offsetFld,
                         mesh->o_x,
                         mesh->o_y,
                         mesh->o_z,
                         o_dist,
                         o_changed);
    oogs::startFinish(o_dist, Nfields, offsetFld, ogsDfloat, ogsMin, mesh->oogs);

    const auto nchange = hsum(mesh, mesh->Nlocal, o_changed, platform->comm.mpiComm);

    // only compute dmax if verbose is true
    double dmax = 0.0;
    if (verbose) {
      dmax = platform->linAlg->amaxMany(mesh->Nlocal, Nfields, offsetFld, o_dist, platform->comm.mpiComm);
      if (platform->comm.mpiRank == 0) {
        std::cout << "distance: " << iter << " " << nchange << " " << dmax << std::endl;
      }
    }

    bool converged = (nchange == 0);
    if (converged) {
      break;
    }
  }

  o_changed.free();

  return o_dist;
}

} // namespace

occa::memory
mesh_t::distance(int nbID, const occa::memory &o_bID, dlong offsetFld, std::string type, int maxIter)
{
  lowerCase(type);

  occa::memory o_dist;
  if (type.find("cheap") != std::string::npos) {
    o_dist = cheapDist(this, nbID, o_bID, offsetFld, false, maxIter);
  } else {
    nekrsAbort(platform->comm.mpiComm,
               EXIT_FAILURE,
               "distance function type %s not supported!\n",
               type.c_str());
  }

  return o_dist;
}

occa::memory mesh_t::minDistance(int nbID, const occa::memory &o_bID, std::string type, int maxIter)
{
  lowerCase(type);

  occa::memory o_dist;
  if (type.find("cheap") != std::string::npos) {
    o_dist = cheapDist(this, nbID, o_bID, this->Nlocal, true, maxIter);
  } else {
    nekrsAbort(platform->comm.mpiComm,
               EXIT_FAILURE,
               "distance function type %s not supported!\n",
               type.c_str());
  }

  return o_dist;
}

std::vector<dfloat>
mesh_t::distance(const std::vector<dlong> &bID, dlong offsetFld, std::string type, int maxIter)
{
  auto o_bid = platform->deviceMemoryPool.reserve<dlong>(bID.size());
  o_bid.copyFrom(bID.data());

  auto o_dist = this->distance(bID.size(), o_bid, offsetFld, type, maxIter);

  std::vector<dfloat> dist(this->Nlocal);
  o_dist.copyTo(dist.data(), this->Nlocal);

  o_bid.free();
  o_dist.free();

  return dist;
}

std::vector<dfloat> mesh_t::minDistance(const std::vector<dlong> &bID, std::string type, int maxIter)
{
  auto o_bid = platform->deviceMemoryPool.reserve<dlong>(bID.size());
  o_bid.copyFrom(bID.data());

  auto o_dist = this->minDistance(bID.size(), o_bid, type, maxIter);

  std::vector<dfloat> dist(this->Nlocal);
  o_dist.copyTo(dist.data(), this->Nlocal);

  o_bid.free();
  o_dist.free();

  return dist;
}
#include "mesh.h"
#include "platform.hpp"

// interpolate to M-points
occa::memory mesh_t::intpMatrix(std::vector<dfloat> M) 
{ 
  nekrsCheck(M.size() > mesh_t::maxNqIntp, 
             MPI_COMM_SELF, 
             EXIT_FAILURE, 
             "%s\n", 
             "target N has to be smaller or equal to %d", mesh_t::maxNqIntp - 1);

  static std::array<occa::memory, mesh_t::maxNqIntp> o_J;
  if (o_J[M.size() - 1].isInitialized()) return o_J[M.size() - 1];

  std::vector<dfloat> J(this->Nq * M.size());
  InterpolationMatrix1D(this->N, this->Nq, this->r, M.size(), M.data(), J.data());

  auto transposeJ = [&]()
  {
    std::vector<dfloat> Jt(J.size());
    for (int i = 0; i < this->Nq; i++) {
      for (int j = 0; j < M.size(); j++) {
        Jt[i * M.size() + j] = J[j * this->Nq + i];
      }
    }
    return Jt;
  };

  o_J[M.size() - 1] = platform->device.malloc<dfloat>(J.size());
  o_J[M.size() - 1].copyFrom((M.size() < this->Nq) ? J.data() : transposeJ().data()); 

  return o_J[M.size() - 1];
}

void mesh_t::interpolate(mesh_t *mesh, const occa::memory& o_z, occa::memory& o_zM)
{
  std::vector<dfloat> M(mesh->Nq);
  for(int i = 0; i < M.size(); i++) M[i] = mesh->r[i];

  platform->linAlg->fill(mesh->Nlocal, 0.0, o_zM);

  const dlong nel = std::min(this->Nelements, mesh->Nelements);
  this->intpKernel[mesh->N](nel, intpMatrix(M), o_z, o_zM);
}

void mesh_t::map2Uniform(mesh_t *meshU, const occa::memory& o_z, occa::memory& o_zU)
{
  const auto Nu = meshU->N;

  auto U = [&]()
  {
    std::vector<dfloat> r(Nu + 1);
    r[0] = -1.0;
    r[Nu] = 1.0;
 
    const auto dr = (r[Nu] - r[0]) / Nu;
    for(int i = 1; i < Nu; i++) r[i] = r[i-1] + dr;
    return r;
  }();

  platform->linAlg->fill(this->Nelements * std::pow(U.size(), this->dim), 0.0, o_zU);
  this->intpKernel[Nu](this->Nelements, intpMatrix(U), o_z, o_zU);
}

void mesh_t::map2Uniform(const occa::memory& o_z, occa::memory& o_zU)
{
  map2Uniform(this, o_z, o_zU); 
}

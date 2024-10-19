#include <array>

#include "platform.hpp"
#include "avm.hpp"
#include "udf.hpp"

/**
 * Persson's artificial viscosity method (http://persson.berkeley.edu/pub/persson06shock.pdf) with P1
 * and https://www.mmnp-journal.org/articles/mmnp/pdf/2011/03/mmnp201163p57.pdf
 **/

namespace avm
{

occa::kernel relativeMassAveragedModeKernel;

occa::kernel computeMaxViscKernel;
occa::kernel interpolateP1Kernel;
occa::kernel modesKernel;

occa::memory o_vertexIds;
occa::memory o_r;
occa::memory o_s;
occa::memory o_t;

occa::memory o_modeMap;
occa::memory o_leastSquares1D;
occa::memory o_baseLineDecay;
occa::memory o_invVT;

oogs_t *gsh = nullptr;
mesh_t *mesh = nullptr;

namespace
{

occa::memory modeInfoKlocknerHex3D(int _N)
{
  const int _Np = (_N + 1) * (_N + 1) * (_N + 1);
  const int _Nmodes1D = (_N + 1);
  std::vector<int> _modeMap(_Np);

  int sk = 0, n = 0;
  for (int id = 0; id < _Nmodes1D; id++) {
    for (int j = 0; j < _Nmodes1D; j++) {
      for (int i = 0; i < _Nmodes1D; i++) {
        for (int k = 0; k < _Nmodes1D; k++) {
          if (std::max(std::max(i, j), k) == id) {
            _modeMap[n++] = sk;
          }
          sk++;
        }
      }
    }
    sk = 0;
  }

  auto o_modeMap = platform->device.malloc<int>(_modeMap.size());
  o_modeMap.copyFrom(_modeMap.data());
  return o_modeMap;
}

occa::memory leastSquaresFitKlockner(int _N)
{
  std::vector<dfloat> tmp(2 * _N);
  for (int n = 0; n < _N; n++) {
    tmp[2 * n + 0] = std::log10(n + 1);
    tmp[2 * n + 1] = 1.0;
  }

  std::vector<dfloat> _LSF(_N);
  auto invTmp = platform->linAlg->matrixPseudoInverse(2, tmp);
  for (int n = 0; n < _N; n++) {
    _LSF[n] = invTmp[n];
  }

  auto o_LSF = platform->device.malloc<dfloat>(_LSF.size());
  o_LSF.copyFrom(_LSF.data());
  return o_LSF;
}

occa::memory baseLineDecayKlockner(int _N)
{
  dfloat bsum = 0.0;
  for (int j = 1; j < _N + 1; j++) {
    bsum += 1.0 / std::pow(j, 2 * _N);
  }
  bsum = 1.0 / std::sqrt(bsum);

  std::vector<dfloat> _BLD(_N + 1, 0.0);
  for (int n = 1; n < _N + 1; n++) {
    const dfloat bdecay = bsum / std::pow(n, _N);
    _BLD[n] = bdecay * bdecay;
  }

  auto o_BLD = platform->device.malloc<dfloat>(_BLD.size());
  o_BLD.copyFrom(_BLD.data());
  return o_BLD;
}

} // namespace

void setup(mesh_t *mesh_, oogs_t *gsh_)
{
  mesh = mesh_;
  gsh = gsh_;

  o_vertexIds = platform->device.malloc<int>(mesh->Nverts, mesh->vertexNodes);
  o_r = platform->device.malloc<dfloat>(mesh->Np, mesh->r);
  o_s = platform->device.malloc<dfloat>(mesh->Np, mesh->s);
  o_t = platform->device.malloc<dfloat>(mesh->Np, mesh->t);

  o_modeMap = modeInfoKlocknerHex3D(mesh->N);
  o_leastSquares1D = leastSquaresFitKlockner(mesh->N);
  o_baseLineDecay = baseLineDecayKlockner(mesh->N);
  o_invVT = [&]() {
    std::vector<dfloat> V(mesh->Nq * mesh->Nq);
    Vandermonde1D(mesh->N, mesh->Nq, mesh->r, V.data());
    auto invV = platform->linAlg->matrixInverse(mesh->Nq, V);
    auto invVT = platform->linAlg->matrixTranspose(mesh->Nq, invV);

    auto o_out = platform->device.malloc<dfloat>(invVT.size());
    o_out.copyFrom(invVT.data());
    return o_out;
  }();

  std::string kernelName;

  kernelName = "relativeMassAveragedMode";
  relativeMassAveragedModeKernel = platform->kernelRequests.load(kernelName);

  kernelName = "computeMaxVisc";
  computeMaxViscKernel = platform->kernelRequests.load(kernelName);

  kernelName = "interpolateP1";
  interpolateP1Kernel = platform->kernelRequests.load(kernelName);

  kernelName = "core-tensorProduct1DHex3D";
  modesKernel = platform->kernelRequests.load(kernelName);
}

occa::memory viscosity(dlong UFieldOffset,
                       const occa::memory &o_U,
                       const occa::memory &o_S,
                       dfloat absTol,
                       dfloat scalingCoeff,
                       dfloat logS0,
                       dfloat kappa,
                       bool makeCont)
{
  auto o_nu = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nlocal);
  viscosity(UFieldOffset, o_U, o_S, o_nu, absTol, scalingCoeff, logS0, kappa, makeCont);
  return o_nu;
}

void viscosity(dlong UFieldOffset,
               const occa::memory &o_U,
               const occa::memory &o_S,
               occa::memory &o_nu,
               dfloat absTol,
               dfloat scalingCoeff,
               dfloat logS0,
               dfloat kappa,
               bool C0)
{
  occa::memory o_logSk = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nelements);
  occa::memory o_Shat = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nlocal);

  modesKernel(mesh->Nelements, o_invVT, o_S, o_Shat);
  relativeMassAveragedModeKernel(mesh->Nelements,
                                 absTol,
                                 o_modeMap,
                                 o_invVT,
                                 o_leastSquares1D,
                                 o_baseLineDecay,
                                 0,
                                 o_Shat,
                                 o_logSk);

  dfloat visMaxCoeff = 1.0;
  computeMaxViscKernel(mesh->Nelements,
                       UFieldOffset,
                       logS0,
                       kappa,
                       visMaxCoeff,
                       scalingCoeff,
                       mesh->o_x,
                       mesh->o_y,
                       mesh->o_z,
                       o_U,
                       o_logSk,
                       o_nu);

  if (C0) {
    oogs::startFinish(o_nu, 1, 0, ogsDfloat, ogsMax, gsh);
    interpolateP1Kernel(mesh->Nelements, o_vertexIds, o_r, o_s, o_t, o_nu);
  }
}

} // namespace avm

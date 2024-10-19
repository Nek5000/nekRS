#include "platform.hpp"
#include "mesh.h"

static const std::string section = "core-";
static const std::string suffix = "Hex3D";

namespace opSEM
{

void grad(mesh_t *mesh, dlong offset, const occa::memory &o_in, occa::memory &o_out)
{
  static occa::kernel kernel;
  if (!kernel.isInitialized()) {
    kernel = platform->kernelRequests.load(section + "wGradientVolume" + suffix);
  }
  kernel(mesh->Nelements, mesh->o_vgeo, mesh->o_D, offset, o_in, o_out);
}

occa::memory grad(mesh_t *mesh, dlong offset, const occa::memory &o_in)
{
  auto o_out = platform->deviceMemoryPool.reserve<dfloat>(mesh->dim * offset);
  grad(mesh, offset, o_in, o_out);
  return o_out;
}

void strongGrad(mesh_t *mesh, dlong offset, const occa::memory &o_in, occa::memory &o_out)
{
  static occa::kernel kernel;
  if (!kernel.isInitialized()) {
    kernel = platform->kernelRequests.load(section + "gradientVolume" + suffix);
  }
  kernel(mesh->Nelements, mesh->o_vgeo, mesh->o_D, offset, o_in, o_out);
}

occa::memory strongGrad(mesh_t *mesh, dlong offset, const occa::memory &o_in)
{
  auto o_out = platform->deviceMemoryPool.reserve<dfloat>(mesh->dim * offset);
  strongGrad(mesh, offset, o_in, o_out);
  return o_out;
}

void strongGradVec(mesh_t *mesh, dlong offset, const occa::memory &o_in, occa::memory &o_out)
{
  for (int i = 0; i < mesh->dim; i++) {
    auto o_u = o_in.slice(i * offset, mesh->Nlocal);
    auto o_grad_u = o_out.slice(i * mesh->dim * offset, mesh->dim * offset);
    strongGrad(mesh, offset, o_u, o_grad_u);
  }
}

occa::memory strongGradVec(mesh_t *mesh, dlong offset, const occa::memory &o_in)
{
  poolDeviceMemory<dfloat> o_out(mesh->dim * mesh->dim * offset);
  strongGradVec(mesh, offset, o_in, o_out);
  return o_out;
}

void divergence(mesh_t *mesh, dlong offset, const occa::memory &o_in, occa::memory &o_out)
{
  static occa::kernel kernel;
  if (!kernel.isInitialized()) {
    kernel = platform->kernelRequests.load(section + "wDivergenceVolume" + suffix);
  }
  kernel(mesh->Nelements, mesh->o_vgeo, mesh->o_D, offset, o_in, o_out);
}

occa::memory divergence(mesh_t *mesh, dlong offset, const occa::memory &o_in)
{
  auto o_out = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nlocal);
  divergence(mesh, offset, o_in, o_out);
  return o_out;
}

void strongDivergence(mesh_t *mesh, dlong offset, const occa::memory &o_in, occa::memory &o_out)
{
  static occa::kernel kernel;
  if (!kernel.isInitialized()) {
    kernel = platform->kernelRequests.load(section + "divergenceVolume" + suffix);
  }
  kernel(mesh->Nelements, mesh->o_vgeo, mesh->o_D, offset, o_in, o_out);
}

occa::memory strongDivergence(mesh_t *mesh, dlong offset, const occa::memory &o_in)
{
  auto o_out = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nlocal);
  strongDivergence(mesh, offset, o_in);
  return o_out;
}

void laplacian(mesh_t *mesh,
               dlong offset,
               const occa::memory &o_lambda,
               const occa::memory &o_in,
               occa::memory &o_out)
{
  static occa::memory o_fieldOffsetScan;
  static occa::kernel kernel;
  if (!kernel.isInitialized()) {
    kernel = platform->kernelRequests.load(section + "weakLaplacian" + suffix);
    o_fieldOffsetScan = platform->device.malloc<dlong>(1);
  }
  kernel(mesh->Nelements, 1, o_fieldOffsetScan, mesh->o_ggeo, mesh->o_D, o_lambda, o_in, o_out);
}

occa::memory laplacian(mesh_t *mesh, dlong offset, const occa::memory &o_lambda, const occa::memory &o_in)
{
  auto o_out = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nlocal);
  laplacian(mesh, offset, o_lambda, o_in, o_out);
  return o_out;
}

void strongLaplacian(mesh_t *mesh,
                     dlong offset,
                     const occa::memory &o_lambda,
                     const occa::memory &o_in,
                     occa::memory &o_out)
{
  auto o_grad = strongGrad(mesh, offset, o_in);
  oogs::startFinish(o_grad, mesh->dim, offset, ogsDfloat, ogsAdd, mesh->oogs);

  auto o_tmp = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nlocal);
  platform->linAlg->axmyz(mesh->Nlocal, 1.0, mesh->o_invAJw, o_lambda, o_tmp);
  platform->linAlg->axmyVector(mesh->Nlocal, offset, 0, 1.0, o_tmp, o_grad);

  o_out = strongDivergence(mesh, offset, o_grad);
}

occa::memory
strongLaplacian(mesh_t *mesh, dlong offset, const occa::memory &o_lambda, const occa::memory &o_in)
{
  occa::memory o_out;
  strongLaplacian(mesh, offset, o_lambda, o_in, o_out);
  return o_out;
}

void strongCurl(mesh_t *mesh, dlong offset, const occa::memory &o_in, occa::memory &o_out)
{
  static occa::kernel kernel;
  if (!kernel.isInitialized()) {
    kernel = platform->kernelRequests.load(section + "curl" + suffix);
  }
  const dlong scaleJW = 1;
  kernel(mesh->Nelements, scaleJW, mesh->o_vgeo, mesh->o_D, offset, o_in, o_out);
}

occa::memory strongCurl(mesh_t *mesh, dlong offset, const occa::memory &o_in)
{
  auto o_out = platform->deviceMemoryPool.reserve<dfloat>(mesh->dim * offset);
  strongCurl(mesh, offset, o_in, o_out);
  return o_out;
}

} // namespace opSEM

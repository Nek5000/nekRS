#include "mesh.h"
#include "linAlg.hpp"
#include "platform.hpp"
void mesh_t::move(){
  platform->timer.tic("meshUpdate", 1);
  // update o_x, o_y and o_z based on mesh->o_U using AB formula
  nStagesSumVectorKernel(
      Nelements * Np,
      fieldOffset,
      nAB,
      o_coeffAB,
      o_U,
      o_x,
      o_y,
      o_z
  );
  update();

  double flops = 6 * static_cast<double>(Nlocal) * nAB;
  platform->flopCounter->add("mesh_t::move", flops);
  platform->timer.toc("meshUpdate");
}
void mesh_t::update()
{
  geometricFactorsKernel(Nelements,
                         o_D,
                         o_gllw,
                         o_x,
                         o_y,
                         o_z,
                         o_LMM,
                         o_vgeo,
                         o_ggeo,
                         platform->o_mempool.slice0);

  double flopsGeometricFactors = 18 * Np * Nq + 91 * Np;
  flopsGeometricFactors *= static_cast<double>(Nelements);

  cubatureGeometricFactorsKernel(Nelements, o_cubD, o_x, o_y, o_z, o_cubInterpT, o_cubw, o_cubvgeo);

  double flopsCubatureGeometricFactors = 0.0;
  flopsCubatureGeometricFactors += 18 * Np * Nq;                                             // deriv
  flopsCubatureGeometricFactors += 18 * (cubNq * Np + cubNq * cubNq * Nq * Nq + cubNp * Nq); // c->f interp
  flopsCubatureGeometricFactors += 55 * cubNp; // geometric factor computation
  flopsCubatureGeometricFactors *= static_cast<double>(Nelements);

  const dfloat minJ =
      platform->linAlg->min(Nelements * Np, platform->o_mempool.slice0, platform->comm.mpiComm);
  const dfloat maxJ =
      platform->linAlg->max(Nelements * Np, platform->o_mempool.slice0, platform->comm.mpiComm);

  if (minJ < 0 || maxJ < 0) {
    if (platform->comm.mpiRank == 0) printf("Jacobian < 0 !!!");
  }  

  volume = platform->linAlg->sum(Nelements * Np, o_LMM, platform->comm.mpiComm);

  computeInvLMM();
  surfaceGeometricFactorsKernel(Nelements, o_gllw, o_faceNodes, o_vgeo, o_sgeo);

  double flopsSurfaceGeometricFactors = 32 * Nq * Nq;
  flopsSurfaceGeometricFactors *= static_cast<double>(Nelements);

  double flops = flopsGeometricFactors + flopsCubatureGeometricFactors + flopsSurfaceGeometricFactors;
  platform->flopCounter->add("mesh_t::update", flops);
}

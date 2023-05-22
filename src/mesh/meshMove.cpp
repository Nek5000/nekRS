#include "mesh.h"
#include "linAlg.hpp"
#include "platform.hpp"

void mesh_t::move()
{
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

void mesh_t::update(bool updateHost)
{
  geometricFactors();

  volume = platform->linAlg->sum(Nlocal, o_LMM, platform->comm.mpiComm);

  computeInvLMM();

  surfaceGeometricFactors();

  if (updateHost) {
    o_x.copyTo(x);
    o_y.copyTo(y);
    o_z.copyTo(z);
  }
}

#include "elliptic.h"
#include "platform.hpp"
#include <string>
#include "gslib.h"
#include "fem_amg_preco.hpp"
void ellipticSEMFEMSolve(elliptic_t* elliptic, occa::memory& o_r, occa::memory& o_z)
{
  mesh_t* mesh = elliptic->mesh;

  occa::memory& o_buffer = elliptic->o_SEMFEMBuffer1;
  occa::memory& o_buffer2 = elliptic->o_SEMFEMBuffer2;

  elliptic->preSEMFEMKernel(
    elliptic->numRowsSEMFEM,
    elliptic->o_dofMap,
    o_r,
    o_buffer
  );

  platform->linAlg->fill(elliptic->Nfields * elliptic->Ntotal, 0.0, o_z);

  // TODO: *NOT* device compatible
  hypre_solve((dfloat*)o_buffer2.ptr(), elliptic->hypreData, (dfloat*)o_buffer.ptr());

  elliptic->postSEMFEMKernel(
    elliptic->numRowsSEMFEM,
    elliptic->o_dofMap,
    o_buffer2,
    o_z
  );

  oogs::startFinish(o_z, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, ogsAdd, elliptic->oogs);
}
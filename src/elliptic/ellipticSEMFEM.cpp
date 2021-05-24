#include "elliptic.h"
#include "platform.hpp"
#include <string>
#include "gslib.h"
#include "boomerAMG.h"
#if 0
#include "amgx.h"
#endif
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

  if(elliptic->options.compareArgs("SEMFEM SOLVER", "BOOMERAMG")){
    // TODO: *NOT* device compatible
    boomerAMGSolve(o_buffer2.ptr(), o_buffer.ptr());
  } else {
    //AMGXsolve(o_buffer2.ptr(), o_buffer.ptr());
  }

  elliptic->postSEMFEMKernel(
    elliptic->numRowsSEMFEM,
    elliptic->o_dofMap,
    o_buffer2,
    o_z
  );

  oogs::startFinish(o_z, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, ogsAdd, elliptic->oogs);
}
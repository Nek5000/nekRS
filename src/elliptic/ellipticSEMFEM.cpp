#include "elliptic.h"
#include "platform.hpp"
#include <string>
#include "gslib.h"
extern "C" long long * get_row_start();
extern "C" long long * get_row_end();
extern "C" void fem_amg_solve(double*,double*);
void ellipticSEMFEMSolve(elliptic_t* elliptic, occa::memory& o_r, occa::memory& o_z)
{
  mesh_t* mesh = elliptic->mesh;

  occa::memory& o_buffer = elliptic->o_SEMFEMBuffer1;
  occa::memory& o_buffer2 = elliptic->o_SEMFEMBuffer2;

  long long rowStart = *get_row_start();
  long long rowEnd = *get_row_end();
  const dlong numRows = rowEnd - rowStart + 1;
  elliptic->preSEMFEMKernel(
    numRows,
    elliptic->o_dofMap,
    o_r,
    o_buffer
  );

  platform->linAlg->fill(elliptic->Nfields * elliptic->Ntotal, 0.0, o_z);

  // TODO: *NOT* device compatible
  fem_amg_solve((dfloat*)o_buffer2.ptr(), (dfloat*)o_buffer.ptr());

  elliptic->postSEMFEMKernel(
    numRows,
    elliptic->o_dofMap,
    o_buffer2,
    o_z
  );

  oogs::startFinish(o_z, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, ogsAdd, elliptic->oogs);
}
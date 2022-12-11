/*

   The MIT License (MIT)

   Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

 */

#include "elliptic.h"
#include "ellipticPrecon.h"
#include "timer.hpp"
#include "platform.hpp"
#include "linAlg.hpp"

void ellipticPreconditioner(elliptic_t *elliptic, occa::memory &o_r, occa::memory &o_z)
{
  mesh_t *mesh = elliptic->mesh;
  precon_t* precon = (precon_t*) elliptic->precon;
  setupAide &options = elliptic->options;

  const dlong Nlocal = mesh->Np * mesh->Nelements;

  occa::memory &o_rPfloat = elliptic->o_rPfloat;
  occa::memory &o_zPfloat = elliptic->o_zPfloat;

  platform->timer.tic(elliptic->name + " preconditioner", 1);
  if (options.compareArgs("PRECONDITIONER", "JACOBI")) {

    const pfloat one = 1.0;
    elliptic->axmyzManyPfloatKernel(Nlocal,
                                    elliptic->Nfields,
                                    elliptic->fieldOffset,
                                    one,
                                    o_r, /* dfloat */
                                    precon->o_invDiagA,
                                    o_z /* dfloat */
    );
    platform->flopCounter->add("jacobiPrecon", static_cast<double>(Nlocal) * elliptic->Nfields);
  }
  else if (options.compareArgs("PRECONDITIONER", "MULTIGRID")) {
    platform->linAlg->pfill(elliptic->fieldOffset * elliptic->Nfields, 0.0, o_zPfloat);
    platform->copyDfloatToPfloatKernel(elliptic->fieldOffset * elliptic->Nfields, o_r, o_rPfloat);
    precon->MGSolver->Run(o_rPfloat, o_zPfloat);
    platform->copyPfloatToDfloatKernel(elliptic->fieldOffset * elliptic->Nfields, o_zPfloat, o_z);
  }
  else if (options.compareArgs("PRECONDITIONER", "SEMFEM")) {
    platform->linAlg->pfill(elliptic->fieldOffset * elliptic->Nfields, 0.0, o_zPfloat);
    platform->copyDfloatToPfloatKernel(elliptic->fieldOffset * elliptic->Nfields, o_r, o_rPfloat);
    precon->SEMFEMSolver->run(o_rPfloat, o_zPfloat);
    platform->copyPfloatToDfloatKernel(elliptic->fieldOffset * elliptic->Nfields, o_zPfloat, o_z);
  }
  else if (options.compareArgs("PRECONDITIONER", "NONE")) {
    o_z.copyFrom(o_r, elliptic->fieldOffset * elliptic->Nfields * sizeof(dfloat));
  }
  else if (options.compareArgs("PRECONDITIONER", "USER")) {
    elliptic->userPreconditioner(o_r, o_z);
  }
  else {
    if (platform->comm.mpiRank == 0)
      printf("ERROR: Unknown preconditioner\n");
    MPI_Abort(platform->comm.mpiComm, 1);
  }
  platform->timer.toc(elliptic->name + " preconditioner");

  if (elliptic->allNeumann)
    ellipticZeroMean(elliptic, o_z);
}

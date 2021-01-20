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
#include "timer.hpp"
#include "linAlg.hpp"

void ellipticPreconditioner(elliptic_t* elliptic, occa::memory &o_r, occa::memory &o_z)
{
  mesh_t* mesh = elliptic->mesh;
  linAlg_t* linAlg = linAlg_t::getSingleton();
  precon_t* precon = elliptic->precon;
  setupAide options = elliptic->options;

  const dlong Nlocal = mesh->Np * mesh->Nelements;
  const dfloat one = 1.0;
  if(options.compareArgs("PRECONDITIONER", "JACOBI")) {
    linAlg->axmyzMany(Nlocal, elliptic->Nfields, elliptic->Ntotal, 1, one, o_r, precon->o_invDiagA, o_z);
  }else if (options.compareArgs("PRECONDITIONER", "MULTIGRID")) {
    timer::tic("preconditioner", 1);
    parAlmond::Precon(precon->parAlmond, o_z, o_r);
    timer::toc("preconditioner");
  }else {
    if(mesh->rank == 0) printf("ERRROR: Unknown preconditioner\n");
    MPI_Abort(mesh->comm, 1);
  }

  if(elliptic->allNeumann) // zero mean of RHS
    ellipticZeroMean(elliptic, o_z);
}

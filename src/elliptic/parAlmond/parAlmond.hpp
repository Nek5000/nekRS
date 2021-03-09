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

#ifndef PARALMOND_HPP
#define PARALMOND_HPP

#include <math.h>
#include <stdlib.h>
#include <cstring>

#include "nrssys.hpp"

#include "defines.hpp"
#include "utils.hpp"
#include "kernels.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "level.hpp"
#include "agmg.hpp"
#include "coarse.hpp"
#include "solver.hpp"


namespace parAlmond {

solver_t *Init(occa::device device, MPI_Comm comm, setupAide options);

void AMGSetup(solver_t* M,
             hlong* rowStarts,
             dlong nnz,
             hlong* Ai,
             hlong* Aj,
             dfloat* Avals,
             bool nullSpace,
             dfloat nullSpacePenalty);

void Precon(solver_t* M, occa::memory o_x, occa::memory o_rhs);

void Report(solver_t *M);

void Free(solver_t* M);

} //namespace parAlmond

#endif

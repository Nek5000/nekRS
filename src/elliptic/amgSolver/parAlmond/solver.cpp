/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#include "parAlmond.hpp"
#include "platform.hpp"

namespace parAlmond {

solver_t::solver_t(occa::device device_, MPI_Comm comm_,
                   setupAide options_) {

  device = device_;

  comm = comm_;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  levels = (multigridLevel **) calloc(MAX_LEVELS,sizeof(multigridLevel *));
  numLevels = 0;

  options = options_;

  if(options.compareArgs("PARALMOND CYCLE", "VCYCLE")) {
    ctype = VCYCLE;
    additive = false;
    if(options.compareArgs("PARALMOND CYCLE", "ADDITIVE")) {
      if (options.compareArgs("PARALMOND SMOOTHER", "CHEBYSHEV")) {
        if(rank==0) printf("Additive vcycle is not supported for Chebyshev acceleration!\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }
      additive = true;
      overlapCrsGridSolve = false;
      if(options.compareArgs("PARALMOND CYCLE", "OVERLAPCRS")){
        if(platform->device.mode() == "Serial" || platform->device.mode() == "OpenMP"){
          overlapCrsGridSolve = false;
        } else {
          overlapCrsGridSolve = true;
          int provided;
          MPI_Query_thread(&provided);
          if(provided != MPI_THREAD_MULTIPLE) {
            overlapCrsGridSolve = false;
            if(rank ==0 && size > 1) printf("disable overlapping coarse solve as MPI_THREAD_MULTIPLE is not supported!\n");
          }
          if(size == 1) overlapCrsGridSolve = true;
        }
        if(rank ==0 && overlapCrsGridSolve) printf("overlapping coarse grid solve enabled\n");
      }
    } else {
      if (options.compareArgs("PARALMOND SMOOTHER", "RAS") || options.compareArgs("PARALMOND SMOOTHER", "ASM")) {
        if(!options.compareArgs("PARALMOND SMOOTHER", "CHEBYSHEV")){
          if(rank==0) 
            printf("Multiplicative vcycle is not supported for RAS/ASM smoother without Chebycehv acceleration!\n");
          MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); 
        }
      }
    }
  } else {
    if(rank==0) printf("Unknown multigrid cycle type!\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); 
  }
}

solver_t::~solver_t() {

  for (int n=0;n<numLevels;n++)
    delete levels[n];

  free(levels);
}

void solver_t::Report() {

}

}

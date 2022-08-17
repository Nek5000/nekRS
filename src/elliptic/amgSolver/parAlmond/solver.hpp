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

#ifndef PARALMOND_SOLVER_HPP
#define PARALMOND_SOLVER_HPP

namespace parAlmond {

class solver_t {

public:
  MPI_Comm comm;
  int rank, size;

  occa::device device;
  setupAide options;

  CycleType ctype;
  SmoothType stype;

  int numLevels;
  int AMGstartLev, baseLevel;
  multigridLevel **levels=NULL;

  coarseSolver *coarseLevel;

  bool additive;
  bool overlapCrsGridSolve;

  solver_t(occa::device otherdevice, MPI_Comm othercomm,
           setupAide otheroptions);

  ~solver_t();

  void Report();

  void device_vcycle(int k);
  void additiveVcycle();
};

} //namespace parAlmond

#endif

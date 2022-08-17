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

#ifndef PARALMOND_LEVEL_HPP
#define PARALMOND_LEVEL_HPP

namespace parAlmond {

class multigridLevel {

public:
  dlong Nrows, Ncols;

  occa::memory o_x, o_rhs, o_res;

  SmoothType stype;

  MPI_Comm comm;

  multigridLevel(dlong N, dlong M, MPI_Comm comm);
  ~multigridLevel();

  virtual void Ax(occa::memory o_x, occa::memory o_Ax)=0;

  virtual void smooth(occa::memory o_rhs, occa::memory o_x, bool x_is_zero)=0;

  virtual void residual(occa::memory o_rhs, occa::memory o_x, occa::memory o_res)=0;

  virtual void coarsen(occa::memory o_x, occa::memory o_Cx)=0;

  virtual void prolongate(occa::memory o_x, occa::memory o_Px)=0;

  virtual void Report()=0;
};

}

#endif

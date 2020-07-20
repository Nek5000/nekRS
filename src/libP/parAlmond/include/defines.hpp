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

#ifndef PARALMOND_DEFINES_HPP
#define PARALMOND_DEFINES_HPP

#define BLOCKSIZE 512
#define NBLOCKS 128

#define MAX_LEVELS 100
#define GPU_CPU_SWITCH_SIZE 0 //host-device switch threshold

#define NUMKCYCLES 3
#define COARSENTHREASHOLD 0.2
#define KCYCLETOL 0.2

namespace parAlmond {

extern int ChebyshevIterations;

typedef enum {VCYCLE=0,KCYCLE=1,EXACT=3} CycleType;
typedef enum {PCG=0,GMRES=1} KrylovType;
typedef enum {JACOBI=0,DAMPED_JACOBI=1,CHEBYSHEV=2} SmoothType;

} //namespace parAlmond

#endif

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

template < int p_Nq >
void ellipticSerialSetScalarKernel(const hlong Nelements,
                                   const dfloat alpha,
                                   dfloat* __restrict__ cpu_a)
{
  cpu_a = (dfloat*)__builtin_assume_aligned(cpu_a, USE_OCCA_MEM_BYTE_ALIGN);

#define p_Np (p_Nq * p_Nq * p_Nq)

  for(hlong e = 0; e < Nelements; ++e)
    for(int i = 0; i < p_Np; ++i)
      cpu_a[e * p_Np + i] = alpha;

#undef p_Np
}

void ellipticSerialSetScalar(const int Nq,
                             const hlong Nelements,
                             const dfloat alpha,
                             occa::memory &o_a)
{
  dfloat* __restrict__ cpu_a =
    (dfloat*)__builtin_assume_aligned(o_a.ptr(), USE_OCCA_MEM_BYTE_ALIGN);

  switch(Nq) {
  case  2: ellipticSerialSetScalarKernel <  2 > (Nelements, alpha, cpu_a);
    break;
  case  3: ellipticSerialSetScalarKernel <  3 > (Nelements, alpha, cpu_a);
    break;
  case  4: ellipticSerialSetScalarKernel <  4 > (Nelements, alpha, cpu_a);
    break;
  case  5: ellipticSerialSetScalarKernel <  5 > (Nelements, alpha, cpu_a);
    break;
  case  6: ellipticSerialSetScalarKernel <  6 > (Nelements, alpha, cpu_a);
    break;
  case  7: ellipticSerialSetScalarKernel <  7 > (Nelements, alpha, cpu_a);
    break;
  case  8: ellipticSerialSetScalarKernel <  8 > (Nelements, alpha, cpu_a);
    break;
  case  9: ellipticSerialSetScalarKernel <  9 > (Nelements, alpha, cpu_a);
    break;
  case 10: ellipticSerialSetScalarKernel < 10 > (Nelements, alpha, cpu_a);
    break;
  }
}

template < int p_Nq, int Nfields >
void ellipticSerialBlockSetScalarKernel(const hlong Nelements,
                                        const hlong offset,
                                        const dfloat alpha,
                                        dfloat* __restrict__ cpu_a)
{
  cpu_a = (dfloat*)__builtin_assume_aligned(cpu_a, USE_OCCA_MEM_BYTE_ALIGN);

#define p_Np (p_Nq * p_Nq * p_Nq)

  for(hlong e = 0; e < Nelements; ++e)
    for(int i = 0; i < p_Np; ++i)
      for(int fld = 0; fld < Nfields; fld++)
        cpu_a[e * p_Np + i + fld * offset] = alpha;

#undef p_Np
}

void ellipticSerialSetScalar(const int Nq,
                             const hlong Nelements,
                             const int Nfields,
                             const hlong offset,
                             const dfloat alpha,
                             occa::memory &o_a)
{
  dfloat* __restrict__ cpu_a =
    (dfloat*)__builtin_assume_aligned(o_a.ptr(), USE_OCCA_MEM_BYTE_ALIGN);
  if(Nfields == 2) {
    switch(Nq) {
    case  2: ellipticSerialBlockSetScalarKernel <  2, 2 > (Nelements, offset, alpha, cpu_a);
      break;
    case  3: ellipticSerialBlockSetScalarKernel <  3, 2 > (Nelements, offset, alpha, cpu_a);
      break;
    case  4: ellipticSerialBlockSetScalarKernel <  4, 2 > (Nelements, offset, alpha, cpu_a);
      break;
    case  5: ellipticSerialBlockSetScalarKernel <  5, 2 > (Nelements, offset, alpha, cpu_a);
      break;
    case  6: ellipticSerialBlockSetScalarKernel <  6, 2 > (Nelements, offset, alpha, cpu_a);
      break;
    case  7: ellipticSerialBlockSetScalarKernel <  7, 2 > (Nelements, offset, alpha, cpu_a);
      break;
    case  8: ellipticSerialBlockSetScalarKernel <  8, 2 > (Nelements, offset, alpha, cpu_a);
      break;
    case  9: ellipticSerialBlockSetScalarKernel <  9, 2 > (Nelements, offset, alpha, cpu_a);
      break;
    case 10: ellipticSerialBlockSetScalarKernel < 10, 2 > (Nelements, offset, alpha, cpu_a);
      break;
    }
  }else if(Nfields == 3) {
    switch(Nq) {
    case  2: ellipticSerialBlockSetScalarKernel <  2, 3 > (Nelements, offset, alpha, cpu_a);
      break;
    case  3: ellipticSerialBlockSetScalarKernel <  3, 3 > (Nelements, offset, alpha, cpu_a);
      break;
    case  4: ellipticSerialBlockSetScalarKernel <  4, 3 > (Nelements, offset, alpha, cpu_a);
      break;
    case  5: ellipticSerialBlockSetScalarKernel <  5, 3 > (Nelements, offset, alpha, cpu_a);
      break;
    case  6: ellipticSerialBlockSetScalarKernel <  6, 3 > (Nelements, offset, alpha, cpu_a);
      break;
    case  7: ellipticSerialBlockSetScalarKernel <  7, 3 > (Nelements, offset, alpha, cpu_a);
      break;
    case  8: ellipticSerialBlockSetScalarKernel <  8, 3 > (Nelements, offset, alpha, cpu_a);
      break;
    case  9: ellipticSerialBlockSetScalarKernel <  9, 3 > (Nelements, offset, alpha, cpu_a);
      break;
    case 10: ellipticSerialBlockSetScalarKernel < 10, 3 > (Nelements, offset, alpha, cpu_a);
      break;
    }
  }else if(Nfields == 4) {
    switch(Nq) {
    case  2: ellipticSerialBlockSetScalarKernel <  2, 4 > (Nelements, offset, alpha, cpu_a);
      break;
    case  3: ellipticSerialBlockSetScalarKernel <  3, 4 > (Nelements, offset, alpha, cpu_a);
      break;
    case  4: ellipticSerialBlockSetScalarKernel <  4, 4 > (Nelements, offset, alpha, cpu_a);
      break;
    case  5: ellipticSerialBlockSetScalarKernel <  5, 4 > (Nelements, offset, alpha, cpu_a);
      break;
    case  6: ellipticSerialBlockSetScalarKernel <  6, 4 > (Nelements, offset, alpha, cpu_a);
      break;
    case  7: ellipticSerialBlockSetScalarKernel <  7, 4 > (Nelements, offset, alpha, cpu_a);
      break;
    case  8: ellipticSerialBlockSetScalarKernel <  8, 4 > (Nelements, offset, alpha, cpu_a);
      break;
    case  9: ellipticSerialBlockSetScalarKernel <  9, 4 > (Nelements, offset, alpha, cpu_a);
      break;
    case 10: ellipticSerialBlockSetScalarKernel < 10, 4 > (Nelements, offset, alpha, cpu_a);
      break;
    }
  }
}

void ellipticSetScalar(elliptic_t* elliptic, dfloat alpha, occa::memory &o_a)
{
  mesh_t* mesh = elliptic->mesh;
  setupAide &options = elliptic->options;

  int continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");
  int serial = options.compareArgs("THREAD MODEL", "SERIAL");

  dlong Nlocal = mesh->Nelements * mesh->Np;

  // a[n] = alpha n\in [0,Ntotal)
  if(serial == 1 && continuous == 1) {
    if(elliptic->blockSolver)
      ellipticSerialSetScalar(mesh->Nq,
                              mesh->Nelements,
                              elliptic->Nfields,
                              elliptic->Ntotal,
                              alpha,
                              o_a);
    else
      ellipticSerialSetScalar(mesh->Nq, mesh->Nelements, alpha, o_a);
    return;
  }

  // if not Serial
  if(elliptic->blockSolver)
    elliptic->setScalarKernel(Nlocal, elliptic->Ntotal, alpha, o_a);
  else
    elliptic->setScalarKernel(Nlocal, alpha, o_a);
}

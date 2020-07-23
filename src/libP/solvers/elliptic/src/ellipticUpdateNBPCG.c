
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
dfloat ellipticSerialUpdate1NBPCGKernel(const hlong Nelements, const int useWeight,
                                        const dfloat* __restrict__ cpu_invDegree,
                                        const dfloat* __restrict__ cpu_z,
                                        const dfloat* __restrict__ cpu_Z,
                                        const dfloat beta,
                                        dfloat* __restrict__ cpu_p,
                                        dfloat* __restrict__ cpu_s)
{
#define p_Np (p_Nq * p_Nq * p_Nq)

  cpu_p  = (dfloat*)__builtin_assume_aligned(cpu_p,  USE_OCCA_MEM_BYTE_ALIGN);
  cpu_s  = (dfloat*)__builtin_assume_aligned(cpu_s,  USE_OCCA_MEM_BYTE_ALIGN);
  cpu_z  = (dfloat*)__builtin_assume_aligned(cpu_z,  USE_OCCA_MEM_BYTE_ALIGN);
  cpu_Z  = (dfloat*)__builtin_assume_aligned(cpu_Z,  USE_OCCA_MEM_BYTE_ALIGN);

  cpu_invDegree = (dfloat*)__builtin_assume_aligned(cpu_invDegree,  USE_OCCA_MEM_BYTE_ALIGN);

  dfloat pdots = 0;

  for(hlong e = 0; e < Nelements; ++e)
    for(int i = 0; i < p_Np; ++i) {
      const hlong n = e * p_Np + i;

      dfloat pn = cpu_p[n];
      dfloat zn = cpu_z[n];
      dfloat sn = cpu_s[n];
      dfloat Zn = cpu_Z[n];

      pn = zn + beta * pn;
      sn = Zn + beta * sn;

      dfloat invDeg = (useWeight) ? cpu_invDegree[n]: 1.0;

      pdots += pn * sn * invDeg;

      cpu_p[n] = pn;
      cpu_s[n] = sn;
    }

#undef p_Np

  return pdots;
}

template < int p_Nq >
dfloat ellipticBlockSerialUpdate1NBPCGKernel(const int Nfields, const hlong offset,
                                             const hlong Nelements, const int useWeight,
                                             const dfloat* __restrict__ cpu_invDegree,
                                             const dfloat* __restrict__ cpu_z,
                                             const dfloat* __restrict__ cpu_Z,
                                             const dfloat beta,
                                             dfloat* __restrict__ cpu_p,
                                             dfloat* __restrict__ cpu_s)
{
#define p_Np (p_Nq * p_Nq * p_Nq)

  cpu_p  = (dfloat*)__builtin_assume_aligned(cpu_p,  USE_OCCA_MEM_BYTE_ALIGN);
  cpu_s  = (dfloat*)__builtin_assume_aligned(cpu_s,  USE_OCCA_MEM_BYTE_ALIGN);
  cpu_z  = (dfloat*)__builtin_assume_aligned(cpu_z,  USE_OCCA_MEM_BYTE_ALIGN);
  cpu_Z  = (dfloat*)__builtin_assume_aligned(cpu_Z,  USE_OCCA_MEM_BYTE_ALIGN);

  cpu_invDegree = (dfloat*)__builtin_assume_aligned(cpu_invDegree,  USE_OCCA_MEM_BYTE_ALIGN);

  dfloat pdots = 0;
  for(int fld = 0; fld < Nfields; fld++)
    for(hlong e = 0; e < Nelements; ++e)
      for(int i = 0; i < p_Np; ++i) {
        const hlong n = e * p_Np + i + fld * offset;

        dfloat pn = cpu_p[n];
        dfloat zn = cpu_z[n];
        dfloat sn = cpu_s[n];
        dfloat Zn = cpu_Z[n];

        pn = zn + beta * pn;
        sn = Zn + beta * sn;

        dfloat invDeg = (useWeight) ? cpu_invDegree[n]: 1.0;

        pdots += pn * sn * invDeg;

        cpu_p[n] = pn;
        cpu_s[n] = sn;
      }

#undef p_Np

  return pdots;
}

dfloat ellipticSerialUpdate1NBPCG(const int Nq,
                                  const hlong Nelements,
                                  const int useWeight,
                                  occa::memory &o_invDegree,
                                  occa::memory &o_z,
                                  occa::memory &o_Z,
                                  const dfloat beta,
                                  occa::memory &o_p,
                                  occa::memory &o_s)
{
  const dfloat* __restrict__ cpu_z  = (dfloat*)__builtin_assume_aligned(o_z.ptr(),
                                                                        USE_OCCA_MEM_BYTE_ALIGN);
  const dfloat* __restrict__ cpu_Z  = (dfloat*)__builtin_assume_aligned(o_Z.ptr(),
                                                                        USE_OCCA_MEM_BYTE_ALIGN);
  dfloat* __restrict__ cpu_p  =
    (dfloat*)__builtin_assume_aligned(o_p.ptr(), USE_OCCA_MEM_BYTE_ALIGN);
  dfloat* __restrict__ cpu_s  =
    (dfloat*)__builtin_assume_aligned(o_s.ptr(), USE_OCCA_MEM_BYTE_ALIGN);

  const dfloat* __restrict__ cpu_invDegree = (dfloat*)__builtin_assume_aligned(
    o_invDegree.ptr(),
    USE_OCCA_MEM_BYTE_ALIGN);

  dfloat pdots = 0;

  switch(Nq) {
  case  2: pdots = ellipticSerialUpdate1NBPCGKernel <  2 >
                   (Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s);
    break;
  case  3: pdots = ellipticSerialUpdate1NBPCGKernel <  3 >
                   (Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s);
    break;
  case  4: pdots = ellipticSerialUpdate1NBPCGKernel <  4 >
                   (Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s);
    break;
  case  5: pdots = ellipticSerialUpdate1NBPCGKernel <  5 >
                   (Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s);
    break;
  case  6: pdots = ellipticSerialUpdate1NBPCGKernel <  6 >
                   (Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s);
    break;
  case  7: pdots = ellipticSerialUpdate1NBPCGKernel <  7 >
                   (Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s);
    break;
  case  8: pdots = ellipticSerialUpdate1NBPCGKernel <  8 >
                   (Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s);
    break;
  case  9: pdots = ellipticSerialUpdate1NBPCGKernel <  9 >
                   (Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s);
    break;
  case 10: pdots = ellipticSerialUpdate1NBPCGKernel < 10 >
                   (Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s);
    break;
  case 11: pdots = ellipticSerialUpdate1NBPCGKernel < 11 >
                   (Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s);
    break;
  case 12: pdots = ellipticSerialUpdate1NBPCGKernel < 12 >
                   (Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p, cpu_s);
    break;
  }

  return pdots;
}

dfloat ellipticBlockSerialUpdate1NBPCG(const int Nfields,
                                       const hlong offset,
                                       const int Nq,
                                       const hlong Nelements,
                                       const int useWeight,
                                       occa::memory &o_invDegree,
                                       occa::memory &o_z,
                                       occa::memory &o_Z,
                                       const dfloat beta,
                                       occa::memory &o_p,
                                       occa::memory &o_s)
{
  const dfloat* __restrict__ cpu_z  = (dfloat*)__builtin_assume_aligned(o_z.ptr(),
                                                                        USE_OCCA_MEM_BYTE_ALIGN);
  const dfloat* __restrict__ cpu_Z  = (dfloat*)__builtin_assume_aligned(o_Z.ptr(),
                                                                        USE_OCCA_MEM_BYTE_ALIGN);
  dfloat* __restrict__ cpu_p  =
    (dfloat*)__builtin_assume_aligned(o_p.ptr(), USE_OCCA_MEM_BYTE_ALIGN);
  dfloat* __restrict__ cpu_s  =
    (dfloat*)__builtin_assume_aligned(o_s.ptr(), USE_OCCA_MEM_BYTE_ALIGN);

  const dfloat* __restrict__ cpu_invDegree = (dfloat*)__builtin_assume_aligned(
    o_invDegree.ptr(),
    USE_OCCA_MEM_BYTE_ALIGN);

  dfloat pdots = 0;

  switch(Nq) {
  case  2: pdots = ellipticBlockSerialUpdate1NBPCGKernel <  2 >
                   (Nfields, offset, Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p,
                    cpu_s);
    break;
  case  3: pdots = ellipticBlockSerialUpdate1NBPCGKernel <  3 >
                   (Nfields, offset, Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p,
                    cpu_s);
    break;
  case  4: pdots = ellipticBlockSerialUpdate1NBPCGKernel <  4 >
                   (Nfields, offset, Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p,
                    cpu_s);
    break;
  case  5: pdots = ellipticBlockSerialUpdate1NBPCGKernel <  5 >
                   (Nfields, offset, Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p,
                    cpu_s);
    break;
  case  6: pdots = ellipticBlockSerialUpdate1NBPCGKernel <  6 >
                   (Nfields, offset, Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p,
                    cpu_s);
    break;
  case  7: pdots = ellipticBlockSerialUpdate1NBPCGKernel <  7 >
                   (Nfields, offset, Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p,
                    cpu_s);
    break;
  case  8: pdots = ellipticBlockSerialUpdate1NBPCGKernel <  8 >
                   (Nfields, offset, Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p,
                    cpu_s);
    break;
  case  9: pdots = ellipticBlockSerialUpdate1NBPCGKernel <  9 >
                   (Nfields, offset, Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p,
                    cpu_s);
    break;
  case 10: pdots = ellipticBlockSerialUpdate1NBPCGKernel < 10 >
                   (Nfields, offset, Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p,
                    cpu_s);
    break;
  case 11: pdots = ellipticBlockSerialUpdate1NBPCGKernel < 11 >
                   (Nfields, offset, Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p,
                    cpu_s);
    break;
  case 12: pdots = ellipticBlockSerialUpdate1NBPCGKernel < 12 >
                   (Nfields, offset, Nelements, useWeight, cpu_invDegree, cpu_z, cpu_Z, beta, cpu_p,
                    cpu_s);
    break;
  }

  return pdots;
}

void ellipticNonBlockingUpdate1NBPCG(elliptic_t* elliptic,
                                     occa::memory &o_z, occa::memory &o_Z, const dfloat beta,
                                     occa::memory &o_p, occa::memory &o_s,
                                     dfloat* localpdots, dfloat* globalpdots, MPI_Request* request)
{
  setupAide &options = elliptic->options;

  int fixedIterationCountFlag = 0;
  int enableGatherScatters = 1;
  int enableReductions = 1;
  int serial = options.compareArgs("THREAD MODEL", "SERIAL");
  int continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");

  options.getArgs("DEBUG ENABLE REDUCTIONS", enableReductions);
  options.getArgs("DEBUG ENABLE OGS", enableGatherScatters);
  if(options.compareArgs("FIXED ITERATION COUNT", "TRUE"))
    fixedIterationCountFlag = 1;

  // register scalars

  mesh_t* mesh = elliptic->mesh;
  const dlong Nlocal = mesh->Np * mesh->Nelements;

  int useWeight = (continuous == 1);

  if(serial == 1) {
    if(elliptic->blockSolver)
      localpdots[0] = ellipticBlockSerialUpdate1NBPCG(elliptic->Nfields, elliptic->Ntotal,
                                                      mesh->Nq, mesh->Nelements, useWeight,
                                                      elliptic->o_invDegree,
                                                      o_z, o_Z, beta, o_p, o_s);
    else
      localpdots[0] = ellipticSerialUpdate1NBPCG(mesh->Nq, mesh->Nelements, useWeight,
                                                 elliptic->o_invDegree,
                                                 o_z, o_Z, beta, o_p, o_s);
  }else {
    // p <= z + beta*p
    // s <= Z + beta*s
    // dot(p,s)
    if(elliptic->blockSolver)
      elliptic->update1NBPCGKernel(Nlocal,
                                   elliptic->Ntotal,
                                   elliptic->NblocksUpdatePCG,
                                   useWeight,
                                   elliptic->o_invDegree,
                                   o_z,
                                   o_Z,
                                   beta,
                                   o_p,
                                   o_s,
                                   elliptic->o_tmppdots);
    else
      elliptic->update1NBPCGKernel(mesh->Nelements * mesh->Np,
                                   elliptic->NblocksUpdatePCG,
                                   useWeight,
                                   elliptic->o_invDegree,
                                   o_z,
                                   o_Z,
                                   beta,
                                   o_p,
                                   o_s,
                                   elliptic->o_tmppdots);

    elliptic->o_tmppdots.copyTo(elliptic->tmppdots);

    *localpdots = 0;
    for(int n = 0; n < elliptic->NblocksUpdatePCG; ++n)
      *localpdots += elliptic->tmppdots[n];
  }

  *globalpdots = 0;
  if(enableReductions)
    MPI_Iallreduce(localpdots, globalpdots, 1, MPI_DFLOAT, MPI_SUM, mesh->comm, request);
  else
    *globalpdots = 1;
}

template < int p_Nq >
void ellipticSerialUpdate2NBPCGKernel(const hlong Nelements, const int useWeight,
                                      const dfloat* __restrict__ cpu_invDegree,
                                      const dfloat* __restrict__ cpu_s,
                                      const dfloat* __restrict__ cpu_S,
                                      const dfloat alpha,
                                      dfloat* __restrict__ cpu_r,
                                      dfloat* __restrict__ cpu_z,
                                      dfloat* localdots)
{
#define p_Np (p_Nq * p_Nq * p_Nq)

  cpu_r  = (dfloat*)__builtin_assume_aligned(cpu_r,  USE_OCCA_MEM_BYTE_ALIGN);
  cpu_z  = (dfloat*)__builtin_assume_aligned(cpu_z,  USE_OCCA_MEM_BYTE_ALIGN);
  cpu_s  = (dfloat*)__builtin_assume_aligned(cpu_s,  USE_OCCA_MEM_BYTE_ALIGN);
  cpu_S  = (dfloat*)__builtin_assume_aligned(cpu_S,  USE_OCCA_MEM_BYTE_ALIGN);

  cpu_invDegree = (dfloat*)__builtin_assume_aligned(cpu_invDegree,  USE_OCCA_MEM_BYTE_ALIGN);

  dfloat rdotz = 0;
  dfloat zdotz = 0;
  dfloat rdotr = 0;

  for(hlong e = 0; e < Nelements; ++e)
    for(int i = 0; i < p_Np; ++i) {
      const hlong n = e * p_Np + i;

      dfloat rn = cpu_r[n];
      dfloat zn = cpu_z[n];
      dfloat sn = cpu_s[n];
      dfloat Sn = cpu_S[n];

      dfloat invDeg = (useWeight) ? cpu_invDegree[n]: 1.0;

      rn = rn - alpha * sn;
      zn = zn - alpha * Sn;

      rdotz += rn * zn * invDeg;
      zdotz += zn * zn * invDeg;
      rdotr += rn * rn * invDeg;

      cpu_r[n] = rn;
      cpu_z[n] = zn;
    }

  localdots[0] = rdotz;
  localdots[1] = zdotz;
  localdots[2] = rdotr;

#undef p_Np
}

template < int p_Nq >
void ellipticBlockSerialUpdate2NBPCGKernel(const int Nfields, const hlong offset,
                                           const hlong Nelements, const int useWeight,
                                           const dfloat* __restrict__ cpu_invDegree,
                                           const dfloat* __restrict__ cpu_s,
                                           const dfloat* __restrict__ cpu_S,
                                           const dfloat alpha,
                                           dfloat* __restrict__ cpu_r,
                                           dfloat* __restrict__ cpu_z,
                                           dfloat* localdots)
{
#define p_Np (p_Nq * p_Nq * p_Nq)

  cpu_r  = (dfloat*)__builtin_assume_aligned(cpu_r,  USE_OCCA_MEM_BYTE_ALIGN);
  cpu_z  = (dfloat*)__builtin_assume_aligned(cpu_z,  USE_OCCA_MEM_BYTE_ALIGN);
  cpu_s  = (dfloat*)__builtin_assume_aligned(cpu_s,  USE_OCCA_MEM_BYTE_ALIGN);
  cpu_S  = (dfloat*)__builtin_assume_aligned(cpu_S,  USE_OCCA_MEM_BYTE_ALIGN);

  cpu_invDegree = (dfloat*)__builtin_assume_aligned(cpu_invDegree,  USE_OCCA_MEM_BYTE_ALIGN);

  dfloat rdotz = 0;
  dfloat zdotz = 0;
  dfloat rdotr = 0;

  for(int fld = 0; fld < Nfields; fld++)
    for(hlong e = 0; e < Nelements; ++e)
      for(int i = 0; i < p_Np; ++i) {
        const hlong n = e * p_Np + i + fld * offset;

        dfloat rn = cpu_r[n];
        dfloat zn = cpu_z[n];
        dfloat sn = cpu_s[n];
        dfloat Sn = cpu_S[n];

        dfloat invDeg = (useWeight) ? cpu_invDegree[n]: 1.0;

        rn = rn - alpha * sn;
        zn = zn - alpha * Sn;

        rdotz += rn * zn * invDeg;
        zdotz += zn * zn * invDeg;
        rdotr += rn * rn * invDeg;

        cpu_r[n] = rn;
        cpu_z[n] = zn;
      }

  localdots[0] = rdotz;
  localdots[1] = zdotz;
  localdots[2] = rdotr;

#undef p_Np
}

void ellipticSerialUpdate2NBPCG(const int Nq,
                                const hlong Nelements,
                                const int useWeight,
                                occa::memory &o_invDegree,
                                occa::memory &o_s,
                                occa::memory &o_S,
                                const dfloat alpha,
                                occa::memory &o_r,
                                occa::memory &o_z,
                                dfloat* localdots)
{
  const dfloat* __restrict__ cpu_s  = (dfloat*)__builtin_assume_aligned(o_s.ptr(),
                                                                        USE_OCCA_MEM_BYTE_ALIGN);
  const dfloat* __restrict__ cpu_S  = (dfloat*)__builtin_assume_aligned(o_S.ptr(),
                                                                        USE_OCCA_MEM_BYTE_ALIGN);

  dfloat* __restrict__ cpu_r  =
    (dfloat*)__builtin_assume_aligned(o_r.ptr(), USE_OCCA_MEM_BYTE_ALIGN);
  dfloat* __restrict__ cpu_z  =
    (dfloat*)__builtin_assume_aligned(o_z.ptr(), USE_OCCA_MEM_BYTE_ALIGN);

  const dfloat* __restrict__ cpu_invDegree = (dfloat*)__builtin_assume_aligned(
    o_invDegree.ptr(),
    USE_OCCA_MEM_BYTE_ALIGN);

  switch(Nq) {
  case  2: ellipticSerialUpdate2NBPCGKernel <  2 >
    (Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, localdots);
    break;
  case  3: ellipticSerialUpdate2NBPCGKernel <  3 >
    (Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, localdots);
    break;
  case  4: ellipticSerialUpdate2NBPCGKernel <  4 >
    (Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, localdots);
    break;
  case  5: ellipticSerialUpdate2NBPCGKernel <  5 >
    (Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, localdots);
    break;
  case  6: ellipticSerialUpdate2NBPCGKernel <  6 >
    (Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, localdots);
    break;
  case  7: ellipticSerialUpdate2NBPCGKernel <  7 >
    (Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, localdots);
    break;
  case  8: ellipticSerialUpdate2NBPCGKernel <  8 >
    (Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, localdots);
    break;
  case  9: ellipticSerialUpdate2NBPCGKernel <  9 >
    (Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, localdots);
    break;
  case 10: ellipticSerialUpdate2NBPCGKernel < 10 >
    (Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, localdots);
    break;
  case 11: ellipticSerialUpdate2NBPCGKernel < 11 >
    (Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, localdots);
    break;
  case 12: ellipticSerialUpdate2NBPCGKernel < 12 >
    (Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z, localdots);
    break;
  }
}

void ellipticBlockSerialUpdate2NBPCG(const int Nfields,
                                     const hlong offset,
                                     const int Nq,
                                     const hlong Nelements,
                                     const int useWeight,
                                     occa::memory &o_invDegree,
                                     occa::memory &o_s,
                                     occa::memory &o_S,
                                     const dfloat alpha,
                                     occa::memory &o_r,
                                     occa::memory &o_z,
                                     dfloat* localdots)
{
  const dfloat* __restrict__ cpu_s  = (dfloat*)__builtin_assume_aligned(o_s.ptr(),
                                                                        USE_OCCA_MEM_BYTE_ALIGN);
  const dfloat* __restrict__ cpu_S  = (dfloat*)__builtin_assume_aligned(o_S.ptr(),
                                                                        USE_OCCA_MEM_BYTE_ALIGN);

  dfloat* __restrict__ cpu_r  =
    (dfloat*)__builtin_assume_aligned(o_r.ptr(), USE_OCCA_MEM_BYTE_ALIGN);
  dfloat* __restrict__ cpu_z  =
    (dfloat*)__builtin_assume_aligned(o_z.ptr(), USE_OCCA_MEM_BYTE_ALIGN);

  const dfloat* __restrict__ cpu_invDegree = (dfloat*)__builtin_assume_aligned(
    o_invDegree.ptr(),
    USE_OCCA_MEM_BYTE_ALIGN);

  switch(Nq) {
  case  2: ellipticBlockSerialUpdate2NBPCGKernel <  2 >
    (Nfields, offset, Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z,
     localdots);
    break;
  case  3: ellipticBlockSerialUpdate2NBPCGKernel <  3 >
    (Nfields, offset, Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z,
     localdots);
    break;
  case  4: ellipticBlockSerialUpdate2NBPCGKernel <  4 >
    (Nfields, offset, Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z,
     localdots);
    break;
  case  5: ellipticBlockSerialUpdate2NBPCGKernel <  5 >
    (Nfields, offset, Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z,
     localdots);
    break;
  case  6: ellipticBlockSerialUpdate2NBPCGKernel <  6 >
    (Nfields, offset, Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z,
     localdots);
    break;
  case  7: ellipticBlockSerialUpdate2NBPCGKernel <  7 >
    (Nfields, offset, Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z,
     localdots);
    break;
  case  8: ellipticBlockSerialUpdate2NBPCGKernel <  8 >
    (Nfields, offset, Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z,
     localdots);
    break;
  case  9: ellipticBlockSerialUpdate2NBPCGKernel <  9 >
    (Nfields, offset, Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z,
     localdots);
    break;
  case 10: ellipticBlockSerialUpdate2NBPCGKernel < 10 >
    (Nfields, offset, Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z,
     localdots);
    break;
  case 11: ellipticBlockSerialUpdate2NBPCGKernel < 11 >
    (Nfields, offset, Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z,
     localdots);
    break;
  case 12: ellipticBlockSerialUpdate2NBPCGKernel < 12 >
    (Nfields, offset, Nelements, useWeight,  cpu_invDegree,cpu_s, cpu_S, alpha, cpu_r, cpu_z,
     localdots);
    break;
  }
}

void ellipticNonBlockingUpdate2NBPCG(elliptic_t* elliptic,
                                     occa::memory &o_s, occa::memory &o_S, const dfloat alpha,
                                     occa::memory &o_r, occa::memory &o_z,
                                     dfloat* localdots, dfloat* globaldots, MPI_Request* request)
{
  setupAide &options = elliptic->options;

  int fixedIterationCountFlag = 0;
  int enableGatherScatters = 1;
  int enableReductions = 1;
  int verbose = verbose  = options.compareArgs("VERBOSE", "TRUE");
  int serial = options.compareArgs("THREAD MODEL", "SERIAL");
  int continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");

  options.getArgs("DEBUG ENABLE REDUCTIONS", enableReductions);
  options.getArgs("DEBUG ENABLE OGS", enableGatherScatters);
  if(options.compareArgs("FIXED ITERATION COUNT", "TRUE"))
    fixedIterationCountFlag = 1;


  mesh_t* mesh = elliptic->mesh;
  const dlong Nlocal = mesh->Np * mesh->Nelements;

  localdots[0] = 0; // rdotz
  localdots[1] = 0; // zdotz
  localdots[2] = 0; // rdotr

  int useWeight = (continuous != 0);

  if(serial == 1) {
    if(elliptic->blockSolver)
      ellipticBlockSerialUpdate2NBPCG(elliptic->Nfields, elliptic->Ntotal,
                                      mesh->Nq, mesh->Nelements, useWeight,
                                      elliptic->o_invDegree,
                                      o_s, o_S, alpha, o_r, o_z,
                                      localdots);

    else
      ellipticSerialUpdate2NBPCG(mesh->Nq, mesh->Nelements, useWeight,
                                 elliptic->o_invDegree,
                                 o_s, o_S, alpha, o_r, o_z,
                                 localdots);

  }else {
    // r <= r - alpha*s
    // z <= z - alpha*S
    // dot(r,z)
    // dot(z,z)
    // dot(r,r)
    if(elliptic->blockSolver)
      elliptic->update2NBPCGKernel(Nlocal,
                                   elliptic->Ntotal,
                                   elliptic->NblocksUpdatePCG,
                                   useWeight,
                                   elliptic->o_invDegree,
                                   o_s,
                                   o_S,
                                   alpha,
                                   o_r,
                                   o_z,
                                   elliptic->o_tmprdotz,
                                   elliptic->o_tmpzdotz,
                                   elliptic->o_tmprdotr);
    else
      elliptic->update2NBPCGKernel(mesh->Nelements * mesh->Np,
                                   elliptic->NblocksUpdatePCG,
                                   useWeight,
                                   elliptic->o_invDegree,
                                   o_s,
                                   o_S,
                                   alpha,
                                   o_r,
                                   o_z,
                                   elliptic->o_tmprdotz,
                                   elliptic->o_tmpzdotz,
                                   elliptic->o_tmprdotr);

    elliptic->o_tmprdotz.copyTo(elliptic->tmprdotz);
    elliptic->o_tmpzdotz.copyTo(elliptic->tmpzdotz);
    elliptic->o_tmprdotr.copyTo(elliptic->tmprdotr);

    for(int n = 0; n < elliptic->NblocksUpdatePCG; ++n) {
      localdots[0] += elliptic->tmprdotz[n];
      localdots[1] += elliptic->tmpzdotz[n];
      localdots[2] += elliptic->tmprdotr[n];
    }
  }

  globaldots[0] = 1;
  globaldots[1] = 1;
  globaldots[2] = 1;
  if(enableReductions)
    MPI_Iallreduce(localdots, globaldots, 3, MPI_DFLOAT, MPI_SUM, mesh->comm, request);
}

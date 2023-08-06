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

#include <stdio.h>
#include <stdlib.h>

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"

void meshPhysicalNodesHex3D(mesh_t *mesh)
{
  mesh->x = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs) * mesh->Np,sizeof(dfloat));
  mesh->y = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs) * mesh->Np,sizeof(dfloat));
  mesh->z = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs) * mesh->Np,sizeof(dfloat));

  const int nx1 = nekData.nx1;

  dfloat* xm1 = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat* ym1 = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat* zm1 = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
 
  dlong cnt = 0;
  for(dlong e = 0; e < mesh->Nelements; ++e) { /* for each element */
    hlong offset = e * nx1 * nx1 * nx1;
    nek::map_m_to_n(xm1, mesh->Nq, &nekData.xm1[offset], nx1);
    nek::map_m_to_n(ym1, mesh->Nq, &nekData.ym1[offset], nx1);
    nek::map_m_to_n(zm1, mesh->Nq, &nekData.zm1[offset], nx1);
 
    for(int n = 0; n < mesh->Np; ++n) { /* for each node */
      /* physical coordinate of interpolation node */
      mesh->x[cnt] = xm1[n];
      mesh->y[cnt] = ym1[n];
      mesh->z[cnt] = zm1[n];
      cnt++;
    }
  }
 
  free(xm1);
  free(ym1);
  free(zm1);
}

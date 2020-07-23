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

#include <stdlib.h>
#include <stdio.h>
#include "mesh3D.h"

void meshPrint3D(mesh3D* mesh)
{
  printf("EToV:\n");
  for(dlong e = 0; e < mesh->Nelements; ++e) {
    for(int v = 0; v < mesh->Nverts; ++v)
      printf(hlongFormat " ", mesh->EToV[e * mesh->Nverts + v]);
    printf("\n");
  }

  printf("EToE:\n");
  for(dlong e = 0; e < mesh->Nelements; ++e) {
    for(int f = 0; f < mesh->Nfaces; ++f)
      printf(dlongFormat " ",  mesh->EToE[e * mesh->Nfaces + f]);
    printf("\n");
  }

  printf("EToB:\n");
  for(dlong e = 0; e < mesh->Nelements; ++e) {
    for(int f = 0; f < mesh->Nfaces; ++f)
      printf("%d ",  mesh->EToB[e * mesh->Nfaces + f]);
    printf("\n");
  }

  printf("EToP:\n");
  for(dlong e = 0; e < mesh->Nelements; ++e) {
    for(int f = 0; f < mesh->Nfaces; ++f)
      printf("%d ",  mesh->EToP[e * mesh->Nfaces + f]);
    printf("\n");
  }
}

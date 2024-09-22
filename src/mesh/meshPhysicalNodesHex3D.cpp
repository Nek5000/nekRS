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

#include "platform.hpp"
#include "nekInterfaceAdapter.hpp"

void meshPhysicalNodesHex3D(mesh_t *mesh)
{
  std::vector<dfloat> xm(mesh->Nlocal);
  std::vector<dfloat> ym(mesh->Nlocal);
  std::vector<dfloat> zm(mesh->Nlocal);

  nek::xm1N(xm.data(), ym.data(), zm.data(), mesh->N, mesh->Nelements);

  mesh->o_x =
    platform->device.malloc<dfloat>(mesh->Nlocal);
  mesh->o_y =
    platform->device.malloc<dfloat>(mesh->Nlocal);
  mesh->o_z =
    platform->device.malloc<dfloat>(mesh->Nlocal);

  mesh->o_x.copyFrom(xm.data());
  mesh->o_y.copyFrom(ym.data());
  mesh->o_z.copyFrom(zm.data());
}

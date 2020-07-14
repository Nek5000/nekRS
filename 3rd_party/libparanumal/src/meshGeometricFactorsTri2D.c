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
#include "mesh2D.h"

void meshGeometricFactorsTri2D(mesh2D *mesh){

  /* unified storage array for geometric factors */
  mesh->Nvgeo = 5;
  mesh->vgeo = (dfloat*) calloc(mesh->Nelements*mesh->Nvgeo, 
				sizeof(dfloat));
  
  /* number of second order geometric factors */
  mesh->Nggeo = 4;
  mesh->ggeo = (dfloat*) calloc(mesh->Nelements*mesh->Nggeo, sizeof(dfloat));
  

  for(dlong e=0;e<mesh->Nelements;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    dlong id = e*mesh->Nverts+0;

    dfloat xe1 = mesh->EX[id+0];
    dfloat xe2 = mesh->EX[id+1];
    dfloat xe3 = mesh->EX[id+2];

    dfloat ye1 = mesh->EY[id+0];
    dfloat ye2 = mesh->EY[id+1];
    dfloat ye3 = mesh->EY[id+2];

    /* compute geometric factors for affine coordinate transform*/
    dfloat J = 0.25*((xe2-xe1)*(ye3-ye1) - (xe3-xe1)*(ye2-ye1));

    if(J<0) printf("bugger: got negative geofac\n");
    dfloat rx =  (0.5/J)*(ye3-ye1);
    dfloat ry = -(0.5/J)*(xe3-xe1);
    dfloat sx = -(0.5/J)*(ye2-ye1);
    dfloat sy =  (0.5/J)*(xe2-xe1);
    
    /* store geometric factors */
    mesh->vgeo[mesh->Nvgeo*e + RXID] = rx;
    mesh->vgeo[mesh->Nvgeo*e + RYID] = ry;
    mesh->vgeo[mesh->Nvgeo*e + SXID] = sx;
    mesh->vgeo[mesh->Nvgeo*e + SYID] = sy;
    mesh->vgeo[mesh->Nvgeo*e +  JID] = J;

    /* store second order geometric factors */
    mesh->ggeo[mesh->Nggeo*e + G00ID] = J*(rx*rx + ry*ry);
    mesh->ggeo[mesh->Nggeo*e + G01ID] = J*(rx*sx + ry*sy);
    mesh->ggeo[mesh->Nggeo*e + G11ID] = J*(sx*sx + sy*sy);
    mesh->ggeo[mesh->Nggeo*e + GWJID]  = J;
  }
}

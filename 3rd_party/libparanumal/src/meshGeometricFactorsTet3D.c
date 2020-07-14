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
#include "mesh3D.h"

void meshGeometricFactorsTet3D(mesh3D *mesh){

  /* unified storage array for geometric factors */
  mesh->Nvgeo = 12;
  mesh->vgeo = (dfloat*) calloc(mesh->Nelements*mesh->Nvgeo, 
        sizeof(dfloat));

  /* number of second order geometric factors */
  mesh->Nggeo = 7;
  mesh->ggeo = (dfloat*) calloc(mesh->Nelements*mesh->Nggeo, sizeof(dfloat));
  

  dfloat minJ = 1e9, maxJ = -1e9;
  for(dlong e=0;e<mesh->Nelements;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    dlong id = e*mesh->Nverts;

    /* vertex coordinates */
    dfloat xe1 = mesh->EX[id+0], ye1 = mesh->EY[id+0], ze1 = mesh->EZ[id+0];
    dfloat xe2 = mesh->EX[id+1], ye2 = mesh->EY[id+1], ze2 = mesh->EZ[id+1];
    dfloat xe3 = mesh->EX[id+2], ye3 = mesh->EY[id+2], ze3 = mesh->EZ[id+2];
    dfloat xe4 = mesh->EX[id+3], ye4 = mesh->EY[id+3], ze4 = mesh->EZ[id+3];

    /* Jacobian matrix */
    dfloat xr = 0.5*(xe2-xe1), xs = 0.5*(xe3-xe1), xt = 0.5*(xe4-xe1);
    dfloat yr = 0.5*(ye2-ye1), ys = 0.5*(ye3-ye1), yt = 0.5*(ye4-ye1);
    dfloat zr = 0.5*(ze2-ze1), zs = 0.5*(ze3-ze1), zt = 0.5*(ze4-ze1);

    /* compute geometric factors for affine coordinate transform*/
    dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);
    
    dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
    dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
    dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;

    if(J<0) printf("bugger: got negative geofac\n");
    minJ = mymin(minJ,J);
    maxJ = mymax(maxJ,J);
    
    /* store geometric factors */
    mesh->vgeo[mesh->Nvgeo*e + RXID] = rx;
    mesh->vgeo[mesh->Nvgeo*e + RYID] = ry;
    mesh->vgeo[mesh->Nvgeo*e + RZID] = rz;
    mesh->vgeo[mesh->Nvgeo*e + SXID] = sx;
    mesh->vgeo[mesh->Nvgeo*e + SYID] = sy;
    mesh->vgeo[mesh->Nvgeo*e + SZID] = sz;
    mesh->vgeo[mesh->Nvgeo*e + TXID] = tx;
    mesh->vgeo[mesh->Nvgeo*e + TYID] = ty;
    mesh->vgeo[mesh->Nvgeo*e + TZID] = tz;
    mesh->vgeo[mesh->Nvgeo*e +  JID] = J;
    //    printf("geo: %g,%g,%g - %g,%g,%g - %g,%g,%g\n",
    //     rx,ry,rz, sx,sy,sz, tx,ty,tz);

    /* store second order geometric factors */
    mesh->ggeo[mesh->Nggeo*e + G00ID] = J*(rx*rx + ry*ry + rz*rz);
    mesh->ggeo[mesh->Nggeo*e + G01ID] = J*(rx*sx + ry*sy + rz*sz);
    mesh->ggeo[mesh->Nggeo*e + G02ID] = J*(rx*tx + ry*ty + rz*tz);
    mesh->ggeo[mesh->Nggeo*e + G11ID] = J*(sx*sx + sy*sy + sz*sz);
    mesh->ggeo[mesh->Nggeo*e + G12ID] = J*(sx*tx + sy*ty + sz*tz);
    mesh->ggeo[mesh->Nggeo*e + G22ID] = J*(tx*tx + ty*ty + tz*tz);
    mesh->ggeo[mesh->Nggeo*e + GWJID] = J;
  }

  //printf("minJ = %g, maxJ = %g\n", minJ, maxJ);
}

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

void meshGeometricFactorsQuad2D(mesh2D *mesh){

  /* unified storage array for geometric factors */
  mesh->Nvgeo = 7;
  
  /* note that we have volume geometric factors for each node */
  mesh->vgeo = (dfloat*) calloc(mesh->Nelements*mesh->Nvgeo*mesh->Np, sizeof(dfloat));
  
  mesh->cubvgeo = (dfloat*) calloc(mesh->Nelements*mesh->Nvgeo*mesh->cubNp, sizeof(dfloat));

  /* number of second order geometric factors */
  mesh->Nggeo = 4;
  mesh->ggeo = (dfloat*) calloc(mesh->Nelements*mesh->Nggeo*mesh->Np, sizeof(dfloat));
  
  for(dlong e=0;e<mesh->Nelements;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    dlong id = e*mesh->Nverts;

    dfloat *xe = mesh->EX + id;
    dfloat *ye = mesh->EY + id;

    for(int n=0;n<mesh->Np;++n){

      /* local node coordinates */
      dfloat rn = mesh->r[n]; 
      dfloat sn = mesh->s[n];

      /* Jacobian matrix */
      dfloat xr = 0.25*( (1-sn)*(xe[1]-xe[0]) + (1+sn)*(xe[2]-xe[3]) );
      dfloat xs = 0.25*( (1-rn)*(xe[3]-xe[0]) + (1+rn)*(xe[2]-xe[1]) );
      dfloat yr = 0.25*( (1-sn)*(ye[1]-ye[0]) + (1+sn)*(ye[2]-ye[3]) );
      dfloat ys = 0.25*( (1-rn)*(ye[3]-ye[0]) + (1+rn)*(ye[2]-ye[1]) );
      
      /* compute geometric factors for affine coordinate transform*/
      dfloat J = xr*ys - xs*yr;

      if(J<1e-8) { printf("Negative or small Jacobian: %g\n", J); exit(-1);}
      dfloat rx =  ys/J;
      dfloat ry = -xs/J;
      dfloat sx = -yr/J;
      dfloat sy =  xr/J;

      int i = n%mesh->Nq;
      int j = n/mesh->Nq;
      dfloat JW = J*mesh->gllw[i]*mesh->gllw[j];
      
      /* store geometric factors */
      mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*RXID] = rx;
      mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*RYID] = ry;
      mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*SXID] = sx;
      mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*SYID] = sy;
      mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*JID]  = J;
      mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*JWID] = JW;
      mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*IJWID] = 1./JW;
      
      /* store second order geometric factors */
      mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*G00ID] = JW*(rx*rx + ry*ry);
      mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*G01ID] = JW*(rx*sx + ry*sy);
      mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*G11ID] = JW*(sx*sx + sy*sy);
      mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*GWJID] = JW;
    }

    //geometric data for quadrature
    for(int j=0;j<mesh->cubNq;++j){
      for(int i=0;i<mesh->cubNq;++i){
        
        dfloat rn = mesh->cubr[i];
        dfloat sn = mesh->cubr[j];

        /* Jacobian matrix */
        dfloat xr = 0.25*( (1-sn)*(xe[1]-xe[0]) + (1+sn)*(xe[2]-xe[3]) );
        dfloat xs = 0.25*( (1-rn)*(xe[3]-xe[0]) + (1+rn)*(xe[2]-xe[1]) );
        dfloat yr = 0.25*( (1-sn)*(ye[1]-ye[0]) + (1+sn)*(ye[2]-ye[3]) );
        dfloat ys = 0.25*( (1-rn)*(ye[3]-ye[0]) + (1+rn)*(ye[2]-ye[1]) );
        
        /* compute geometric factors for affine coordinate transform*/
        dfloat J = xr*ys - xs*yr;

        if(J<1e-8) { printf("Negative or small Jacobian: %g\n", J); exit(-1);}
        dfloat rx =  ys/J;
        dfloat ry = -xs/J;
        dfloat sx = -yr/J;
        dfloat sy =  xr/J;

        dfloat JW = J*mesh->cubw[i]*mesh->cubw[j];
        
        /* store geometric factors */
        dlong base = mesh->Nvgeo*mesh->cubNp*e + i + j*mesh->cubNq;
        mesh->cubvgeo[base + mesh->cubNp*RXID] = rx;
        mesh->cubvgeo[base + mesh->cubNp*RYID] = ry;
        mesh->cubvgeo[base + mesh->cubNp*SXID] = sx;
        mesh->cubvgeo[base + mesh->cubNp*SYID] = sy;
        mesh->cubvgeo[base + mesh->cubNp*JID]  = J;
        mesh->cubvgeo[base + mesh->cubNp*JWID] = JW;
        mesh->cubvgeo[base + mesh->cubNp*IJWID] = 1./JW;
      }
    }
  }
}

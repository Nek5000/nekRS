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

// custom geometric factors specialized for 3D quad on sphere

void meshGeometricFactorsQuad3D(mesh_t *mesh){

  /* unified storage array for geometric factors */
  mesh->Nvgeo = 12; // 
  
  /* note that we have volume geometric factors for each node */
  mesh->vgeo = (dfloat*) calloc(mesh->Nelements*mesh->Nvgeo*mesh->Np, sizeof(dfloat));

  mesh->cubvgeo = (dfloat*) calloc(mesh->Nelements*mesh->Nvgeo*mesh->cubNp, sizeof(dfloat));

  // Can be computed on the fly
  mesh->Nggeo = 7; 
  mesh->ggeo  = (dfloat *) calloc(mesh->Nelements*mesh->Np*mesh->Nggeo, sizeof(dfloat));

  dfloat *cxr = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));
  dfloat *cxs = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));
  dfloat *cyr = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));
  dfloat *cys = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));
  dfloat *czr = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));
  dfloat *czs = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));
  dfloat *cx  = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));
  dfloat *cy  = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));
  dfloat *cz  = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));
  
  for(int e=0;e<mesh->Nelements;++e){ /* for each element */
    
    for(int n=0;n<mesh->cubNq*mesh->cubNq;++n){
      cxr[n] = 0; cyr[n] = 0; czr[n] = 0;
      cxs[n] = 0; cys[n] = 0; czs[n] = 0;
      cx[n] = 0;  cy[n] = 0;  cz[n] = 0;
    }
    
    for(int j=0;j<mesh->Nq;++j){
      for(int i=0;i<mesh->Nq;++i){
  
  dfloat xij = mesh->x[i+j*mesh->Nq+e*mesh->Np];
  dfloat yij = mesh->y[i+j*mesh->Nq+e*mesh->Np];
  dfloat zij = mesh->z[i+j*mesh->Nq+e*mesh->Np];

  dfloat xr = 0, yr = 0, zr = 0;
  dfloat xs = 0, ys = 0, zs = 0;
  
  for(int n=0;n<mesh->Nq;++n){

    dfloat Din = mesh->D[i*mesh->Nq+n];
    dfloat Djn = mesh->D[j*mesh->Nq+n];

    xr += Din*mesh->x[n+j*mesh->Nq+e*mesh->Np];
    yr += Din*mesh->y[n+j*mesh->Nq+e*mesh->Np];
    zr += Din*mesh->z[n+j*mesh->Nq+e*mesh->Np];

    xs += Djn*mesh->x[i+n*mesh->Nq+e*mesh->Np];
    ys += Djn*mesh->y[i+n*mesh->Nq+e*mesh->Np];
    zs += Djn*mesh->z[i+n*mesh->Nq+e*mesh->Np];

  }

  {
    dfloat rx = ys*zij - zs*yij; // dXds x X
    dfloat ry = zs*xij - xs*zij;
    dfloat rz = xs*yij - ys*xij;
    
    dfloat sx = zr*yij - yr*zij; // -dXdr x X
    dfloat sy = xr*zij - zr*xij;
    dfloat sz = yr*xij - xr*yij;
    
    dfloat tx = yr*zs - zr*ys; // dXdr x dXds ~ X*|dXdr x dXds|/|X|
    dfloat ty = zr*xs - xr*zs;
    dfloat tz = xr*ys - yr*xs;
    
    dfloat Gx = tx, Gy = ty, Gz = tz;
    
    dfloat J = xij*tx + yij*ty + zij*tz;
    
    if(J<1e-8) { printf("Negative or small Jacobian: %g\n", J); exit(-1);}
    
    rx /= J;      sx /= J;      tx /= J;
    ry /= J;      sy /= J;      ty /= J;
    rz /= J;      sz /= J;      tz /= J;

    // use this for "volume" Jacobian
    dfloat Jnew = sqrt(Gx*Gx+Gy*Gy+Gz*Gz);  //(difference between actual Jacobian and sphere Jac)
    J = Jnew;
    
    if(J<1e-8) { printf("Negative or small Jacobian: %g\n", J); exit(-1);}
    //    printf("before: grad r = %g,%g,%g\n", rx, ry, rz);
  }

  dfloat GG00 = xr*xr+yr*yr+zr*zr;
  dfloat GG11 = xs*xs+ys*ys+zs*zs;
  dfloat GG01 = xr*xs+yr*ys+zr*zs;
  dfloat detGG = GG00*GG11 - GG01*GG01;

  // are these tangential
  dfloat rx = (xr*GG11-xs*GG01)/detGG;
  dfloat ry = (yr*GG11-ys*GG01)/detGG;
  dfloat rz = (zr*GG11-zs*GG01)/detGG;

  dfloat sx = (-xr*GG01+xs*GG00)/detGG;
  dfloat sy = (-yr*GG01+ys*GG00)/detGG;
  dfloat sz = (-zr*GG01+zs*GG00)/detGG;

  dfloat tx = yr*zs - zr*ys; // dXdr x dXds ~ X*|dXdr x dXds|/|X|
  dfloat ty = zr*xs - xr*zs;
  dfloat tz = xr*ys - yr*xs;

  // use this for "volume" Jacobian
  dfloat J = sqrt(tx*tx+ty*ty+tz*tz); // (difference between actual Jacobian and sphere Jac)

  //  printf("after: grad r = %g,%g,%g\n", rx, ry, rz);
  
  dfloat JW = J*mesh->gllw[i]*mesh->gllw[j];
  
  /* store geometric factors */
  int base = mesh->Nvgeo*mesh->Np*e + j*mesh->Nq + i;

  mesh->vgeo[base + mesh->Np*RXID] = rx;
  mesh->vgeo[base + mesh->Np*RYID] = ry;
  mesh->vgeo[base + mesh->Np*RZID] = rz;
  mesh->vgeo[base + mesh->Np*SXID] = sx;
  mesh->vgeo[base + mesh->Np*SYID] = sy;
  mesh->vgeo[base + mesh->Np*SZID] = sz;
  mesh->vgeo[base + mesh->Np*TXID] = tx;
  mesh->vgeo[base + mesh->Np*TYID] = ty;
  mesh->vgeo[base + mesh->Np*TZID] = tz;
  mesh->vgeo[base + mesh->Np*JID]  = J;
  mesh->vgeo[base + mesh->Np*JWID] = JW;
  mesh->vgeo[base + mesh->Np*IJWID] = 1./JW;

  /* store second order geometric factors (can be computed on the fly, later!!!)*/
  int gbase = mesh->Nggeo*mesh->Np*e + j*mesh->Nq + i;
  mesh->ggeo[gbase + mesh->Np*G00ID] = JW*(rx*rx + ry*ry + rz*rz);
  mesh->ggeo[gbase + mesh->Np*G01ID] = JW*(rx*sx + ry*sy + rz*sz);
  mesh->ggeo[gbase + mesh->Np*G02ID] = JW*(rx*tx + ry*ty + rz*tz); 

  mesh->ggeo[gbase + mesh->Np*G11ID] = JW*(sx*sx + sy*sy + sz*sz);
  mesh->ggeo[gbase + mesh->Np*G12ID] = JW*(sx*tx + sy*ty + sz*tz);

  mesh->ggeo[gbase + mesh->Np*G22ID] = JW*(tx*tx + ty*ty + tz*tz);
  mesh->ggeo[gbase + mesh->Np*GWJID] = JW;

  // now do for cubvgeo
  // 1. interpolate Jacobian matrix to cubature nodes
  for(int m=0;m<mesh->cubNq;++m){
    for(int n=0;n<mesh->cubNq;++n){
      dfloat cIni = mesh->cubInterp[n*mesh->Nq+i];
      dfloat cImj = mesh->cubInterp[m*mesh->Nq+j];
      cxr[n+m*mesh->cubNq] += cIni*cImj*xr;
      cxs[n+m*mesh->cubNq] += cIni*cImj*xs;
      cyr[n+m*mesh->cubNq] += cIni*cImj*yr;
      cys[n+m*mesh->cubNq] += cIni*cImj*ys;
      czr[n+m*mesh->cubNq] += cIni*cImj*zr;
      czs[n+m*mesh->cubNq] += cIni*cImj*zs;
      cx[n+m*mesh->cubNq] += cIni*cImj*xij;
      cy[n+m*mesh->cubNq] += cIni*cImj*yij;
      cz[n+m*mesh->cubNq] += cIni*cImj*zij;
    }
  }
      }
    }

    
    for(int n=0;n<mesh->cubNq*mesh->cubNq;++n){
      
      dfloat rx = cys[n]*cz[n] - czs[n]*cy[n]; // dXds x X
      dfloat ry = czs[n]*cx[n] - cxs[n]*cz[n];
      dfloat rz = cxs[n]*cy[n] - cys[n]*cx[n];
      
      dfloat sx = czr[n]*cy[n] - cyr[n]*cz[n]; // -dXdr x X
      dfloat sy = cxr[n]*cz[n] - czr[n]*cx[n];
      dfloat sz = cyr[n]*cx[n] - cxr[n]*cy[n];
      
      dfloat tx = cyr[n]*czs[n] - czr[n]*cys[n]; // dXdr x dXds ~ X*|dXdr x dXds|/|X|
      dfloat ty = czr[n]*cxs[n] - cxr[n]*czs[n];
      dfloat tz = cxr[n]*cys[n] - cyr[n]*cxs[n];
      
      dfloat Gx = tx, Gy = ty, Gz = tz;
      
      dfloat J = cx[n]*tx + cy[n]*ty + cz[n]*tz;
      
      if(J<1e-8) { printf("Negative or small Jacobian: %g\n", J); exit(-1);}
      
      rx /= J;      sx /= J;      tx /= J;
      ry /= J;      sy /= J;      ty /= J;
      rz /= J;      sz /= J;      tz /= J;
      
      // use this for "volume" Jacobian
      J = sqrt(Gx*Gx+Gy*Gy+Gz*Gz);
      
      if(J<1e-8) { printf("Negative or small cubature Jacobian: %g (Gx,y,z=%g,%g,%g)\n",
        J, Gx, Gy, Gz); exit(-1);}
      
      dfloat JW = J*mesh->cubw[n%mesh->cubNq]*mesh->cubw[n/mesh->cubNq];
      
      /* store geometric factors */
      int base = mesh->Nvgeo*mesh->cubNp*e + n;
      
      mesh->cubvgeo[base + mesh->cubNp*RXID] = rx;
      mesh->cubvgeo[base + mesh->cubNp*RYID] = ry;
      mesh->cubvgeo[base + mesh->cubNp*RZID] = rz;
      mesh->cubvgeo[base + mesh->cubNp*SXID] = sx;
      mesh->cubvgeo[base + mesh->cubNp*SYID] = sy;
      mesh->cubvgeo[base + mesh->cubNp*SZID] = sz;
      mesh->cubvgeo[base + mesh->cubNp*TXID] = tx;
      mesh->cubvgeo[base + mesh->cubNp*TYID] = ty;
      mesh->cubvgeo[base + mesh->cubNp*TZID] = tz;
      mesh->cubvgeo[base + mesh->cubNp*JID]  = J;
      mesh->cubvgeo[base + mesh->cubNp*JWID] = JW;
      mesh->cubvgeo[base + mesh->cubNp*IJWID] = 1./JW;
    } 
  }
}

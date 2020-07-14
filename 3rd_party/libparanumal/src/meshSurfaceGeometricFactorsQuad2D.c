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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

/* compute outwards facing normals, surface Jacobian, and volume Jacobian for all face nodes */
void meshSurfaceGeometricFactorsQuad2D(mesh2D *mesh){

  /* unified storage array for geometric factors */
  mesh->Nsgeo = 7;
  mesh->sgeo = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*
                                mesh->Nsgeo*mesh->Nfp*mesh->Nfaces, 
                                sizeof(dfloat));

  mesh->cubsgeo = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*
                                mesh->Nsgeo*mesh->cubNq*mesh->Nfaces, 
                                sizeof(dfloat));
  
  for(dlong e=0;e<mesh->Nelements+mesh->totalHaloPairs;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    dlong id = e*mesh->Nverts;

    dfloat *xe = mesh->EX + id;
    dfloat *ye = mesh->EY + id;
    
    for(int f=0;f<mesh->Nfaces;++f){ // for each face
      
      for(int i=0;i<mesh->Nfp;++i){  // for each node on face

        /* volume index of face node */
        int n = mesh->faceNodes[f*mesh->Nfp+i];

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
        
        /* face f normal and length */
        dfloat nx =   ye[(f+1)%mesh->Nverts]-ye[f];
        dfloat ny = -(xe[(f+1)%mesh->Nverts]-xe[f]);
        dfloat  d = norm2(nx,ny);

        /* output index */
        dlong base = mesh->Nsgeo*(mesh->Nfaces*mesh->Nfp*e + mesh->Nfp*f + i);

        /* store normal, surface Jacobian, and reciprocal of volume Jacobian */
        mesh->sgeo[base+NXID] = nx/d;
        mesh->sgeo[base+NYID] = ny/d;
        mesh->sgeo[base+SJID] = d/2.;
        mesh->sgeo[base+IJID] = 1./J;

        mesh->sgeo[base+WIJID] = 1./(J*mesh->gllw[0]);
        mesh->sgeo[base+WSJID] = (d/2.)*mesh->gllw[i];
      }

      //geometric data for quadrature
      for(int i=0;i<mesh->cubNq;++i){  // for each quadrature node on face

        dfloat rn = 0., sn = 0.;

        /* interpolate local node coordinates */
        for (int j=0;j<mesh->Nfp;j++) {
          /* volume index of face node */
          int n = mesh->faceNodes[f*mesh->Nfp+j];

          rn += mesh->cubInterp[i*mesh->Nfp+j]*mesh->r[n];
          sn += mesh->cubInterp[i*mesh->Nfp+j]*mesh->s[n];
        }
        
        /* Jacobian matrix */
        dfloat xr = 0.25*( (1-sn)*(xe[1]-xe[0]) + (1+sn)*(xe[2]-xe[3]) );
        dfloat xs = 0.25*( (1-rn)*(xe[3]-xe[0]) + (1+rn)*(xe[2]-xe[1]) );
        dfloat yr = 0.25*( (1-sn)*(ye[1]-ye[0]) + (1+sn)*(ye[2]-ye[3]) );
        dfloat ys = 0.25*( (1-rn)*(ye[3]-ye[0]) + (1+rn)*(ye[2]-ye[1]) );
        
        /* compute geometric factors for affine coordinate transform*/
        dfloat J = xr*ys - xs*yr;
        
        /* face f normal and length */
        dfloat nx =   ye[(f+1)%mesh->Nverts]-ye[f];
        dfloat ny = -(xe[(f+1)%mesh->Nverts]-xe[f]);
        dfloat  d = norm2(nx,ny);

        /* output index */
        dlong base = mesh->Nsgeo*(mesh->Nfaces*mesh->cubNq*e + mesh->cubNq*f + i);

        /* store normal, surface Jacobian, and reciprocal of volume Jacobian */
        mesh->cubsgeo[base+NXID] = nx/d;
        mesh->cubsgeo[base+NYID] = ny/d;
        mesh->cubsgeo[base+SJID] = d/2.;
        mesh->cubsgeo[base+IJID] = 1./J;

        mesh->cubsgeo[base+WIJID] = 1./(J*mesh->cubw[0]);
        mesh->cubsgeo[base+WSJID] = (d/2.)*mesh->cubw[i];
      }
    }
  }

  for(dlong e=0;e<mesh->Nelements;++e){ /* for each non-halo element */
    for(int n=0;n<mesh->Nfp*mesh->Nfaces;++n){
      dlong baseM = e*mesh->Nfp*mesh->Nfaces + n;
      dlong baseP = mesh->mapP[baseM];
      if(baseP<0) baseP = baseM;
      
      // rescaling - missing factor of 2 ? (only impacts penalty and thus stiffness)
      dfloat hinvM = mesh->sgeo[baseM*mesh->Nsgeo + SJID]*mesh->sgeo[baseM*mesh->Nsgeo + IJID];
      dfloat hinvP = mesh->sgeo[baseP*mesh->Nsgeo + SJID]*mesh->sgeo[baseP*mesh->Nsgeo + IJID];
      
      mesh->sgeo[baseM*mesh->Nsgeo+IHID] = mymax(hinvM,hinvP);
      mesh->sgeo[baseP*mesh->Nsgeo+IHID] = mymax(hinvM,hinvP);
    }
  }  
}

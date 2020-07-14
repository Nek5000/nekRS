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
#include "mesh3D.h"

void computeFrame(dfloat nx, dfloat ny, dfloat nz,
		  dfloat &tanx, dfloat &tany, dfloat &tanz,
		  dfloat &binx, dfloat &biny, dfloat &binz){

  dfloat rdotn, ranx, rany, ranz;
  do{
    ranx = drand48();
    rany = drand48();
    ranz = drand48();
    
    dfloat magran = sqrt(ranx*ranx+rany*rany+ranz*ranz);
    
    ranx /= magran;
    rany /= magran;
    ranz /= magran;
    
    rdotn = nx*ranx+ny*rany+nz*ranz;
  }while(fabs(rdotn)<1e-4);
  
  tanx = ny*ranz - nz*rany;
  tany = nz*ranx - nx*ranz;
  tanz = nx*rany - ny*ranx;

  dfloat magtan = sqrt(tanx*tanx+tany*tany+tanz*tanz);

  tanx /= magtan;
  tany /= magtan;
  tanz /= magtan;

  binx = ny*tanz - nz*tany;
  biny = nz*tanx - nx*tanz;
  binz = nx*tany - ny*tanx;

  dfloat magbin = sqrt(binx*binx+biny*biny+binz*binz);

  binx /= magbin;
  biny /= magbin;
  binz /= magbin;

  //  printf("nor = %g,%g,%g; tan = %g,%g,%g; bin = %g,%g,%g\n", nx, ny, nz, tanx, tany, tanz, binx, biny, binz);
}

void interpolateFaceHex3D(int *faceNodes, dfloat *I, dfloat *x, int N, dfloat *Ix, int M){
  
  dfloat *Ix0 = (dfloat*) calloc(N*N, sizeof(dfloat));
  dfloat *Ix1 = (dfloat*) calloc(N*M, sizeof(dfloat));
  
  for(int j=0;j<N;++j){
    for(int i=0;i<N;++i){
      Ix0[j*N+i] = x[faceNodes[j*N+i]];
    }
  }
  
  for(int j=0;j<N;++j){
    for(int i=0;i<M;++i){
      dfloat tmp = 0;
      for(int n=0;n<N;++n){
	tmp += I[i*N + n]*Ix0[j*N+n];
      }
      Ix1[j*M+i] = tmp;
    }
  }

  for(int j=0;j<M;++j){
    for(int i=0;i<M;++i){
      dfloat tmp = 0;
      for(int n=0;n<N;++n){
	tmp += I[j*N + n]*Ix1[n*M+i];
      }
      Ix[j*M+i] = tmp;
    }
  }

  free(Ix0);
  free(Ix1);
  
}


/* compute outwards facing normals, surface Jacobian, and volume Jacobian for all face nodes */
void meshSurfaceGeometricFactorsHex3D(mesh3D *mesh){

  /* unified storage array for geometric factors */
  mesh->Nsgeo = 17;
  mesh->sgeo = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*
                                mesh->Nsgeo*mesh->Nfp*mesh->Nfaces, 
                                sizeof(dfloat));

  mesh->cubsgeo = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*
				   mesh->Nsgeo*mesh->cubNfp*mesh->Nfaces, 
				   sizeof(dfloat));
  
  dfloat *xre = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *xse = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *xte = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *yre = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *yse = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *yte = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *zre = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *zse = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *zte = (dfloat*) calloc(mesh->Np, sizeof(dfloat));

  dfloat *cubxre = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));
  dfloat *cubxse = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));
  dfloat *cubxte = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));
  dfloat *cubyre = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));
  dfloat *cubyse = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));
  dfloat *cubyte = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));
  dfloat *cubzre = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));
  dfloat *cubzse = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));
  dfloat *cubzte = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));

  dfloat *cubxe = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));
  dfloat *cubye = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));
  dfloat *cubze = (dfloat*) calloc(mesh->cubNq*mesh->cubNq, sizeof(dfloat));
  
  for(dlong e=0;e<mesh->Nelements+mesh->totalHaloPairs;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    dlong id = e*mesh->Nverts;

    dfloat *xe = mesh->EX + id;
    dfloat *ye = mesh->EY + id;
    dfloat *ze = mesh->EZ + id;

    for(int n=0;n<mesh->Np;++n){
      xre[n] = 0; xse[n] = 0; xte[n] = 0;
      yre[n] = 0; yse[n] = 0; yte[n] = 0; 
      zre[n] = 0; zse[n] = 0; zte[n] = 0; 
    }

    for(int k=0;k<mesh->Nq;++k){
      for(int j=0;j<mesh->Nq;++j){
        for(int i=0;i<mesh->Nq;++i){
	  
          int n = i + j*mesh->Nq + k*mesh->Nq*mesh->Nq;

          /* local node coordinates */
          dfloat rn = mesh->r[n]; 
          dfloat sn = mesh->s[n];
          dfloat tn = mesh->t[n];

	  for(int m=0;m<mesh->Nq;++m){
	    int idr = e*mesh->Np + k*mesh->Nq*mesh->Nq + j*mesh->Nq + m;
	    int ids = e*mesh->Np + k*mesh->Nq*mesh->Nq + m*mesh->Nq + i;
	    int idt = e*mesh->Np + m*mesh->Nq*mesh->Nq + j*mesh->Nq + i;
	    xre[n] += mesh->D[i*mesh->Nq+m]*mesh->x[idr];
	    xse[n] += mesh->D[j*mesh->Nq+m]*mesh->x[ids];
	    xte[n] += mesh->D[k*mesh->Nq+m]*mesh->x[idt];
	    yre[n] += mesh->D[i*mesh->Nq+m]*mesh->y[idr];
	    yse[n] += mesh->D[j*mesh->Nq+m]*mesh->y[ids];
	    yte[n] += mesh->D[k*mesh->Nq+m]*mesh->y[idt];
	    zre[n] += mesh->D[i*mesh->Nq+m]*mesh->z[idr];
	    zse[n] += mesh->D[j*mesh->Nq+m]*mesh->z[ids];
	    zte[n] += mesh->D[k*mesh->Nq+m]*mesh->z[idt];
	  }
	}
      }
    }
	  
    for(int f=0;f<mesh->Nfaces;++f){ // for each face
      
      for(int i=0;i<mesh->Nfp;++i){  // for each node on face

        /* volume index of face node */
        int n = mesh->faceNodes[f*mesh->Nfp+i];

        /* local node coordinates */
        dfloat rn = mesh->r[n]; 
        dfloat sn = mesh->s[n];
        dfloat tn = mesh->t[n];
	
	dfloat xr = xre[n], xs = xse[n], xt = xte[n];
	dfloat yr = yre[n], ys = yse[n], yt = yte[n];
	dfloat zr = zre[n], zs = zse[n], zt = zte[n];
	
        /* determinant of Jacobian matrix */
        dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);
        
        dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
        dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
        dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;
        
        /* face f normal and length */
        dfloat nx, ny, nz, d;
        switch(f){
        case 0: nx = -tx; ny = -ty; nz = -tz; break;
        case 1: nx = -sx; ny = -sy; nz = -sz; break;
        case 2: nx = +rx; ny = +ry; nz = +rz; break;
        case 3: nx = +sx; ny = +sy; nz = +sz; break;
        case 4: nx = -rx; ny = -ry; nz = -rz; break;
        case 5: nx = +tx; ny = +ty; nz = +tz; break;
        }

        dfloat sJ = sqrt(nx*nx+ny*ny+nz*nz);
        nx /= sJ; ny /= sJ; nz /= sJ;
        sJ *= J;
        
        /* output index */
        dlong base = mesh->Nsgeo*(mesh->Nfaces*mesh->Nfp*e + mesh->Nfp*f + i);

        /* store normal, surface Jacobian, and reciprocal of volume Jacobian */
        mesh->sgeo[base+NXID] = nx;
        mesh->sgeo[base+NYID] = ny;
        mesh->sgeo[base+NZID] = nz;
        mesh->sgeo[base+SJID] = sJ;
        mesh->sgeo[base+IJID] = 1./J;

        mesh->sgeo[base+WIJID] = 1./(J*mesh->gllw[0]);
        mesh->sgeo[base+WSJID] = sJ*mesh->gllw[i%mesh->Nq]*mesh->gllw[i/mesh->Nq];

	computeFrame(nx, ny, nz,
		     mesh->sgeo[base+STXID], mesh->sgeo[base+STYID], mesh->sgeo[base+STZID],
		     mesh->sgeo[base+SBXID], mesh->sgeo[base+SBYID], mesh->sgeo[base+SBZID]);
      }
    
      // now interpolate geofacs to cubature
      interpolateFaceHex3D(mesh->faceNodes+f*mesh->Nfp, mesh->cubInterp, xre, mesh->Nq, cubxre, mesh->cubNq);
      interpolateFaceHex3D(mesh->faceNodes+f*mesh->Nfp, mesh->cubInterp, xse, mesh->Nq, cubxse, mesh->cubNq);
      interpolateFaceHex3D(mesh->faceNodes+f*mesh->Nfp, mesh->cubInterp, xte, mesh->Nq, cubxte, mesh->cubNq);
      interpolateFaceHex3D(mesh->faceNodes+f*mesh->Nfp, mesh->cubInterp, yre, mesh->Nq, cubyre, mesh->cubNq);
      interpolateFaceHex3D(mesh->faceNodes+f*mesh->Nfp, mesh->cubInterp, yse, mesh->Nq, cubyse, mesh->cubNq);
      interpolateFaceHex3D(mesh->faceNodes+f*mesh->Nfp, mesh->cubInterp, yte, mesh->Nq, cubyte, mesh->cubNq);
      interpolateFaceHex3D(mesh->faceNodes+f*mesh->Nfp, mesh->cubInterp, zre, mesh->Nq, cubzre, mesh->cubNq);
      interpolateFaceHex3D(mesh->faceNodes+f*mesh->Nfp, mesh->cubInterp, zse, mesh->Nq, cubzse, mesh->cubNq);
      interpolateFaceHex3D(mesh->faceNodes+f*mesh->Nfp, mesh->cubInterp, zte, mesh->Nq, cubzte, mesh->cubNq);

      interpolateFaceHex3D(mesh->faceNodes+f*mesh->Nfp, mesh->cubInterp, mesh->x+e*mesh->Np, mesh->Nq, cubxe, mesh->cubNq);
      interpolateFaceHex3D(mesh->faceNodes+f*mesh->Nfp, mesh->cubInterp, mesh->y+e*mesh->Np, mesh->Nq, cubye, mesh->cubNq);
      interpolateFaceHex3D(mesh->faceNodes+f*mesh->Nfp, mesh->cubInterp, mesh->z+e*mesh->Np, mesh->Nq, cubze, mesh->cubNq);
      
      //geometric data for quadrature
      for(int i=0;i<mesh->cubNfp;++i){  // for each quadrature node on face

	dfloat xr = cubxre[i], xs = cubxse[i], xt = cubxte[i];
	dfloat yr = cubyre[i], ys = cubyse[i], yt = cubyte[i];
	dfloat zr = cubzre[i], zs = cubzse[i], zt = cubzte[i];

        /* determinant of Jacobian matrix */
        dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);
        
        dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
        dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
        dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;
        
        /* face f normal and length */
        dfloat nx, ny, nz, d;
        switch(f){
        case 0: nx = -tx; ny = -ty; nz = -tz; break;
        case 1: nx = -sx; ny = -sy; nz = -sz; break;
        case 2: nx = +rx; ny = +ry; nz = +rz; break;
        case 3: nx = +sx; ny = +sy; nz = +sz; break;
        case 4: nx = -rx; ny = -ry; nz = -rz; break;
        case 5: nx = +tx; ny = +ty; nz = +tz; break;
        }

        dfloat sJ = sqrt(nx*nx+ny*ny+nz*nz);
        nx /= sJ; ny /= sJ; nz /= sJ;
        sJ *= J;
        

        /* output index */
        dlong base = mesh->Nsgeo*(mesh->Nfaces*mesh->cubNfp*e + mesh->cubNfp*f + i);

        /* store normal, surface Jacobian, and reciprocal of volume Jacobian */
        mesh->cubsgeo[base+NXID] = nx;
        mesh->cubsgeo[base+NYID] = ny;
        mesh->cubsgeo[base+NZID] = nz;
        mesh->cubsgeo[base+SJID] = sJ;
        mesh->cubsgeo[base+IJID] = 1./J;

        mesh->cubsgeo[base+WIJID] = 1./(J*mesh->cubw[0]);
        mesh->cubsgeo[base+WSJID] = sJ*mesh->cubw[i%mesh->cubNq]*mesh->cubw[i/mesh->cubNq];

        mesh->cubsgeo[base+SURXID] = cubxe[i];
	mesh->cubsgeo[base+SURYID] = cubye[i];
	mesh->cubsgeo[base+SURZID] = cubze[i];
	
	computeFrame(nx, ny, nz,
		     mesh->cubsgeo[base+STXID], mesh->cubsgeo[base+STYID], mesh->cubsgeo[base+STZID],
		     mesh->cubsgeo[base+SBXID], mesh->cubsgeo[base+SBYID], mesh->cubsgeo[base+SBZID]);
	
	
      }
    }
  }

  for(dlong e=0;e<mesh->Nelements;++e){ /* for each non-halo element */
    for(int n=0;n<mesh->Nfp*mesh->Nfaces;++n){
      dlong baseM = e*mesh->Nfp*mesh->Nfaces + n;
      dlong baseP = mesh->mapP[baseM];
      // rescaling - missing factor of 2 ? (only impacts penalty and thus stiffness)
      dfloat hinvM = mesh->sgeo[baseM*mesh->Nsgeo + SJID]*mesh->sgeo[baseM*mesh->Nsgeo + IJID];
      dfloat hinvP = mesh->sgeo[baseP*mesh->Nsgeo + SJID]*mesh->sgeo[baseP*mesh->Nsgeo + IJID];
      mesh->sgeo[baseM*mesh->Nsgeo+IHID] = mymax(hinvM,hinvP);
      mesh->sgeo[baseP*mesh->Nsgeo+IHID] = mymax(hinvM,hinvP);
    }
  }


  free(xre); free(xse); free(xte);
  free(yre); free(yse); free(yte);
  free(zre); free(zse); free(zte);

  free(cubxre); free(cubxse); free(cubxte);
  free(cubyre); free(cubyse); free(cubyte);
  free(cubzre); free(cubzse); free(cubzte);
  
}

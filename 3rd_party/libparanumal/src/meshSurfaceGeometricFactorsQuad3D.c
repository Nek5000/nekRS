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

/* compute outwards facing normals, surface Jacobian, and volume Jacobian for all face nodes */
void meshSurfaceGeometricFactorsQuad3D(mesh_t *mesh){

  /* unified storage array for geometric factors */
  mesh->Nsgeo = 14; // fix later
  mesh->sgeo = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*
				mesh->Nsgeo*mesh->Nfp*mesh->Nfaces, 
				sizeof(dfloat));

  mesh->cubsgeo = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*
                                mesh->Nsgeo*mesh->cubNq*mesh->Nfaces, 
                                sizeof(dfloat));
  
  dfloat *cubx = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*
				  mesh->cubNq*mesh->Nfaces, sizeof(dfloat));

  dfloat *cuby = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*
				  mesh->cubNq*mesh->Nfaces, sizeof(dfloat));

  dfloat *cubz = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*
				  mesh->cubNq*mesh->Nfaces, sizeof(dfloat));
  


  dfloat *xr = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *yr = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *zr = (dfloat*) calloc(mesh->Np, sizeof(dfloat));

  dfloat *xs = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *ys = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *zs = (dfloat*) calloc(mesh->Np, sizeof(dfloat));

  dfloat *J  = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  
  for(int e=0;e<mesh->Nelements+mesh->totalHaloPairs;++e){ /* for each element */

    for(int j=0;j<mesh->Nq;++j){
      for(int i=0;i<mesh->Nq;++i){

	dfloat x = mesh->x[i+j*mesh->Nq+e*mesh->Np];
	dfloat y = mesh->y[i+j*mesh->Nq+e*mesh->Np];
	dfloat z = mesh->z[i+j*mesh->Nq+e*mesh->Np];
	
	dfloat xrij = 0, yrij = 0, zrij = 0;
	dfloat xsij = 0, ysij = 0, zsij = 0;

	for(int n=0;n<mesh->Nq;++n){

	  dfloat Din = mesh->D[i*mesh->Nq+n];
	  dfloat Djn = mesh->D[j*mesh->Nq+n];

	  xrij += Din*mesh->x[n+j*mesh->Nq+e*mesh->Np];
	  yrij += Din*mesh->y[n+j*mesh->Nq+e*mesh->Np];
	  zrij += Din*mesh->z[n+j*mesh->Nq+e*mesh->Np];

	  xsij += Djn*mesh->x[i+n*mesh->Nq+e*mesh->Np];
	  ysij += Djn*mesh->y[i+n*mesh->Nq+e*mesh->Np];
	  zsij += Djn*mesh->z[i+n*mesh->Nq+e*mesh->Np];

	}

	dfloat txij = yrij*zsij - zrij*ysij;
	dfloat tyij = zrij*xsij - xrij*zsij;
	dfloat tzij = xrij*ysij - yrij*xsij;

	dfloat Gx = txij, Gy = tyij, Gz = tzij;
	
	dfloat Jij = x*txij + y*tyij + z*tzij;
	
	xr[i+j*mesh->Nq] = xrij;
	yr[i+j*mesh->Nq] = yrij;
	zr[i+j*mesh->Nq] = zrij;
			       
	xs[i+j*mesh->Nq] = xsij;
	ys[i+j*mesh->Nq] = ysij;
	zs[i+j*mesh->Nq] = zsij;

	J[i+j*mesh->Nq] = sqrt(Gx*Gx+Gy*Gy+Gz*Gz);
      }
    }

    for(int f=0;f<mesh->Nfaces;++f){
      for(int n=0;n<mesh->Nq;++n){
	int id = mesh->faceNodes[n+f*mesh->Nq];

	dfloat xid = mesh->x[id+e*mesh->Np];
	dfloat yid = mesh->y[id+e*mesh->Np];
	dfloat zid = mesh->z[id+e*mesh->Np];
	dfloat Jid = J[id];
	
	dfloat nx, ny, nz;

	if(f==0){
	  nx = yr[id]*zid - zr[id]*yid;
	  ny = zr[id]*xid - xr[id]*zid;
	  nz = xr[id]*yid - yr[id]*xid;
	}

	if(f==1){
	  nx = ys[id]*zid - zs[id]*yid;
	  ny = zs[id]*xid - xs[id]*zid;
	  nz = xs[id]*yid - ys[id]*xid;
	}

	if(f==2){
	  nx = -yr[id]*zid + zr[id]*yid;
	  ny = -zr[id]*xid + xr[id]*zid;
	  nz = -xr[id]*yid + yr[id]*xid;
	}

	if(f==3){
	  nx = -ys[id]*zid + zs[id]*yid;
	  ny = -zs[id]*xid + xs[id]*zid;
	  nz = -xs[id]*yid + ys[id]*xid;
	}
	
	dfloat R = sqrt(xid*xid+yid*yid+zid*zid);
	
	nx /= R;
	ny /= R;
	nz /= R;

	dfloat sJ = sqrt(nx*nx+ny*ny+nz*nz);

	nx /= sJ;
	ny /= sJ;
	nz /= sJ;

	if(sJ<1e-8) { printf("Negative or small surface Jacobian: %g\n", sJ); exit(-1);}
	
	int base = mesh->Nsgeo*(e*mesh->Nq*mesh->Nfaces + n + f*mesh->Nq);
	
	mesh->sgeo[base+NXID] = nx;
	mesh->sgeo[base+NYID] = ny;
	mesh->sgeo[base+NZID] = nz;
	mesh->sgeo[base+SJID] = sJ;

	mesh->sgeo[base+IJID] = 1./Jid;

	mesh->sgeo[base+WIJID] = 1./(Jid*mesh->gllw[0]);
	mesh->sgeo[base+WSJID] = sJ*mesh->gllw[n];
      }
    }

    // interpolate geofacs to surface quadrature
    for(int f=0;f<mesh->Nfaces;++f){

      for(int n=0;n<mesh->cubNq;++n){
	dfloat cxr = 0, cxs = 0, cx = 0;
	dfloat cyr = 0, cys = 0, cy = 0;
	dfloat czr = 0, czs = 0, cz = 0;
	
	for(int i=0;i<mesh->Nq;++i){
	  int id = mesh->faceNodes[i+f*mesh->Nq];
	  dfloat cIni = mesh->cubInterp[n*mesh->Nq+i];
	  cxr += cIni*xr[id];
	  cxs += cIni*xs[id];
	  cyr += cIni*yr[id];
	  cys += cIni*ys[id];
	  czr += cIni*zr[id];
	  czs += cIni*zs[id];
	  cx  += cIni*mesh->x[id+e*mesh->Np];
	  cy  += cIni*mesh->y[id+e*mesh->Np];
	  cz  += cIni*mesh->z[id+e*mesh->Np];
	}

	cubx[e*mesh->cubNq*mesh->Nfaces+f*mesh->cubNq + n] = cx;
	cuby[e*mesh->cubNq*mesh->Nfaces+f*mesh->cubNq + n] = cy;
	cubz[e*mesh->cubNq*mesh->Nfaces+f*mesh->cubNq + n] = cz;
	
	dfloat Gx = cyr*czs - czr*cys;
	dfloat Gy = czr*cxs - cxr*czs;
	dfloat Gz = cxr*cys - cyr*cxs;	
	dfloat cJ = sqrt(Gx*Gx+Gy*Gy+Gz*Gz);
	dfloat volJ = cx*Gx + cy*Gy + cz*Gz; // xij*tx + yij*ty + zij*tz;	
	dfloat nx, ny, nz;

	if(f==0){
	  nx = cyr*cz - czr*cy;
	  ny = czr*cx - cxr*cz;
	  nz = cxr*cy - cyr*cx;
	}
	
	if(f==1){
	  nx = cys*cz - czs*cy;
	  ny = czs*cx - cxs*cz;
	  nz = cxs*cy - cys*cx;
	}

	if(f==2){
	  nx = -cyr*cz + czr*cy;
	  ny = -czr*cx + cxr*cz;
	  nz = -cxr*cy + cyr*cx;
	}

	if(f==3){
	  nx = -cys*cz + czs*cy;
	  ny = -czs*cx + cxs*cz;
	  nz = -cxs*cy + cys*cx;
	}

	dfloat R = sqrt(cx*cx+cy*cy+cz*cz);
	
	nx /= R;
	ny /= R;
	nz /= R;
	
	dfloat sJ = sqrt(nx*nx+ny*ny+nz*nz);

	nx /= sJ;
	ny /= sJ;
	nz /= sJ;

	if(sJ<1e-8) { printf("Negative or small surface Jacobian: %g\n", sJ); exit(-1);}
	
	int base = mesh->Nsgeo*(e*mesh->cubNq*mesh->Nfaces + n + f*mesh->cubNq);
	
	mesh->cubsgeo[base+NXID] = nx;
	mesh->cubsgeo[base+NYID] = ny;
	mesh->cubsgeo[base+NZID] = nz;
	mesh->cubsgeo[base+SJID] = sJ;
	mesh->cubsgeo[base+IHID] = sJ/volJ;
	//	mesh->cubsgeo[base+WSJID] = sJ*mesh->cubw[n];
      }
    }
  }


#if 0
  for(int e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      for(int n=0;n<mesh->Nq;++n){
	int idM = n+f*mesh->Nq+e*mesh->Nfaces*mesh->Nq;
	int idP = mesh->mapP[idM];
	int eP = idP/(mesh->Nq*mesh->Nfaces);
	int fP = (idP%(mesh->Nq*mesh->Nfaces))/mesh->Nq;
	int nP = (idP%mesh->Nq);
	int baseM = e*mesh->Nq*mesh->Nfaces*mesh->Nsgeo + f*mesh->Nq*mesh->Nsgeo + n;
	int baseP = eP*mesh->Nq*mesh->Nfaces*mesh->Nsgeo + fP*mesh->Nq*mesh->Nsgeo + nP;
	printf("e,f,n=(%d,%d,%d)-(%d,%d,%d): xP-xM=(%g,%g,%g) : norP+norM=%g,%g,%g\n",
	       e,f,n,eP,fP,nP,
	       mesh->x[mesh->vmapP[idM]]-mesh->x[mesh->vmapM[idM]],
	       mesh->y[mesh->vmapP[idM]]-mesh->y[mesh->vmapM[idM]],
	       mesh->z[mesh->vmapP[idM]]-mesh->z[mesh->vmapM[idM]],
	       mesh->sgeo[baseM+NXID*mesh->Nq]+mesh->sgeo[baseP+NXID*mesh->Nq],
	       mesh->sgeo[baseM+NYID*mesh->Nq]+mesh->sgeo[baseP+NYID*mesh->Nq],
	       mesh->sgeo[baseM+NZID*mesh->Nq]+mesh->sgeo[baseP+NZID*mesh->Nq]);

      }
    }
  }
#endif  
  // TW: omit 1/min(h) calculation

  for(dlong e=0;e<mesh->Nelements;++e){ /* for each non-halo element */
    for(int n=0;n<mesh->Nfp*mesh->Nfaces;++n){
      dlong baseM = e*mesh->Nfp*mesh->Nfaces + n;
      dlong baseP = mesh->mapP[baseM];
      if(baseP<0) baseP = baseM;
      
      // rescaling - missing factor of 2 ? (only impacts penalty and thus stiffness)
      dfloat hinvM = mesh->sgeo[baseM*mesh->Nsgeo + SJID]*mesh->sgeo[baseM*mesh->Nsgeo + IJID];
      dfloat hinvP = mesh->sgeo[baseP*mesh->Nsgeo + SJID]*mesh->sgeo[baseP*mesh->Nsgeo + IJID];

      //      printf("hinvM/P = %g,%g\n", hinvM, hinvP);
      
      mesh->sgeo[baseM*mesh->Nsgeo+IHID] = mymax(hinvM,hinvP);
      mesh->sgeo[baseP*mesh->Nsgeo+IHID] = mymax(hinvM,hinvP);
    }
  }  

  for(dlong e=0;e<mesh->Nelements;++e){ /* for each non-halo element */
    for(int f=0;f<mesh->Nfaces;++f){
      dlong eP = mesh->EToE[e*mesh->Nfaces+f];
      dlong fP = mesh->EToF[e*mesh->Nfaces+f];

      dfloat maxhinv  = 0;
      for(int n=0;n<mesh->cubNq;++n){
	dlong idM = e*mesh->cubNq*mesh->Nfaces+f*mesh->cubNq+n;
	dfloat cxM = cubx[idM];
	dfloat cyM = cuby[idM];
	dfloat czM = cubz[idM];
	
	dfloat mindist2;
	int minidP = 0;
	// jump through hoops to find neighbor cubature node
	// [ not needed elsewhere since we interpolate consistently ]
	dlong idP;
	for(int m=0;m<mesh->cubNq;++m){
	  idP = eP*mesh->cubNq*mesh->Nfaces+fP*mesh->cubNq+m;

	  dfloat cxP = cubx[idP];
	  dfloat cyP = cuby[idP];
	  dfloat czP = cubz[idP];

	  dfloat dist2 = pow(cxP-cxM,2)+pow(cyP-cyM,2)+pow(czP-czM,2);

	  if(m==0 || dist2<mindist2){
	    mindist2 = dist2;
	    minidP = m;
	  }
	}

	if(mindist2>1e-12)
	printf("mindist2 = %g\n", mindist2);
	
	idM = mesh->Nsgeo*( e*mesh->cubNq*mesh->Nfaces+ f*mesh->cubNq+n)+IHID;
	idP = mesh->Nsgeo*(eP*mesh->cubNq*mesh->Nfaces+fP*mesh->cubNq+minidP)+IHID;
	
	dfloat hinv = mymax(mesh->cubsgeo[idM],mesh->cubsgeo[idP]);
	mesh->cubsgeo[idM] = hinv;
	mesh->cubsgeo[idP] = hinv;

      }
    }  
  } 

  
}

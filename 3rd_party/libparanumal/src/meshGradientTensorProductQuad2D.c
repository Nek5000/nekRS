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

// baseline tensor product mesh Gradient for quadrilateral elements
void meshGradientTensorProductQuad2D(mesh2D *mesh,
				     dfloat * q, 
				     dfloat * dqdx, 
				     dfloat * dqdy){
  
  // loop over elements
  for(int e=0;e<mesh->Nelements;++e){
    
    // compute gradient at each node
    for(int j=0;j<mesh->N+1;++j){
      for(int i=0;i<mesh->N+1;++i){
	
	// local node index 
	int n = i + (mesh->N+1)*j;
	
	// load geometric factors
	int gid = mesh->Np*mesh->Nvgeo*e + n;
	float drdx = vgeo[gid + mesh->Np*RXID];
	float drdy = vgeo[gid + mesh->Np*RYID];
	float dsdx = vgeo[gid + mesh->Np*SXID];
	float dsdy = vgeo[gid + mesh->Np*SYID];
	
	// matrix-vector multiplies
	dfloat dqdr = 0, dqds = 0;
	for(int m=0;m<mesh->N+1;++m){
	  dqdr += mesh->D[i*(mesh->N+1) + m]*q[m + j*(mesh->N+1) + e*mesh->Np];
	  dqds += mesh->D[j*(mesh->N+1) + m]*q[i + m*(mesh->N+1) + e*mesh->Np];
	}
	
	// chain rule
	dqdx[n+e*mesh->Np] = drdx*dqdr + dsdx*dqds;
	dqdy[n+e*mesh->Np] = drdy*dqdr + dsdy*dqds;
      }
    }
  }
}


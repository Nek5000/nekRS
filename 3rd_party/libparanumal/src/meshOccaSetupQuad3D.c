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
#include <string.h>
#include  "mpi.h"

#include "mesh3D.h"

void meshOccaSetupQuad3D(mesh_t *mesh, setupAide &newOptions, occa::properties &kernelInfo){

  // conigure device
  occaDeviceConfig(mesh, newOptions);

  //make seperate stream for halo exchange
  mesh->defaultStream = mesh->device.getStream();
  mesh->dataStream = mesh->device.createStream();
  mesh->device.setStream(mesh->defaultStream);

  // find elements that have all neighbors on this process
  dlong *internalElementIds = (dlong*) calloc(mesh->Nelements, sizeof(dlong));
  dlong *notInternalElementIds = (dlong*) calloc(mesh->Nelements, sizeof(dlong));

  dlong Ninterior = 0, NnotInterior = 0;
  for(dlong e=0;e<mesh->Nelements;++e){
    int flag = 0;
    for(int f=0;f<mesh->Nfaces;++f)
      if(mesh->EToP[e*mesh->Nfaces+f]!=-1)
        flag = 1;
    if(!flag)
      internalElementIds[Ninterior++] = e;
    else
      notInternalElementIds[NnotInterior++] = e;
  }

  //printf("NinteriorElements = %d, NnotInternalElements = %d\n", Ninterior, NnotInterior);

  mesh->NinternalElements = Ninterior;
  mesh->NnotInternalElements = NnotInterior;
  if(Ninterior)
    mesh->o_internalElementIds    = mesh->device.malloc(Ninterior*sizeof(dlong), internalElementIds);

  if(NnotInterior>0)
    mesh->o_notInternalElementIds = mesh->device.malloc(NnotInterior*sizeof(dlong), notInternalElementIds);

  // OCCA allocate device memory (remember to go back for halo)
  mesh->o_q =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->q);
  mesh->o_rhsq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);
  mesh->o_resq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->resq);
  
  size_t bytes = mesh->device.memoryAllocated();
  printf("bytes allocated: %lg\n", bytes/1.e9);
  
  //lumped mass matrix
  mesh->MM = (dfloat *) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  for (int j=0;j<mesh->Nq;j++) {
    for (int i=0;i<mesh->Nq;i++) {
      int n = i+j*mesh->Nq;
      mesh->MM[n+n*mesh->Np] = mesh->gllw[i]*mesh->gllw[j];
    }
  }
  
  dfloat *cubDWT = (dfloat*) calloc(mesh->cubNq*mesh->Nq, sizeof(dfloat));
  dfloat *cubProjectT = (dfloat*) calloc(mesh->cubNq*mesh->Nq, sizeof(dfloat));
  dfloat *cubInterpT = (dfloat*) calloc(mesh->cubNq*mesh->Nq, sizeof(dfloat));
  for(int n=0;n<mesh->Nq;++n){
    for(int m=0;m<mesh->cubNq;++m){
      cubDWT[n+m*mesh->Nq] = mesh->cubDW[n*mesh->cubNq+m];
      cubProjectT[n+m*mesh->Nq] = mesh->cubProject[n*mesh->cubNq+m];
      cubInterpT[m+n*mesh->cubNq] = mesh->cubInterp[m*mesh->Nq+n];
    }
  }

  mesh->intx = (dfloat*) calloc(mesh->Nelements*mesh->Nfaces*mesh->cubNq, sizeof(dfloat));
  mesh->inty = (dfloat*) calloc(mesh->Nelements*mesh->Nfaces*mesh->cubNq, sizeof(dfloat));
  mesh->intz = (dfloat*) calloc(mesh->Nelements*mesh->Nfaces*mesh->cubNq, sizeof(dfloat));
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      for(int n=0;n<mesh->cubNq;++n){
	dfloat ix = 0, iy = 0, iz = 0;
	for(int m=0;m<mesh->Nq;++m){
	  dlong vid = mesh->vmapM[m+f*mesh->Nfp+e*mesh->Nfp*mesh->Nfaces];
	  dfloat xm = mesh->x[vid];
	  dfloat ym = mesh->y[vid];
	  dfloat zm = mesh->z[vid];

	  dfloat Inm = mesh->cubInterp[m+n*mesh->Nq];
	  ix += Inm*xm;
	  iy += Inm*ym;
	  iz += Inm*ym;
	}
	dlong id = n + f*mesh->cubNq + e*mesh->Nfaces*mesh->cubNq;
	mesh->intx[id] = ix;
	mesh->inty[id] = iy;
	mesh->intz[id] = iz;
      }
    }
  }
  
  mesh->o_D = mesh->device.malloc(mesh->Nq*mesh->Nq*sizeof(dfloat), mesh->D);
  
  bytes = mesh->device.memoryAllocated();
  printf("bytes allocated: %lg\n", bytes/1.e9);
  
  // bundle D and (W^{-1} D^t W)
  dfloat *Dmatrices = (dfloat*) calloc(mesh->Nq*mesh->Nq*2, sizeof(dfloat));
  for(int n=0;n<mesh->Nq*mesh->Nq;++n){
    Dmatrices[n] = mesh->D[n];
  }
  for(int j=0;j<mesh->Nq;++j){
    for(int i=0;i<mesh->Nq;++i){
      // note minus
      Dmatrices[mesh->Nq*mesh->Nq + i + j*mesh->Nq] = -mesh->D[i*mesh->Nq + j]*mesh->gllw[i]/mesh->gllw[j];
    }
  }
  
  mesh->o_Dmatrices = mesh->device.malloc(2*mesh->Nq*mesh->Nq*sizeof(dfloat), Dmatrices);
  
  free(Dmatrices);
  
  mesh->o_Smatrices = mesh->device.malloc(mesh->Nq*mesh->Nq*sizeof(dfloat), mesh->D); //dummy
  
  mesh->o_vgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nvgeo*mesh->Np*sizeof(dfloat),
			mesh->vgeo);
  mesh->o_sgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->Nfp*mesh->Nsgeo*sizeof(dfloat),
			mesh->sgeo);
  mesh->o_ggeo =
    mesh->device.malloc(mesh->Nelements*mesh->Np*mesh->Nggeo*sizeof(dfloat),
			mesh->ggeo);
  
  bytes = mesh->device.memoryAllocated();
  printf("bytes allocated: %lg\n", bytes/1.e9);
  
  mesh->o_cubvgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nvgeo*mesh->cubNp*sizeof(dfloat),
			mesh->cubvgeo);
  
  mesh->o_cubsgeo =
      mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->cubNq*mesh->Nsgeo*sizeof(dfloat),
			  mesh->cubsgeo);
  
  mesh->o_cubInterpT =
    mesh->device.malloc(mesh->Nq*mesh->cubNq*sizeof(dfloat),
			cubInterpT);
  
  mesh->o_cubProjectT =
    mesh->device.malloc(mesh->Nq*mesh->cubNq*sizeof(dfloat),
			cubProjectT);
  
  mesh->o_cubDWT =
    mesh->device.malloc(mesh->Nq*mesh->cubNq*sizeof(dfloat),
			cubDWT);
  
  mesh->o_cubDWmatrices = mesh->device.malloc(mesh->cubNq*mesh->Nq*sizeof(dfloat), cubDWT);
  
  dfloat *LIFTT = (dfloat*) calloc(mesh->Np*mesh->Nfaces*mesh->Nfp, sizeof(dfloat));
  
  mesh->o_LIFTT =
    mesh->device.malloc(1*sizeof(dfloat)); // dummy
  
  mesh->o_intx =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->cubNq*sizeof(dfloat),
			mesh->intx);
  
  mesh->o_inty =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->cubNq*sizeof(dfloat),
			mesh->inty);
  
  mesh->o_intz =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->cubNq*sizeof(dfloat),
			mesh->intz);
  
  //dummy quadrature lifter operators
  mesh->o_intInterpT = mesh->device.malloc(mesh->cubNq*mesh->Nq*sizeof(dfloat));
  mesh->o_intInterpT.copyFrom(mesh->o_cubInterpT);
  
  mesh->o_intLIFTT = mesh->device.malloc(mesh->cubNq*mesh->Nq*sizeof(dfloat));
  mesh->o_intLIFTT.copyFrom(mesh->o_cubProjectT);
 
  
  mesh->o_MM =
    mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
			mesh->MM);
  
  mesh->o_vmapM =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(dlong),
			mesh->vmapM);
  
  mesh->o_vmapP =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(dlong),
			mesh->vmapP);
  
  mesh->o_EToB =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(int),
			mesh->EToB);
  
  mesh->o_x =
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat),
			mesh->x);
  
  mesh->o_y =
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat),
			mesh->y);
  
  // dummy z variables (note used y)
  mesh->o_z =
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat),
			mesh->z);
  
  if(mesh->totalHaloPairs>0){
    // copy halo element list to DEVICE
    mesh->o_haloElementList =
      mesh->device.malloc(mesh->totalHaloPairs*sizeof(dlong), mesh->haloElementList);
    
    // temporary DEVICE buffer for halo (maximum size Nfields*Np for dfloat)
    mesh->o_haloBuffer =
      mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat));
  }
  
  
  //-------------------------------------
  // NBN: 2 streams for async MPI updates
  // {Vol, Surf, update}  run on q[0]
  // {halo-get, copy} run on q[1]
  //-------------------------------------
  mesh->stream0 = mesh->device.getStream();
#ifdef USE_2_STREAMS
  mesh->stream1 = mesh->device.createStream();  // NBN: second stream
#else
  mesh->stream1 = mesh->stream0;                // NBN: stream1 == stream0
#endif
  mesh->device.setStream(mesh->stream0);
  //-------------------------------------
  
  kernelInfo["defines/" "p_Nfields"]= mesh->Nfields;
  kernelInfo["defines/" "p_N"]= mesh->N;
  kernelInfo["defines/" "p_Nq"]= mesh->N+1;
  kernelInfo["defines/" "p_Np"]= mesh->Np;
  kernelInfo["defines/" "p_Nfp"]= mesh->Nfp;
  kernelInfo["defines/" "p_Nfaces"]= mesh->Nfaces;
  kernelInfo["defines/" "p_NfacesNfp"]= mesh->Nfp*mesh->Nfaces;
  kernelInfo["defines/" "p_Nvgeo"]= mesh->Nvgeo;
  kernelInfo["defines/" "p_Nsgeo"]= mesh->Nsgeo;
  kernelInfo["defines/" "p_Nggeo"]= mesh->Nggeo;

  kernelInfo["defines/" "p_NXID"]= NXID;
  kernelInfo["defines/" "p_NYID"]= NYID;
  kernelInfo["defines/" "p_NZID"]= NZID;
  kernelInfo["defines/" "p_SJID"]= SJID;
  kernelInfo["defines/" "p_IJID"]= IJID;
  kernelInfo["defines/" "p_IHID"]= IHID;
  kernelInfo["defines/" "p_WIJID"]= WIJID;
  kernelInfo["defines/" "p_WSJID"]= WSJID;

  kernelInfo["defines/" "p_max_EL_nnz"]= mesh->max_EL_nnz; // for Bernstein Bezier lift

  kernelInfo["defines/" "p_cubNq"]= mesh->cubNq;
  kernelInfo["defines/" "p_cubNp"]= mesh->cubNp;
  kernelInfo["defines/" "p_intNfp"]= mesh->intNfp;
  kernelInfo["defines/" "p_intNfpNfaces"]= mesh->intNfp*mesh->Nfaces;

  if(sizeof(dfloat)==4){
    kernelInfo["defines/" "dfloat"]="float";
    kernelInfo["defines/" "dfloat2"]="float2";
    kernelInfo["defines/" "dfloat4"]="float4";
    kernelInfo["defines/" "dfloat8"]="float8";
  }
  if(sizeof(dfloat)==8){
    kernelInfo["defines/" "dfloat"]="double";
    kernelInfo["defines/" "dfloat2"]="double2";
    kernelInfo["defines/" "dfloat4"]="double4";
    kernelInfo["defines/" "dfloat8"]="double8";
  }

  if(sizeof(dlong)==4){
    kernelInfo["defines/" "dlong"]="int";
  }
  if(sizeof(dlong)==8){
    kernelInfo["defines/" "dlong"]="long long int";
  }

  if(mesh->device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
    kernelInfo["compiler_flags"] += " --ftz=true ";
    kernelInfo["compiler_flags"] += " --prec-div=false ";
    kernelInfo["compiler_flags"] += " --prec-sqrt=false ";
    kernelInfo["compiler_flags"] += " --use_fast_math ";
    kernelInfo["compiler_flags"] += " --fmad=true "; // compiler option for cuda
  }

  kernelInfo["defines/" "p_G00ID"]= G00ID;
  kernelInfo["defines/" "p_G01ID"]= G01ID;
  kernelInfo["defines/" "p_G02ID"]= G02ID;
  kernelInfo["defines/" "p_G11ID"]= G11ID;
  kernelInfo["defines/" "p_G12ID"]= G12ID;
  kernelInfo["defines/" "p_G22ID"]= G22ID;
  kernelInfo["defines/" "p_GWJID"]= GWJID;


  kernelInfo["defines/" "p_RXID"]= RXID;
  kernelInfo["defines/" "p_SXID"]= SXID;
  kernelInfo["defines/" "p_TXID"]= TXID;

  kernelInfo["defines/" "p_RYID"]= RYID;
  kernelInfo["defines/" "p_SYID"]= SYID;
  kernelInfo["defines/" "p_TYID"]= TYID;

  kernelInfo["defines/" "p_RZID"]= RZID;
  kernelInfo["defines/" "p_SZID"]= SZID;
  kernelInfo["defines/" "p_TZID"]= TZID;

  kernelInfo["defines/" "p_JID"]= JID;
  kernelInfo["defines/" "p_JWID"]= JWID;
  kernelInfo["defines/" "p_IJWID"]= IJWID;

}

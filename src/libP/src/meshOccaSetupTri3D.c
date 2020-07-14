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

void meshOccaSetupTri3D(mesh_t *mesh, setupAide &newOptions, occa::properties &kernelInfo){

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


  // build Dr, Ds, LIFT transposes
  dfloat *DrT = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *DsT = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Np;++m){
      DrT[n+m*mesh->Np] = mesh->Dr[n*mesh->Np+m];
      DsT[n+m*mesh->Np] = mesh->Ds[n*mesh->Np+m];
    }
  }
  
  // build Dr, Ds transposes
  dfloat *DrsT = (dfloat*) calloc(2*mesh->Np*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Np;++m){
      DrsT[n+m*mesh->Np] = mesh->Dr[n*mesh->Np+m];
      DrsT[n+m*mesh->Np+mesh->Np*mesh->Np] = mesh->Ds[n*mesh->Np+m];
    }
  }
  
  dfloat *LIFTT = (dfloat*) calloc(mesh->Np*mesh->Nfaces*mesh->Nfp, sizeof(dfloat));
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Nfaces*mesh->Nfp;++m){
      LIFTT[n+m*mesh->Np] = mesh->LIFT[n*mesh->Nfp*mesh->Nfaces+m];
    }
  }
  
  // build volume cubature matrix transposes
  int cubNpBlocked = mesh->Np*((mesh->cubNp+mesh->Np-1)/mesh->Np);
  dfloat *cubDrWT = (dfloat*) calloc(cubNpBlocked*mesh->Np, sizeof(dfloat));
  dfloat *cubDsWT = (dfloat*) calloc(cubNpBlocked*mesh->Np, sizeof(dfloat));
  dfloat *cubDrsWT = (dfloat*) calloc(2*mesh->cubNp*mesh->Np, sizeof(dfloat));
  dfloat *cubProjectT = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  dfloat *cubInterpT = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->cubNp;++m){
      cubDrWT[n+m*mesh->Np] = mesh->cubDrW[n*mesh->cubNp+m];
      cubDsWT[n+m*mesh->Np] = mesh->cubDsW[n*mesh->cubNp+m];
      
      cubDrsWT[n+m*mesh->Np] = mesh->cubDrW[n*mesh->cubNp+m];
      cubDrsWT[n+m*mesh->Np+mesh->cubNp*mesh->Np] = mesh->cubDsW[n*mesh->cubNp+m];
      
      cubProjectT[n+m*mesh->Np] = mesh->cubProject[n*mesh->cubNp+m];
      cubInterpT[m+n*mesh->cubNp] = mesh->cubInterp[m*mesh->Np+n];
    }
  }
  
  // build surface integration matrix transposes
  dfloat *intLIFTT = (dfloat*) calloc(mesh->Np*mesh->Nfaces*mesh->intNfp, sizeof(dfloat));
  dfloat *intInterpT = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces*mesh->intNfp, sizeof(dfloat));
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Nfaces*mesh->intNfp;++m){
      intLIFTT[n+m*mesh->Np] = mesh->intLIFT[n*mesh->intNfp*mesh->Nfaces+m];
    }
  }
  for(int n=0;n<mesh->intNfp*mesh->Nfaces;++n){
    for(int m=0;m<mesh->Nfp;++m){
      intInterpT[n+m*mesh->Nfaces*mesh->intNfp] = mesh->intInterp[n*mesh->Nfp + m];
    }
  }

  // =============== BB operators [added by JC] ===============
  // deriv operators: transpose from row major to column major
  int *D1ids = (int*) calloc(mesh->Np*3,sizeof(int));
  int *D2ids = (int*) calloc(mesh->Np*3,sizeof(int));
  int *D3ids = (int*) calloc(mesh->Np*3,sizeof(int));
  dfloat *Dvals = (dfloat*) calloc(mesh->Np*3,sizeof(dfloat));  
  
  dfloat *VBq = (dfloat*) calloc(mesh->Np*mesh->cubNp,sizeof(dfloat));
  dfloat *PBq = (dfloat*) calloc(mesh->Np*mesh->cubNp,sizeof(dfloat));
  
  dfloat *L0vals = (dfloat*) calloc(mesh->Nfp*3,sizeof(dfloat)); // tridiag
  int *ELids = (int*) calloc(1+mesh->Np*mesh->max_EL_nnz,sizeof(int));
  dfloat *ELvals = (dfloat*) calloc(1+mesh->Np*mesh->max_EL_nnz,sizeof(dfloat));
  
  for (int i = 0; i < mesh->Np; ++i){
    for (int j = 0; j < 3; ++j){
      D1ids[i+j*mesh->Np] = mesh->D1ids[j+i*3];
      D2ids[i+j*mesh->Np] = mesh->D2ids[j+i*3];
      D3ids[i+j*mesh->Np] = mesh->D3ids[j+i*3];      
      Dvals[i+j*mesh->Np] = mesh->Dvals[j+i*3];    
    }
  }
  
  for (int i = 0; i < mesh->cubNp; ++i){
    for (int j = 0; j < mesh->Np; ++j){
      VBq[i+j*mesh->cubNp] = mesh->VBq[j+i*mesh->Np];
      PBq[j+i*mesh->Np] = mesh->PBq[i+j*mesh->cubNp];
    }
  }

  
  for (int i = 0; i < mesh->Nfp; ++i){
    for (int j = 0; j < 3; ++j){
      L0vals[i+j*mesh->Nfp] = mesh->L0vals[j+i*3];
    }
  }
  
  for (int i = 0; i < mesh->Np; ++i){
    for (int j = 0; j < mesh->max_EL_nnz; ++j){
      ELids[i + j*mesh->Np] = mesh->ELids[j+i*mesh->max_EL_nnz];
      ELvals[i + j*mesh->Np] = mesh->ELvals[j+i*mesh->max_EL_nnz]; // ???
    }
  }
  
  //BB mass matrix
  mesh->BBMM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  for (int n = 0; n < mesh->Np; ++n){
    for (int m = 0; m < mesh->Np; ++m){
      for (int i = 0; i < mesh->Np; ++i){
	for (int j = 0; j < mesh->Np; ++j){
	  mesh->BBMM[n+m*mesh->Np] += mesh->VB[m+j*mesh->Np]*mesh->MM[i+j*mesh->Np]*mesh->VB[n+i*mesh->Np];
	}
      } 
    }
  }
  
  // =============== end BB stuff =============================
  

  //build element stiffness matrices
  dfloat *SrrT, *SrsT, *SsrT, *SssT;
  if (mesh->Nverts == 3) {
    mesh->Srr = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Srs = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Ssr = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Sss = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    for (int n=0;n<mesh->Np;n++) {
      for (int m=0;m<mesh->Np;m++) {
	for (int k=0;k<mesh->Np;k++) {
	  for (int l=0;l<mesh->Np;l++) {
	    mesh->Srr[m+n*mesh->Np] += mesh->Dr[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Dr[m+k*mesh->Np];
	    mesh->Srs[m+n*mesh->Np] += mesh->Dr[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Ds[m+k*mesh->Np];
	    mesh->Ssr[m+n*mesh->Np] += mesh->Ds[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Dr[m+k*mesh->Np];
	    mesh->Sss[m+n*mesh->Np] += mesh->Ds[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Ds[m+k*mesh->Np];
	  }
	} 
      }
    }
    SrrT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    SrsT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    SsrT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    SssT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    for (int n=0;n<mesh->Np;n++) {
      for (int m=0;m<mesh->Np;m++) {  
	SrrT[m+n*mesh->Np] = mesh->Srr[n+m*mesh->Np];
	SrsT[m+n*mesh->Np] = mesh->Srs[n+m*mesh->Np];
	SsrT[m+n*mesh->Np] = mesh->Ssr[n+m*mesh->Np];
	SssT[m+n*mesh->Np] = mesh->Sss[n+m*mesh->Np];     
      }
    }
  }

  dfloat *ST = (dfloat*) calloc(3*mesh->Np*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Np;++m){
      ST[n+m*mesh->Np+0*mesh->Np*mesh->Np] = mesh->Srr[n*mesh->Np+m];
      ST[n+m*mesh->Np+1*mesh->Np*mesh->Np] = mesh->Srs[n*mesh->Np+m]+mesh->Ssr[n*mesh->Np+m];
      ST[n+m*mesh->Np+2*mesh->Np*mesh->Np] = mesh->Sss[n*mesh->Np+m];
    }
  }
    
  mesh->intx = (dfloat*) calloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp, sizeof(dfloat));
  mesh->inty = (dfloat*) calloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp, sizeof(dfloat));
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      for(int n=0;n<mesh->intNfp;++n){
	dfloat ix = 0, iy = 0;
	for(int m=0;m<mesh->Nfp;++m){
	  dlong vid = mesh->vmapM[m+f*mesh->Nfp+e*mesh->Nfp*mesh->Nfaces];
	  dfloat xm = mesh->x[vid];
	  dfloat ym = mesh->y[vid];
	  //dfloat Inm = mesh->intInterp[n+f*mesh->intNfp+m*mesh->intNfp*mesh->Nfaces];
	  dfloat Inm = mesh->intInterp[m+n*mesh->Nfp+f*mesh->intNfp*mesh->Nfp]; // Fixed
	  ix += Inm*xm;
	  iy += Inm*ym;
	}
	dlong id = n + f*mesh->intNfp + e*mesh->Nfaces*mesh->intNfp;
	mesh->intx[id] = ix;
	mesh->inty[id] = iy;
      }
    }
  }

  mesh->o_Dr = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
				   mesh->Dr);

  mesh->o_Ds = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
				   mesh->Ds);

  mesh->o_DrT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
				    DrT);

  mesh->o_DsT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
				    DsT);

  mesh->o_DtT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
				    DsT); // note: dummy allocated with DsT

  mesh->o_Dmatrices = mesh->device.malloc(2*mesh->Np*mesh->Np*sizeof(dfloat), DrsT);

  mesh->o_LIFT =
    mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat),
			mesh->LIFT);

  mesh->o_LIFTT =
    mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat),
			LIFTT);

  mesh->o_vgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nvgeo*sizeof(dfloat),
			mesh->vgeo);

  mesh->o_cubvgeo =   mesh->device.malloc(sizeof(dfloat));// dummy
    
  mesh->o_sgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->Nsgeo*sizeof(dfloat),
			mesh->sgeo);

  mesh->o_ggeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nggeo*sizeof(dfloat),
			mesh->ggeo);

  mesh->o_cubsgeo = mesh->o_sgeo; //dummy cubature geo factors

  mesh->o_SrrT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), SrrT);
  mesh->o_SrsT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), SrsT);
  mesh->o_SsrT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), SsrT);
  mesh->o_SssT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), SssT);
  mesh->o_Smatrices = mesh->device.malloc(3*mesh->Np*mesh->Np*sizeof(dfloat), ST);

  mesh->o_D1ids = mesh->device.malloc(mesh->Np*3*sizeof(int),D1ids);
  mesh->o_D2ids = mesh->device.malloc(mesh->Np*3*sizeof(int),D2ids);
  mesh->o_D3ids = mesh->device.malloc(mesh->Np*3*sizeof(int),D3ids);
  mesh->o_Dvals = mesh->device.malloc(mesh->Np*3*sizeof(dfloat),Dvals);

  mesh->o_BBMM = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),mesh->BBMM);

  mesh->o_VBq = mesh->device.malloc(mesh->Np*mesh->cubNp*sizeof(dfloat),VBq);
  mesh->o_PBq = mesh->device.malloc(mesh->Np*mesh->cubNp*sizeof(dfloat),PBq);

  mesh->o_L0vals = mesh->device.malloc(mesh->Nfp*3*sizeof(dfloat),L0vals);
  mesh->o_ELids =
    mesh->device.malloc(mesh->Np*mesh->max_EL_nnz*sizeof(int),ELids);
  mesh->o_ELvals =
    mesh->device.malloc(mesh->Np*mesh->max_EL_nnz*sizeof(dfloat),ELvals);

  mesh->o_cubInterpT =
    mesh->device.malloc(mesh->Np*mesh->cubNp*sizeof(dfloat),
			cubInterpT);

  mesh->o_cubProjectT =
    mesh->device.malloc(mesh->Np*mesh->cubNp*sizeof(dfloat),
			cubProjectT);


  mesh->o_cubDrWT =
    mesh->device.malloc(mesh->Np*mesh->cubNp*sizeof(dfloat),
			cubDrWT);

  mesh->o_cubDsWT =
    mesh->device.malloc(mesh->Np*mesh->cubNp*sizeof(dfloat),
			cubDsWT);

  mesh->o_cubDtWT =
    mesh->device.malloc(mesh->Np*mesh->cubNp*sizeof(dfloat),
			cubDsWT); // dummy to align with 3d

  mesh->o_cubDWmatrices = mesh->device.malloc(2*mesh->cubNp*mesh->Np*sizeof(dfloat), cubDrsWT);

  mesh->o_intInterpT =
    mesh->device.malloc(mesh->Nfp*mesh->Nfaces*mesh->intNfp*sizeof(dfloat),
			intInterpT);

  mesh->o_intLIFTT =
    mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->intNfp*sizeof(dfloat),
			intLIFTT);

  mesh->o_intx =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp*sizeof(dfloat),
			mesh->intx);

  mesh->o_inty =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp*sizeof(dfloat),
			mesh->inty);

  mesh->o_intz =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp*sizeof(dfloat),
			mesh->inty); // dummy to align with 3d

  free(DrsT); free(ST);
    
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
  kernelInfo["defines/" "p_G11ID"]= G11ID;
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

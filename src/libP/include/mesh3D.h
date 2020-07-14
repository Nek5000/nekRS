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

#ifndef MESH3D_H 
#define MESH3D_H 1

// generic mesh structure 
#include "mesh.h"

extern "C" { // Begin C Linkage
#define mesh3D mesh_t

// mesh readers
mesh3D* meshParallelReaderTri3D(char *fileName);
mesh3D* meshParallelReaderQuad3D(char *fileName);
mesh3D* meshParallelReaderTet3D(char *fileName);
mesh3D* meshParallelReaderHex3D(char *fileName);

// build connectivity in serial
void meshConnect3D(mesh3D *mesh);

// build element-boundary connectivity
void meshConnectBoundary3D(mesh3D *mesh);

// build connectivity in parallel
void meshParallelConnect3D(mesh3D *mesh);

// repartition elements in parallel
void meshGeometricPartition3D(mesh3D *mesh);

// print out mesh 
void meshPrint3D(mesh3D *mesh);

// print out mesh in parallel from the root process
void meshParallelPrint3D(mesh3D *mesh);

// print out mesh partition in parallel
void meshVTU3D(mesh3D *mesh, char *fileName);

// print out mesh field
void meshPlotVTU3D(mesh3D *mesh, char *fileNameBase, int fld);
void meshPlotContour3D(mesh_t *mesh, char *fname, dfloat *u, int Nlevels, dfloat *levels);
void meshPlotAdaptiveContour3D(mesh_t *mesh, char *fname, dfloat *u, int Nlevels, dfloat *levels, dfloat tol);

// compute geometric factors for local to physical map
void meshGeometricFactorsTri3D(mesh3D *mesh);
void meshGeometricFactorsQuad3D(mesh3D *mesh);
void meshGeometricFactorsTet3D(mesh3D *mesh);
void meshGeometricFactorsHex3D(mesh3D *mesh);

void meshSurfaceGeometricFactorsTri3D(mesh3D *mesh);
void meshSurfaceGeometricFactorsQuad3D(mesh3D *mesh);
void meshSurfaceGeometricFactorsTet3D(mesh3D *mesh);
void meshSurfaceGeometricFactorsHex3D(mesh3D *mesh);

void meshPhysicalNodesTri3D(mesh3D *mesh);
void meshPhysicalNodesQuad3D(mesh3D *mesh);
void meshPhysicalNodesTet3D(mesh3D *mesh);
void meshPhysicalNodesHex3D(mesh3D *mesh);

void meshLoadReferenceNodesTet3D(mesh3D *mesh, int N);
void meshLoadReferenceNodesHex3D(mesh3D *mesh, int N);

void meshGradientTet3D(mesh3D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy, dfloat *dqdz);
void meshGradientHex3D(mesh3D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy, dfloat *dqdz);

// print out parallel partition i
void meshPartitionStatistics3D(mesh3D *mesh);

// default occa set up
void meshOccaSetup3D(mesh3D *mesh, setupAide &newOptions, occa::properties &kernelInfo);
void meshOccaSetupQuad3D(mesh_t *mesh, setupAide &newOptions, occa::properties &kernelInfo);
void meshOccaSetupTri3D(mesh_t *mesh, setupAide &newOptions, occa::properties &kernelInfo);

void meshOccaPopulateDevice3D(mesh3D *mesh, setupAide &newOptions, occa::properties &kernelInfo);
void meshOccaCloneDevice(mesh_t *donorMesh, mesh_t *mesh);

// functions that call OCCA kernels
void occaTest3D(mesh3D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy, dfloat *dqdz);

// 
void occaOptimizeGradientTet3D(mesh3D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy, dfloat *dqdz);
void occaOptimizeGradientHex3D(mesh3D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy, dfloat *dqdz);

// serial face-node to face-node connection
void meshConnectFaceNodes3D(mesh3D *mesh);

//
mesh3D *meshSetupTri3D(char *filename, int N, dfloat sphereRadius);
mesh3D *meshSetupQuad3D(char *filename, int N, dfloat sphereRadius);
mesh3D *meshSetupTet3D(char *filename, int N);
mesh3D *meshSetupHex3D(char *filename, int N);

void meshParallelConnectNodesHex3D(mesh3D *mesh);

// halo connectivity information
void meshHaloSetup3D(mesh3D *mesh);

// perform halo exchange
void meshHaloExchange3D(mesh3D *mesh,
			size_t Nbytes,  // number of bytes per element
			void *sourceBuffer, 
			void *sendBuffer, 
			void *recvBuffer);

void meshHaloExchangeStart3D(mesh3D *mesh,
			     size_t Nbytes,       // message size per element
			     void *sendBuffer,    // temporary buffer
			     void *recvBuffer);

void meshHaloExchangeFinish3D(mesh3D *mesh);

// build list of nodes on each face of the reference element
void meshBuildFaceNodes3D(mesh3D *mesh);
void meshBuildFaceNodesHex3D(mesh3D *mesh);



dfloat meshMRABSetup3D(mesh3D *mesh, dfloat *EToDT, int maxLevels, dfloat finalTime); 

//MRAB weighted mesh partitioning
void meshMRABWeightedPartition3D(mesh3D *mesh, dfloat *weights,
                                      int numLevels, int *levels);

void interpolateHex3D(dfloat *Inter, dfloat *x, int N, dfloat *Ix, int M);

#define norm3(a,b,c) ( sqrt((a)*(a)+(b)*(b)+(c)*(c)) )

/* offsets for geometric factors */
#define RXID 0  
#define RYID 1  
#define SXID 2  
#define SYID 3
#define  JID 4
#define JWID 5
#define IJWID 6
#define RZID 7
#define SZID 8  
#define TXID 9  
#define TYID 10  
#define TZID 11  



/* offsets for second order geometric factors */
#define G00ID 0  
#define G01ID 1  
#define G11ID 2
#define GWJID 3  
#define G12ID 4
#define G02ID 5
#define G22ID 6  


/* offsets for nx, ny, sJ, 1/J */
#define NXID 0  
#define NYID 1  
#define SJID 2
#define IJID 3
#define IHID 4
#define WSJID 5
#define WIJID 6
#define NZID 7
#define STXID 8
#define STYID 9  
#define STZID 10 
#define SBXID 11 
#define SBYID 12 
#define SBZID 13
#define SURXID 14
#define SURYID 15
#define SURZID 16
//
//offsets for boltzmann PML variables
#define QXID1 0  
#define QXID2 1  
#define QXID3 2
#define QXID4 3  
#define QXID5 4  
#define QXID6 5  
#define QXID8 6 
//
#define QYID1 7  
#define QYID2 8  
#define QYID3 9
#define QYID4 10  
#define QYID5 11  
#define QYID7 12  
#define QYID9 13 
//
#define QZID1 14  
#define QZID2 15  
#define QZID3 16
#define QZID4 17  
#define QZID6 18  
#define QZID7 19  
#define QZID10  20   

mesh3D *meshSetupBoxHex3D(int N, setupAide &options);
void meshConnectPeriodicFaceNodes3D(mesh3D *mesh, dfloat xper, dfloat yper, dfloat zper);

} // end C Linkage
#endif


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

#ifndef MESH2D_H 
#define MESH2D_H 1

#include "mesh.h"

extern "C" { // Begin C Linkage
// will eventually rename mesh2D to mesh_t in src
#define mesh2D mesh_t

mesh2D* meshReaderTri2D(char *fileName);
mesh2D* meshReaderQuad2D(char *fileName);

// mesh readers
mesh2D* meshParallelReaderTri2D(char *fileName);
mesh2D* meshParallelReaderQuad2D(char *fileName);

// build connectivity in serial
void meshConnect2D(mesh2D *mesh);

// build element-boundary connectivity
void meshConnectBoundary2D(mesh2D *mesh);

// build connectivity in parallel
void meshParallelConnect2D(mesh2D *mesh);

// build global connectivity in parallel
void meshParallelConnectNodesQuad2D(mesh2D *mesh);

// create global number of nodes
void meshNumberNodes2D(mesh2D *mesh);

// repartition elements in parallel
void meshGeometricPartition2D(mesh2D *mesh);

// print out mesh 
void meshPrint2D(mesh2D *mesh);

// print out mesh in parallel from the root process
void meshParallelPrint2D(mesh2D *mesh);

// print out mesh partition in parallel
void meshVTU2D(mesh2D *mesh, char *fileName);

// print out solution at plot nodes 
void meshPlotVTU2D(mesh2D *mesh, char *fileNameBase, int fld);

// compute geometric factors for local to physical map 
void meshGeometricFactorsTri2D(mesh2D *mesh);
void meshGeometricFactorsQuad2D(mesh2D *mesh);

void meshSurfaceGeometricFactorsTri2D(mesh2D *mesh);
void meshSurfaceGeometricFactorsQuad2D(mesh2D *mesh);

void meshPhysicalNodesTri2D(mesh2D *mesh);
void meshPhysicalNodesQuad2D(mesh2D *mesh);

void meshLoadReferenceNodesTri2D(mesh2D *mesh, int N);
void meshLoadReferenceNodesQuad2D(mesh2D *mesh, int N);

void meshGradientTri2D(mesh2D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy);
void meshGradientQuad2D(mesh2D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy);

// print out parallel partition i
void meshPartitionStatistics2D(mesh2D *mesh);

// functions that call OCCA kernels
void occaTest(mesh2D *mesh);

// 
void occaOptimizeGradientTri2D(mesh2D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy);
void occaOptimizeGradientQuad2D(mesh2D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy);

// serial face-node to face-node connection
void meshConnectFaceNodes2D(mesh2D *mesh);

// serial face-mode to face-mode connection
void meshConnectFaceModes2D(mesh2D *mesh, int *faceModes, dfloat *V);

// halo connectivity information
void meshHaloSetup2D(mesh2D *mesh);

// perform complete halo exchange
void meshHaloExchange2D(mesh2D *mesh,
			size_t Nbytes, // number of bytes per element
			void *sourceBuffer, 
			void *sendBuffer, 
			void *recvBuffer);

// start halo exchange
void meshHaloExchangeStart2D(mesh2D *mesh,
			     size_t Nbytes,       // message size per element
			     void *sendBuffer,    // outgoing halo
			     void *recvBuffer);   // incoming halo

// finish halo exchange
void meshHaloExchangeFinish2D(mesh2D *mesh);

// extract halo data from sourceBuffer and save to sendBuffer
void meshHaloExtract2D(mesh2D *mesh, size_t Nbytes, void *sourceBuffer, void *sendBuffer);

// build list of nodes on each face of the reference element
void meshBuildFaceNodesTri2D(mesh2D *mesh);
void meshBuildFaceNodesQuad2D(mesh2D *mesh);

mesh2D *meshSetupTri2D(char *filename, int N);
mesh2D *meshSetupQuad2D(char *filename, int N);

// set up OCCA device and copy generic element info to device
void meshOccaSetup2D(mesh2D *mesh, setupAide &newOptions, occa::properties &kernelInfo);

// void meshMRABSetup2D(mesh2D *mesh, dfloat *EToDT, int maxLevels); 
dfloat meshMRABSetup2D(mesh2D *mesh, dfloat *EToDT, int maxLevels, dfloat finalTime); 


//MRAB weighted mesh partitioning
void meshMRABWeightedPartition2D(mesh2D *mesh, dfloat *weights,
                                      int numLevels, int *levels);


// Setup probe information
// Probe Setup : AK
void meshProbeSetup2D(mesh2D *mesh, dfloat *pX, dfloat *pY);
void meshVandermonde2D(int N, int sizeR, dfloat *r, dfloat *s, dfloat *V);
dfloat meshSimplex2D(dfloat a, dfloat b, int i, int j);
dfloat meshJacobiP(dfloat a, dfloat alpha, dfloat beta, int N);
dfloat meshFactorial(int n);


#define norm2(a,b) ( sqrt((a)*(a)+(b)*(b)) )


/* offsets for geometric factors */
#define RXID 0  
#define RYID 1  
#define SXID 2  
#define SYID 3  
#define  JID 4
#define JWID 5
#define IJWID 6


/* offsets for second order geometric factors */
#define G00ID 0  
#define G01ID 1  
#define G11ID 2
#define GWJID 3

/* offsets for nx, ny, sJ, 1/J */
#define NXID 0  
#define NYID 1  
#define SJID 2  
#define IJID 3  
#define IHID 4
#define WSJID 5
#define WIJID 6


mesh2D *meshSetupBoxQuad2D(int N, setupAide &options);
void meshConnectPeriodicFaceNodes2D(mesh2D *mesh, dfloat xper, dfloat yper);

} // end C Linkage
#endif


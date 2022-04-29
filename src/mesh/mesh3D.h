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

// build connectivity in serial
void meshConnect3D(mesh_t *mesh);

// build connectivity in parallel
void meshParallelConnect3D(mesh_t *mesh);

// compute geometric factors for local to physical map
void meshGeometricFactorsHex3D(mesh_t *mesh);

void meshSurfaceGeometricFactorsHex3D(mesh_t *mesh);

void meshPhysicalNodesHex3D(mesh_t *mesh);

void meshLoadReferenceNodesHex3D(mesh_t *mesh, int N, int cubN);

// default occa set up
void meshOccaSetup3D(mesh_t *mesh, setupAide &newOptions, occa::properties &kernelInfo);

void meshOccaPopulateDeviceHex3D(mesh_t *mesh, setupAide &newOptions, occa::properties &kernelInfo);
void meshOccaCloneDevice(mesh_t* donorMesh, mesh_t* mesh);

// serial face-node to face-node connection
void meshConnectFaceNodes3D(mesh_t *mesh);

void meshParallelConnectNodesHex3D(mesh_t *mesh);

// halo connectivity information
void meshHaloSetup3D(mesh_t *mesh);

// perform halo exchange
void meshHaloExchange3D(mesh_t *mesh,
                        size_t Nbytes, // number of bytes per element
                        void *sourceBuffer,
                        void *sendBuffer,
                        void *recvBuffer);

void meshHaloExchangeStart3D(mesh_t *mesh,
                             size_t Nbytes,    // message size per element
                             void *sendBuffer, // temporary buffer
                             void *recvBuffer);

void meshHaloExchangeFinish3D(mesh_t *mesh);

// build list of nodes on each face of the reference element
void meshBuildFaceNodes3D(mesh_t *mesh);
void meshBuildFaceNodesHex3D(mesh_t *mesh);
void interpolateHex3D(dfloat* Inter, dfloat* x, int N, dfloat* Ix, int M);

/* offsets for geometric factors */
#define RXID 0
#define RYID 1
#define RZID 2
#define SXID 3
#define SYID 4
#define SZID 5
#define TXID 6
#define TYID 7
#define TZID 8
#define JID 9
#define JWID 10
#define IJWID 11

/* offsets for second order geometric factors */
#define G00ID 0
#define G01ID 1
#define G11ID 2
#define G12ID 3
#define G02ID 4
#define G22ID 5
#define GWJID 6

/* offsets for nx, ny, sJ, 1/J */
#define NXID 0
#define NYID 1
#define NZID 2

// tangentails
#define T1XID 3
#define T1YID 4
#define T1ZID 5

#define T2XID 6
#define T2YID 7
#define T2ZID 8

#define SJID 9
#define IJID 10
#define WIJID 11
#define WSJID 12

// Mesh generation
void NodesHex3D(int _N, dfloat* _r, dfloat* _s, dfloat* _t);
void FaceNodesHex3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t, int *_faceNodes);
#endif

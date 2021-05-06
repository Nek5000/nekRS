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

#ifndef MESH_H
#define MESH_H 1

#include <unistd.h>
#include <assert.h>

#include <math.h>
#include <stdlib.h>

#include "nrssys.hpp"
#include "ogs.hpp"
#include "linAlg.hpp"

#define TRIANGLES 3
#define QUADRILATERALS 4
#define TETRAHEDRA 6
#define HEXAHEDRA 12

struct mesh_t
{
  void move();
  void update();
  void computeInvLMM();

  int nAB;
  dfloat* coeffAB; // coefficients for AB integration
  occa::memory o_coeffAB;
  int dim;
  int Nverts, Nfaces, NfaceVertices;

  int cht;

  mesh_t* fluid;

  hlong Nnodes;
  dfloat* EX; // coordinates of vertices for each element
  dfloat* EY;
  dfloat* EZ;

  dlong Nelements;
  dlong fieldOffset;
  dlong Nlocal;
  hlong* EToV; // element-to-vertex connectivity
  dlong* EToE; // element-to-element connectivity
  int* EToF;   // element-to-(local)face connectivity
  int* EToP;   // element-to-partition/process connectivity
  int* EToB;   // element-to-boundary condition type

  dlong* elementInfo; //type of element
  occa::memory o_elementInfo;

  // boundary faces
  hlong NboundaryFaces; // number of boundary faces
  hlong* boundaryInfo; // list of all boundary faces (type, vertex-1, vertex-2, vertex-3) in the mesh

  // MPI halo exchange info
  dlong totalHaloPairs;   // number of elements to be sent in halo exchange
  dlong* haloElementList; // sorted list of elements to be sent in halo exchange
  int* NhaloPairs;      // number of elements worth of data to send/recv
  int NhaloMessages;      // number of messages to send

  dlong* haloGetNodeIds; // volume node ids of outgoing halo nodes
  dlong* haloPutNodeIds; // volume node ids of incoming halo nodes

  void* haloSendRequests;
  void* haloRecvRequests;

  // CG gather-scatter info
  hlong* globalIds;
  void* gsh, * hostGsh; // gslib struct pointer
  ogs_t* ogs; //occa gs pointer
  oogs_t* oogs; //occa gs pointer

  // list of elements that are needed for global gather-scatter
  dlong NglobalGatherElements;
  dlong* globalGatherElementList;
  occa::memory o_globalGatherElementList;

  // list of elements that are not needed for global gather-scatter
  dlong NlocalGatherElements;
  dlong* localGatherElementList;
  occa::memory o_localGatherElementList;

  // volumeGeometricFactors;
  dlong Nvgeo;
  dfloat* vgeo;

  // second order volume geometric factors
  dlong Nggeo;
  dfloat* ggeo;

  // volume node info
  int N, Np;
  dfloat* r, * s, * t;    // coordinates of local nodes
  dfloat* MM;
  dfloat* LMM, * invLMM;
  dfloat* x, * y, * z;    // coordinates of physical nodes

  dfloat volume;

  // indices of vertex nodes
  int* vertexNodes;

  // indices of edge nodes
  int* edgeNodes;

  int NedgeNodes;

  // quad specific quantity
  int Nq, NqP, NpP;

  dfloat* D; // 1D differentiation matrix (for tensor-product)
  dfloat* DW; // weak 1D differentiation matrix (for tensor-product)
  dfloat* gllz; // 1D GLL quadrature nodes
  dfloat* gllw; // 1D GLL quadrature weights

  // face node info
  int Nfp;        // number of nodes per face
  int* faceNodes; // list of element reference interpolation nodes on element faces
  dlong* vmapM;     // list of volume nodes that are face nodes
  dlong* vmapP;     // list of volume nodes that are paired with face nodes
  dlong* mapP;     // list of surface nodes that are paired with -ve surface  nodes
  int* faceVertices; // list of mesh vertices on each face

  dlong Nsgeo;
  dfloat* sgeo;

  // field info for PDE solver
  int Nfields;
  // cubature
  int cubNp, cubNfp, cubNq;
  dfloat* cubr, * cubs, * cubt, * cubw; // coordinates and weights of local cubature nodes
  dfloat* cubx, * cuby, * cubz;    // coordinates of physical nodes
  dfloat* cubInterp; // interpolate from W&B to cubature nodes
  dfloat* cubProject; // projection matrix from cubature nodes to W&B nodes
  dfloat* cubD;       // 1D differentiation matrix
  dfloat* cubDiffInterp;     // 1D weak differentiation matrix
  dfloat* cubDW;     // 1D weak differentiation matrix
  dfloat* cubDWmatrices;

  dfloat* cubvgeo;  //volume geometric data at cubature points

  dfloat* interpRaise;
  dfloat* interpLower;

  // surface integration node info
  int intNfp;       // number of integration nodes on each face
  dfloat* intInterp; // interp from surface node to integration nodes
  dfloat* intx, * inty, * intz; // coordinates of suface integration nodes

  occa::memory o_LMM, o_invLMM, o_divU;

  // mesh velocity
  occa::memory o_U;
  dfloat* U; // host shadow of mesh velocity

  occa::memory o_D;
  occa::memory o_DPfloat;

  occa::memory o_DW; // tensor product differentiation matrix (for Hexes)
  occa::memory o_DT;
  occa::memory o_DTPfloat;

  occa::memory o_vgeo, o_sgeo;
  occa::memory o_vmapM, o_vmapP, o_mapP;

  occa::memory o_EToB, o_x, o_y, o_z;

  // cubature (for wadg)
  occa::memory o_cubDWT, o_cubD;
  occa::memory o_cubDiffInterpT;
  occa::memory o_cubDWmatrices;
  occa::memory o_cubInterpT, o_cubProjectT;

  occa::memory o_cubvgeo;

  // DG halo exchange info
  occa::memory o_haloElementList;
  occa::memory o_haloBuffer;
  occa::memory o_haloGetNodeIds;
  occa::memory o_haloPutNodeIds;

  occa::memory o_ggeo; // second order geometric factors
  occa::memory o_ggeoPfloat; // second order geometric factors

  occa::memory o_gllw;
  occa::memory o_cubw;
  occa::memory o_faceNodes;

  occa::kernel haloExtractKernel;

  occa::kernel maskKernel;
  occa::kernel maskPfloatKernel;

  occa::kernel geometricFactorsKernel;
  occa::kernel surfaceGeometricFactorsKernel;
  occa::kernel nStagesSumVectorKernel;
  occa::kernel velocityDirichletKernel;
};

occa::properties populateMeshProperties(mesh_t*);
// serial sort
void mysort(hlong* data, int N, const char* order);

// sort entries in an array in parallel
void parallelSort(int size, int rank, MPI_Comm comm,
                  int N, void* vv, size_t sz,
                  int (* compare)(const void*, const void*),
                  void (* match)(void*, void*)
                  );

#define mymax(a,b) (((a) > (b))?(a):(b))
#define mymin(a,b) (((a) < (b))?(a):(b))

/* dimension independent mesh operations */
void meshConnect(mesh_t* mesh);

/* build parallel face connectivity */
void meshParallelConnect(mesh_t* mesh);

/* build global connectivity in parallel */
void meshGlobalIds(mesh_t* mesh);

void meshHaloSetup(mesh_t* mesh);
void meshHaloPhysicalNodes(mesh_t* mesh);

/* extract whole elements for the halo exchange */
void meshHaloExtract(mesh_t* mesh, size_t Nbytes, void* sourceBuffer, void* haloBuffer);

void meshHaloExchange(mesh_t* mesh,
                      size_t Nbytes, // message size per element
                      void* sourceBuffer,
                      void* sendBuffer, // temporary buffer
                      void* recvBuffer);

void meshHaloExchangeStart(mesh_t* mesh,
                           size_t Nbytes, // message size per element
                           void* sendBuffer, // temporary buffer
                           void* recvBuffer);

void meshHaloExchangeFinish(mesh_t* mesh);

void meshHaloExchangeBlocking(mesh_t* mesh,
                              size_t Nbytes, // message size per element
                              void* sendBuffer, // temporary buffer
                              void* recvBuffer);

// print out parallel partition i
void meshPartitionStatistics(mesh_t* mesh);

// build element-boundary connectivity
void meshConnectBoundary(mesh_t* mesh);

void meshParallelGatherScatterSetup(mesh_t* mesh,
                                    dlong N,
                                    hlong* globalIds,
                                    MPI_Comm &comm,
                                    int verbose);

// generic mesh setup
mesh_t* meshSetup(char* filename, int N, setupAide &options);
void meshFree(mesh_t*);

void occaTimerTic(occa::device device,std::string name);
void occaTimerToc(occa::device device,std::string name);

extern "C"
{
void* xxtSetup(uint num_local_rows,
               void* row_ids,
               uint nnz,
               void*   A_i,
               void*   A_j,
               void* A_vals,
               int null_space,
               const char* inttype,
               const char* floattype);

void xxtSolve(void* x,
              void* A,
              void* rhs);

void xxtFree(void* A);
}

extern "C"
{
void dgesv_ ( int* N, int* NRHS, double* A,
              int* LDA,
              int* IPIV,
              double* B,
              int* LDB,
              int* INFO );

// void dgemm_(const char *TRANSA, const char *TRANSB, const int *M,
//             const int *N, const int *K, double *ALPHA, double *A, const int *LDA, double *B,
//             const int *LDB, double *BETA, double *C, const int *LDC);

void dgemm_ (char*, char*, int*, int*, int*,
             const dfloat*, const dfloat* __restrict, int*,
             const dfloat* __restrict, int*,
             const dfloat*, dfloat* __restrict, int*);

void sgesv_(int* N, int* NRHS,float* A, int* LDA, int* IPIV, float* B, int* LDB,int* INFO);

void dgetrf_(int* M, int* N, double* A, int* lda, int* IPIV, int* INFO);
void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
void dgeev_(char* JOBVL, char* JOBVR, int* N, double* A, int* LDA, double* WR, double* WI,
            double* VL, int* LDVL, double* VR, int* LDVR, double* WORK, int* LWORK, int* INFO );

double dlange_(char* NORM, int* M, int* N, double* A, int* LDA, double* WORK);
void dgecon_(char* NORM, int* N, double* A, int* LDA, double* ANORM,
             double* RCOND, double* WORK, int* IWORK, int* INFO );
}

void meshApplyElementMatrix(mesh_t* mesh, dfloat* A, dfloat* q, dfloat* Aq);
void meshApplyVectorElementMatrix(mesh_t* mesh, int Nfield, const dlong offset, dfloat* A, dfloat* q, dfloat* Aq);

void meshRecursiveSpectralBisectionPartition(mesh_t* mesh);

void matrixInverse(int N, dfloat* A);
dfloat matrixConditionNumber(int N, dfloat* A);

void matrixRightSolve(int NrowsA, int NcolsA, dfloat* A, int NrowsB, int NcolsB, dfloat* B, dfloat* C);
void matrixEig(int N, dfloat* A, dfloat* VR, dfloat* WR, dfloat* WI);
void matrixTranspose(const int M, const int N,
                     const dfloat* A, const int LDA,
                     dfloat* AT, const int LDAT);

// 1D mesh basis functions
void Nodes1D(int _N, dfloat* _r);
void EquispacedNodes1D(int _N, dfloat* _r);
void OrthonormalBasis1D(dfloat a, int i, dfloat* P);
void GradOrthonormalBasis1D(dfloat a, int i, dfloat* Pr);
void Vandermonde1D(int _N, int Npoints, dfloat* _r, dfloat* V);
void GradVandermonde1D(int _N, int Npoints, dfloat* _r, dfloat* Vr);
void MassMatrix1D(int _Np, dfloat* V, dfloat* _MM);
void Dmatrix1D(int _N, int NpointsIn, dfloat* _rIn,
               int NpointsOut, dfloat* _rOut, dfloat* _Dr);
void DWmatrix1D(int _N, dfloat* _D, dfloat* _DT);

void InterpolationMatrix1D(int _N,
                           int NpointsIn, dfloat* rIn,
                           int NpointsOut, dfloat* rOut,
                           dfloat* I);
void DegreeRaiseMatrix1D(int Nc, int Nf, dfloat* P);
void CubatureWeakDmatrix1D(int _Nq, int _cubNq,
                           dfloat* _cubProject, dfloat* _cubD, dfloat* _cubPDT);
dfloat JacobiP(dfloat a, dfloat alpha, dfloat beta, int _N);
dfloat GradJacobiP(dfloat a, dfloat alpha, dfloat beta, int _N);
void JacobiGLL(int _N, dfloat* _x, dfloat* _w = nullptr);
void JacobiGQ(dfloat alpha, dfloat beta, int _N, dfloat* _x, dfloat* _w);
#endif

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

#include "nekrsSys.hpp"
#include "ogs.hpp"
#include "linAlg.hpp"

#define TRIANGLES 3
#define QUADRILATERALS 4
#define TETRAHEDRA 6
#define HEXAHEDRA 12

struct mesh_t {

  // Distance pseudo function
  // type refers to the type of distance function
  //   type = "cheap_dist" : nek5000-style cheap_dist
  // other types are not yet supported

  // distance: for each boundary id, compute the distance from the boundary
  // returns on output a vector field with nbID entries, where each entry is the distance to the boundary id
  occa::memory
  distance(int nbID, const occa::memory &o_bID, dlong offsetFld, std::string type, int maxIter = 10000);
  std::vector<dfloat>
  distance(const std::vector<dlong> &bID, dlong offsetFld, std::string type, int maxIter = 10000);

  // minDistance: compute minimum distance across all boundary ids
  // returns a single distance field
  occa::memory minDistance(int nbID, const occa::memory &o_bID, std::string type, int maxIter = 10000);
  std::vector<dfloat> minDistance(const std::vector<dlong> &bID, std::string type, int maxIter = 10000);

  occa::memory intpMatrix(std::vector<dfloat> M);
  void interpolate(const occa::memory& o_z, mesh_t *meshC, occa::memory& o_zC, bool uniform = false, dlong nel = 0);

  void move();
  void update();

  void geometricFactors();
  void surfaceGeometricFactors();

  void computeInvLMM();

  std::vector<dfloat> surfaceAreaMultiplyIntegrate(int nbID, 
                                                   const occa::memory &o_bID, 
                                                   const occa::memory &o_fld);

  std::vector<dfloat> surfaceAreaMultiplyIntegrate(int Nfields,
                                                   dlong fieldOffset,
                                                   int nbID,
                                                   const occa::memory &o_bID,
                                                   const occa::memory &o_fld);

  std::vector<dfloat> surfaceAreaNormalMultiplyIntegrate(dlong fieldOffset,
                                                         int nbID,
                                                         const occa::memory &o_bID,
                                                         const occa::memory &o_fld);

  occa::memory surfaceAreaMultiply(int nbID, const occa::memory &o_bID, const occa::memory &o_fld);

  std::tuple<std::vector<dfloat>, std::vector<dfloat>, std::vector<dfloat>> xyzHost() const 
  {
    std::vector<dfloat> x(Nlocal);
    std::vector<dfloat> y(Nlocal);
    std::vector<dfloat> z(Nlocal);

    o_x.copyTo(x.data(), Nlocal);
    o_y.copyTo(y.data(), Nlocal);
    o_z.copyTo(z.data(), Nlocal);

    return {x, y, z}; 
  };

  int nAB;
  dfloat *coeffAB; // coefficients for AB integration
  occa::memory o_coeffAB;
  int dim;
  int Nverts, Nfaces, NfaceVertices;

  int Nbid;

  int cht;

  hlong Nnodes;
  dfloat *EX; // coordinates of vertices for each element
  dfloat *EY;
  dfloat *EZ;

  dlong fieldOffset;

  dlong Nelements;
  hlong NelementsGlobal;

  dlong Nlocal;
  hlong Nglobal; // global T-vector size 

  hlong NboundaryFaces;
  hlong *EToV; // element-to-vertex connectivity
  dlong *EToE; // element-to-element connectivity
  int *EToF;   // element-to-(local)face connectivity
  int *EToP;   // element-to-partition/process connectivity
  int *EToB;   // element-to-boundary condition type

  dlong *elementInfo; // type of element
  occa::memory o_elementInfo;

  hlong *globalIds;
  ogs_t *ogs;
  oogs_t *oogs;

  // list of all elements
  // elementList[e] = e
  dlong *elementList;
  occa::memory o_elementList;

  // list of elements that are needed for global gather-scatter
  dlong NglobalGatherElements;
  dlong *globalGatherElementList;
  occa::memory o_globalGatherElementList;

  // list of elements that are not needed for global gather-scatter
  dlong NlocalGatherElements;
  dlong *localGatherElementList;
  occa::memory o_localGatherElementList;

  // volumeGeometricFactors;
  dlong Nvgeo;

  // second order volume geometric factors
  dlong Nggeo;

  // volume node info
  static constexpr int maxNqIntp = 16; 
  int N, Np;
  dfloat *r, *s, *t; // coordinates of local nodes
  dfloat *MM;

  dfloat volume;

  // indices of vertex nodes
  int *vertexNodes;

  // indices of edge nodes
  int *edgeNodes;

  int NedgeNodes;

  // quad specific quantity
  int Nq, NqP, NpP;

  dfloat *D;    // 1D differentiation matrix (for tensor-product)
  dfloat *DW;   // weak 1D differentiation matrix (for tensor-product)
  dfloat *gllz; // 1D GLL quadrature nodes
  dfloat *gllw; // 1D GLL quadrature weights

  // face node info
  int Nfp;           // number of nodes per face
  int *faceNodes;    // list of element reference interpolation nodes on element faces
  dlong *vmapM;      // list of volume nodes that are face nodes
  int *faceVertices; // list of mesh vertices on each face

  dlong Nsgeo;

  // field info for PDE solver
  int Nfields = -1;

  // cubature
  int cubNp, cubNfp, cubNq;
  dfloat *cubr, *cubs, *cubt, *cubw; // coordinates and weights of local cubature nodes
  dfloat *cubx, *cuby, *cubz;        // coordinates of physical nodes
  dfloat *cubInterp;                 // interpolate from W&B to cubature nodes
  dfloat *cubProject;                // projection matrix from cubature nodes to W&B nodes
  dfloat *cubD;                      // 1D differentiation matrix
  dfloat *cubDiffInterp;             // 1D weak differentiation matrix
  dfloat *cubDW;                     // 1D weak differentiation matrix
  dfloat *cubDWmatrices;

  dfloat *interpRaise;
  dfloat *interpLower;

  // surface integration node info
  int intNfp;                 // number of integration nodes on each face
  dfloat *intInterp;          // interp from surface node to integration nodes
  dfloat *intx, *inty, *intz; // coordinates of suface integration nodes

  occa::memory o_divU;
  occa::memory o_LMM, o_invLMM;

  occa::memory& o_Jw = o_LMM;
  occa::memory& o_invAJw = o_invLMM;
  
  occa::memory o_invAJwTimesInvDegree;

  occa::memory o_U;
  occa::memory o_Ue;

  occa::memory o_D;

  occa::memory o_DW;
  occa::memory o_DT;

  occa::memory o_vgeo, o_sgeo;
  occa::memory o_vmapM, o_vmapP, o_mapP;

  occa::memory o_EToB, o_x, o_y, o_z;

  // cubature (for wadg)
  occa::memory o_cubDWT, o_cubD;
  occa::memory o_cubDiffInterpT;
  occa::memory o_cubDWmatrices;
  occa::memory o_cubInterpT, o_cubProjectT;

  occa::memory o_cubvgeo;

  occa::memory o_ggeo; // second order geometric factors

  occa::memory o_gllz;
  occa::memory o_gllw;
  occa::memory o_cubw;
  occa::memory o_faceNodes;

  occa::kernel haloExtractKernel;

  std::array<occa::kernel, maxNqIntp> intpKernel;

  occa::kernel geometricFactorsKernel;
  occa::kernel surfaceGeometricFactorsKernel;
  occa::kernel cubatureGeometricFactorsKernel;
  occa::kernel nStagesSumVectorKernel;
  occa::kernel velocityDirichletKernel;

  occa::kernel setBIDKernel;
  occa::kernel distanceKernel;
  occa::kernel hlongSumKernel;
};

std::pair<mesh_t*, mesh_t*> createMesh(MPI_Comm comm, int N, int cubN, bool cht, occa::properties &kernelInfo);
mesh_t *createMeshMG(mesh_t *_mesh, int Nc);

occa::properties meshKernelProperties(int N);
// serial sort
void mysort(hlong *data, int N, const char *order);

// sort entries in an array in parallel
void parallelSort(int size,
                  int rank,
                  MPI_Comm comm,
                  int N,
                  void *vv,
                  size_t sz,
                  int (*compare)(const void *, const void *),
                  void (*match)(void *, void *));

/* dimension independent mesh operations */
void meshConnect(mesh_t *mesh);

/* build parallel face connectivity */
void meshParallelConnect(mesh_t *mesh);

/* build global connectivity in parallel */
void meshGlobalIds(mesh_t *mesh);

void planarAvg(mesh_t *mesh,
               const std::string &dir,
               int NELGX,
               int NELGY,
               int NELGZ,
               int nflds,
               dlong offset,
               occa::memory &o_avg);

void meshPartitionStatistics(mesh_t *mesh);

void meshParallelGatherScatterSetup(mesh_t *mesh,
                                    dlong N,
                                    hlong *globalIds,
                                    MPI_Comm &comm,
                                    oogs_mode gsMode,
                                    int verbose);

void meshFree(mesh_t *);

void printMeshMetrics(mesh_t *mesh);

extern "C" {
void dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);

// void dgemm_(const char *TRANSA, const char *TRANSB, const int *M,
//             const int *N, const int *K, double *ALPHA, double *A, const int *LDA, double *B,
//             const int *LDB, double *BETA, double *C, const int *LDC);

void dgemm_(char *,
            char *,
            int *,
            int *,
            int *,
            const double *,
            const double *__restrict,
            int *,
            const double *__restrict,
            int *,
            const double *,
            double *__restrict,
            int *);

void sgesv_(int *N, int *NRHS, float *A, int *LDA, int *IPIV, float *B, int *LDB, int *INFO);

void dgetrf_(int *M, int *N, double *A, int *lda, int *IPIV, int *INFO);
void dgetri_(int *N, double *A, int *lda, int *IPIV, double *WORK, int *lwork, int *INFO);
void dgeev_(char *JOBVL,
            char *JOBVR,
            int *N,
            double *A,
            int *LDA,
            double *WR,
            double *WI,
            double *VL,
            int *LDVL,
            double *VR,
            int *LDVR,
            double *WORK,
            int *LWORK,
            int *INFO);

double dlange_(char *NORM, int *M, int *N, double *A, int *LDA, double *WORK);
void dgecon_(char *NORM,
             int *N,
             double *A,
             int *LDA,
             double *ANORM,
             double *RCOND,
             double *WORK,
             int *IWORK,
             int *INFO);
}

void meshApplyElementMatrix(mesh_t *mesh, dfloat *A, dfloat *q, dfloat *Aq);
void meshApplyVectorElementMatrix(mesh_t *mesh,
                                  int Nfield,
                                  const dlong offset,
                                  dfloat *A,
                                  dfloat *q,
                                  dfloat *Aq);

void meshRecursiveSpectralBisectionPartition(mesh_t *mesh);

dfloat matrixConditionNumber(int N, dfloat *A);

void matrixRightSolve(int NrowsA, int NcolsA, dfloat *A, int NrowsB, int NcolsB, dfloat *B, dfloat *C);
void matrixEig(int N, dfloat *A, dfloat *VR, dfloat *WR, dfloat *WI);

void meshLoadKernels(mesh_t *mesh);

void OrthonormalBasisHex3D(dfloat a, dfloat b, dfloat c, int i, int j, int k, dfloat* P);
void VandermondeHex3D(int _N, int Npoints, dfloat* _r, dfloat* _s, dfloat* _t, dfloat* V);

// 1D mesh basis functions
void Nodes1D(int _N, dfloat *_r);
void EquispacedNodes1D(int _N, dfloat *_r);
void OrthonormalBasis1D(dfloat a, int i, dfloat *P);
void GradOrthonormalBasis1D(dfloat a, int i, dfloat *Pr);
void Vandermonde1D(int _N, int Npoints, dfloat *_r, dfloat *V);
void GradVandermonde1D(int _N, int Npoints, dfloat *_r, dfloat *Vr);
void MassMatrix1D(int _Np, dfloat *V, dfloat *_MM);
void Dmatrix1D(int _N, int NpointsIn, dfloat *_rIn, int NpointsOut, dfloat *_rOut, dfloat *_Dr);
void DWmatrix1D(int _N, dfloat *_D, dfloat *_DT);

void InterpolationMatrix1D(int _N, int NpointsIn, dfloat *rIn, int NpointsOut, dfloat *rOut, dfloat *I);
void DegreeRaiseMatrix1D(int Nc, int Nf, dfloat *P);
void CubatureWeakDmatrix1D(int _Nq, int _cubNq, dfloat *_cubProject, dfloat *_cubD, dfloat *_cubPDT);
dfloat JacobiP(dfloat a, dfloat alpha, dfloat beta, int _N);
dfloat GradJacobiP(dfloat a, dfloat alpha, dfloat beta, int _N);
void JacobiGLL(int _N, dfloat *_x, dfloat *_w = nullptr);
void JacobiGQ(dfloat alpha, dfloat beta, int _N, dfloat *_x, dfloat *_w);

#endif

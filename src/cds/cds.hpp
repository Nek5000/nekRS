#ifndef CDS_H
#define CDS_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "nrssys.hpp"
#include "mesh3D.h"
#include "elliptic.h"

#define NSCALAR_MAX 100

struct cds_t
{
  static constexpr double targetTimeBenchmark {0.1};
  int dim, elementType;

  mesh_t* mesh[NSCALAR_MAX];
  dlong fieldOffset[NSCALAR_MAX];
  dlong fieldOffsetScan[NSCALAR_MAX];
  dlong fieldOffsetSum;
  mesh_t* meshV;
  elliptic_t* solver[NSCALAR_MAX];

  int NVfields;            // Number of velocity fields
  int NSfields;            // Number of scalar fields

  setupAide options[NSCALAR_MAX];

  oogs_t *gsh, *gshT;

  dlong vFieldOffset;
  dlong vCubatureOffset;
  dfloat idt;
  dfloat *dt;
  int tstep;
  dfloat g0, ig0;

  int nEXT;
  int nBDF;
  int dtAdaptStep;

  int compute[NSCALAR_MAX];

  dfloat *U, *S;

  // filter
  int filterNc;
  dfloat filterS[NSCALAR_MAX];
  dfloat* filterM;
  occa::memory o_filterMT;

  //RK Subcycle Data
  int nRK;
  dfloat* coeffsfRK, * weightsRK, * nodesRK;
  occa::memory o_coeffsfRK, o_weightsRK;

  occa::memory o_relUrst;
  occa::memory o_Urst;

  //EXTBDF data
  dfloat* coeffEXT, * coeffBDF, * coeffSubEXT;

  int* EToB[NSCALAR_MAX];
  occa::memory o_EToB[NSCALAR_MAX];

  occa::memory* o_usrwrk;

  int Nsubsteps;

  dfloat* Ue;
  occa::memory o_Ue;

  dfloat* prop;
  occa::memory o_prop, o_ellipticCoeff;
  occa::memory o_rho, o_diff;

  dfloat* cU, * cSd, * cS, * FS, * BF;
  occa::memory o_cU, o_cSd, o_cS, o_FS, o_BF, o_BFDiag;

  occa::kernel sumMakefKernel;
  occa::kernel subCycleVolumeKernel,  subCycleCubatureVolumeKernel;
  occa::kernel subCycleSurfaceKernel, subCycleCubatureSurfaceKernel;
  occa::kernel subCycleRKUpdateKernel;
  occa::kernel subCycleStrongCubatureVolumeKernel;
  occa::kernel subCycleStrongVolumeKernel;
  occa::kernel subCycleRKKernel;
  occa::kernel subCycleInitU0Kernel;

  occa::kernel nStagesSum3Kernel;

  occa::kernel filterRTKernel; // Relaxation-Term based filtering
  // occa::kernel constrainKernel;

  occa::memory o_U;
  occa::memory o_S, o_Se;

  //EXTBDF data
  occa::memory o_coeffEXT, o_coeffBDF, o_coeffSubEXT;

  occa::kernel advectionVolumeKernel;
  occa::kernel advectionSurfaceKernel;
  occa::kernel advectionCubatureVolumeKernel;
  occa::kernel advectionCubatureSurfaceKernel;
  occa::kernel strongAdvectionVolumeKernel;
  occa::kernel strongAdvectionCubatureVolumeKernel;
  occa::kernel advectMeshVelocityKernel;

  occa::kernel helmholtzRhsIpdgBCKernel;
  occa::kernel helmholtzRhsBCKernel;
  occa::kernel dirichletBCKernel;
  occa::kernel setEllipticCoeffKernel;

  occa::kernel maskCopyKernel;

  occa::properties* kernelInfo;
};

occa::memory cdsSolve(int i, cds_t* cds, dfloat time, int stage);

#endif

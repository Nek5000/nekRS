#ifndef CDS_H
#define CDS_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "nrssys.hpp"
#include "mesh3D.h"
#include "elliptic.h"
#include "neknek.hpp"
#include "cvode.hpp"
#include <vector>

struct cds_t {
  static constexpr double targetTimeBenchmark{0.2};
  int dim, elementType;

  std::vector<mesh_t*> mesh;
  std::vector<dlong> fieldOffset;
  std::vector<dlong> fieldOffsetScan;
  occa::memory o_fieldOffsetScan;
  dlong fieldOffsetSum;
  mesh_t* meshV;
  std::vector<elliptic_t*> solver;
  neknek_t* neknek;
  cvode_t* cvode;

  bool anyCvodeSolver = false;
  bool anyEllipticSolver = false;

  int NVfields; // Number of velocity fields
  int NSfields; // Number of scalar fields

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

  std::vector<int> compute;
  std::vector<int> cvodeSolve;
  occa::memory o_compute;
  occa::memory o_cvodeSolve;

  dfloat *U, *S;

  // filter
  int filterNc;
  std::vector<dfloat> filterS;
  dfloat* filterM;
  occa::memory o_applyFilterRT;
  occa::memory o_filterS;
  occa::memory o_filterRT;
  int applyFilter;

  // RK Subcycle Data
  int nRK;
  dfloat *coeffsfRK, *weightsRK, *nodesRK;
  occa::memory o_coeffsfRK, o_weightsRK;

  occa::memory o_relUrst;
  occa::memory o_Urst;

  // EXTBDF data
  dfloat *coeffEXT, *coeffBDF;

  int *EToB;
  occa::memory o_EToB;
  dlong EToBOffset;

  occa::memory *o_usrwrk;

  int Nsubsteps;

  occa::memory o_Ue;

  occa::memory o_prop, o_ellipticCoeff;
  occa::memory o_rho, o_diff;

  dfloat *cU, *cSd, *cS;
  occa::memory o_cU, o_cSd, o_cS, o_FS, o_BF, o_BFDiag;

  occa::kernel sumMakefKernel;
  occa::kernel subCycleStrongCubatureVolumeKernel;
  occa::kernel subCycleStrongVolumeKernel;

  occa::kernel filterRTKernel;

  occa::memory o_U;
  occa::memory o_S, o_Se;

  occa::memory o_coeffEXT, o_coeffBDF;

  occa::kernel advectionVolumeKernel;
  occa::kernel advectionSurfaceKernel;
  occa::kernel advectionCubatureVolumeKernel;
  occa::kernel advectionCubatureSurfaceKernel;
  occa::kernel strongAdvectionVolumeKernel;
  occa::kernel strongAdvectionCubatureVolumeKernel;
  occa::kernel advectMeshVelocityKernel;

  occa::kernel neumannBCKernel;
  occa::kernel dirichletBCKernel;
  occa::kernel setEllipticCoeffKernel;

  occa::kernel maskCopyKernel;
  occa::kernel maskCopy2Kernel;

  occa::properties *kernelInfo;
};

occa::memory cdsSolve(int i, cds_t *cds, double time, int stage);

#endif

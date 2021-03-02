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
  int dim, elementType;

  mesh_t* meshT[NSCALAR_MAX];
  mesh_t* meshV;
  elliptic_t* solver[NSCALAR_MAX];

  int NVfields;            // Number of velocity fields
  int NSfields;            // Number of scalar fields

  setupAide options[NSCALAR_MAX];

  oogs_t *gsh, *gshT;

  dlong vFieldOffset;
  dlong fieldOffset;
  dlong Ntotal;
  dfloat idt;
  dfloat *dt;
  int tstep;
  dfloat g0, ig0;

  int nEXT;
  int nBDF;
  int dtAdaptStep;

  int compute[NSCALAR_MAX];

  dfloat* U, * S;
  dfloat* rkNS;
  //  dfloat *rhsS;
  dfloat* rkS;

  //RK Subcycle Data
  int SNrk;
  dfloat* Srka, * Srkb, * Srkc;
  occa::memory o_Srka, o_Srkb;

  //EXTBDF data
  dfloat* coeffEXT, * coeffBDF, * coeffSubEXT;
  dfloat* extC;

  int* mapB[NSCALAR_MAX], * EToB[NSCALAR_MAX];
  occa::memory o_mapB[NSCALAR_MAX];
  occa::memory o_EToB[NSCALAR_MAX];

  occa::memory* o_usrwrk;

  //halo data
  dfloat* sendBuffer;
  dfloat* recvBuffer;
  dfloat* haloGatherTmp;
  // //
  dfloat* ssendBuffer;
  dfloat* srecvBuffer;
  dfloat* shaloGatherTmp;

  occa::memory o_sendBuffer, h_sendBuffer;
  occa::memory o_recvBuffer, h_recvBuffer;
  occa::memory o_gatherTmpPinned, h_gatherTmpPinned;

  //
  occa::memory o_ssendBuffer, h_ssendBuffer;
  occa::memory o_srecvBuffer, h_srecvBuffer;
  occa::memory o_sgatherTmpPinned, h_sgatherTmpPinned;

  int Nsubsteps;
  dfloat sdt;
  dfloat* Ue;
  occa::memory o_Ue;

  int var_coeff;
  dfloat* prop, * ellipticCoeff;
  occa::memory o_prop, o_ellipticCoeff;
  occa::memory o_rho, o_diff;

  dfloat* cU, * cSd, * cS, * FS, * BF;
  occa::memory o_cU, o_cSd, o_cS, o_FS, o_BF, o_BFDiag;

  occa::memory o_wrk0, o_wrk1, o_wrk2, o_wrk3, o_wrk4, o_wrk5, o_wrk6;

  occa::kernel sumMakefKernel;
  occa::kernel subCycleVolumeKernel,  subCycleCubatureVolumeKernel;
  occa::kernel subCycleSurfaceKernel, subCycleCubatureSurfaceKernel;
  occa::kernel subCycleRKUpdateKernel;
  occa::kernel subCycleStrongCubatureVolumeKernel;
  occa::kernel subCycleStrongVolumeKernel;
  occa::kernel subCycleRKKernel;
  occa::kernel subCycleExtrapolateFieldKernel;
  occa::kernel subCycleExtrapolateScalarKernel;

  occa::kernel filterRTKernel; // Relaxation-Term based filtering
  // occa::kernel constrainKernel;

  occa::memory o_U;
  occa::memory o_S, o_Se;

  // occa::memory o_Vort, o_Div; // Not sure to keep it
  occa::memory o_haloBuffer;
  occa::memory o_haloGatherTmp;

  occa::memory o_shaloBuffer;
  occa::memory o_shaloGatherTmp;

  //ARK data
  occa::memory o_rkC;

  //EXTBDF data
  occa::memory o_coeffEXT, o_coeffBDF, o_coeffSubEXT;
  occa::memory o_extC;

// Will be depreceated.....AK
  occa::kernel haloExtractKernel;
  occa::kernel haloScatterKernel;
  occa::kernel scalarHaloExtractKernel;
  occa::kernel scalarHaloScatterKernel;

  occa::kernel haloGetKernel;
  occa::kernel haloPutKernel;
  occa::kernel scalarHaloGetKernel;
  occa::kernel scalarHaloPutKernel;

  occa::kernel setFlowFieldKernel;

  occa::kernel advectionVolumeKernel;
  occa::kernel advectionSurfaceKernel;
  occa::kernel advectionCubatureVolumeKernel;
  occa::kernel advectionCubatureSurfaceKernel;
  occa::kernel advectionStrongVolumeKernel;
  occa::kernel advectionStrongCubatureVolumeKernel;
  occa::kernel advectMeshVelocityKernel;

  occa::kernel helmholtzRhsIpdgBCKernel;
  occa::kernel helmholtzRhsBCKernel;
  occa::kernel dirichletBCKernel;
  occa::kernel setEllipticCoeffKernel;

  occa::kernel maskCopyKernel;

  occa::properties* kernelInfo;
};

occa::memory cdsSolve(int i, cds_t* cds, dfloat time);

#endif

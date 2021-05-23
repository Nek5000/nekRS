#include "elliptic.h"
#include "platform.hpp"
#include <string>
#include "gslib.h"
#include "fem_amg_preco.hpp"
#include "crs_hypre.h"
void ellipticSEMFEMSetup(elliptic_t* elliptic)
{
  mesh_t* mesh = elliptic->mesh;
  const int N = mesh->Nq;
  const int Nelements = mesh->Nelements;
  const int nullspace = (int) elliptic->allNeumann;
  double nullParams[] = {0}; // use default parameters
  double* mask = (double*) malloc(N*N*N*Nelements*sizeof(double));
  for(int i = 0; i < N*N*N*Nelements; ++i) mask[i] = 1.0;
  for(dlong n = 0; n < elliptic->Nmasked; n++){
    mask[elliptic->maskIds[n]] = 0.0;
  }
  
  SEMFEMData* data = fem_amg_setup(
    &N,
    &N,
    &N,
    &Nelements,
    mesh->x,
    mesh->y,
    mesh->z,
    mask,
    platform->comm.mpiComm,
    mesh->globalIds
  );

  const long long numRows = data->rowEnd - data->rowStart + 1;
  elliptic->o_dofMap = platform->device.malloc(numRows * sizeof(long long), data->dofMap);
  elliptic->o_SEMFEMBuffer1 = platform->device.malloc(elliptic->Nfields * elliptic->Ntotal,sizeof(dfloat));
  elliptic->o_SEMFEMBuffer2 = platform->device.malloc(elliptic->Nfields * elliptic->Ntotal,sizeof(dfloat));

  elliptic->numRowsSEMFEM = numRows;


  std::string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string oklpath = install_dir + "/okl/elliptic/";
  std::string filename = oklpath + "ellipticPreSEMFEM.okl";
  elliptic->preSEMFEMKernel = platform->device.buildKernel(
    filename,
    "preSEMFEM",
    platform->kernelInfo
  );
  filename = oklpath + "ellipticPostSEMFEM.okl";
  elliptic->postSEMFEMKernel = platform->device.buildKernel(
    filename,
    "postSEMFEM",
    platform->kernelInfo
  );

  elliptic->hypreData =
    hypre_setup(
      numRows,
      data->rowStart,
      data->nnz,
      data->Ai,
      data->Aj,
      data->Av,
      nullspace,
      platform->comm.mpiComm,
      1,
      &nullParams[0]
    );

  delete data;
}
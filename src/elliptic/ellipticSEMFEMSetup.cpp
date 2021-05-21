#include "elliptic.h"
#include "platform.hpp"
#include <string>
#include "gslib.h"
extern "C" long long * get_dof_map();
extern "C" long long * get_row_start();
extern "C" long long * get_row_end();
extern "C" void fem_amg_setup(const sint *n_x_, const sint *n_y_, const sint *n_z_, 
                   const sint *n_elem_, const sint *n_dim_, 
                   double *x_m_, double *y_m_, double *z_m_, 
                   double *pmask_, const sint *nullspace,
                   double *param,
                   MPI_Comm comm,
                   long long int* gatherGlobalNodes
                   );
void ellipticSEMFEMSetup(elliptic_t* elliptic)
{
  mesh_t* mesh = elliptic->mesh;
  const int N = mesh->Nq;
  const int Nelements = mesh->Nelements;
  const int Ndim = 3;
  const int nullspace = (int) elliptic->allNeumann;
  double nullParams[] = {0}; // use default parameters
  double* mask = (double*) malloc(N*N*N*Nelements*sizeof(double));
  for(int i = 0; i < N*N*N*Nelements; ++i) mask[i] = 1.0;
  for(dlong n = 0; n < elliptic->Nmasked; n++){
    mask[elliptic->maskIds[n]] = 0.0;
  }
  
  fem_amg_setup(
    &N,
    &N,
    &N,
    &Nelements,
    &Ndim,
    mesh->x,
    mesh->y,
    mesh->z,
    mask,
    &nullspace,
    &nullParams[0],
    platform->comm.mpiComm,
    mesh->globalIds
  );

  long long row_start = *get_row_start();
  long long row_end = *get_row_end();
  const long long numRows = row_end - row_start + 1;
  long long * h_dofMap = get_dof_map();
  elliptic->o_dofMap = platform->device.malloc(numRows * sizeof(long long), h_dofMap);
  elliptic->o_SEMFEMBuffer1 = platform->device.malloc(elliptic->Nfields * elliptic->Ntotal,sizeof(dfloat));
  elliptic->o_SEMFEMBuffer2 = platform->device.malloc(elliptic->Nfields * elliptic->Ntotal,sizeof(dfloat));


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
}
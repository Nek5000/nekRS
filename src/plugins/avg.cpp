/*
     compute running averages E(X), E(X*X) 
     and E(X*Y) for veloctiy only

     statistics can be obtained by:

     avg(X)   := E(X)
     var(X)   := E(X*X) - E(X)*E(X)
     cov(X,Y) := E(X*Y) - E(X)*E(Y)

     Note: The E-operator is linear, in the sense that the expected
           value is given by E(X) = 1/N * avg[ E(X)_i ], where E(X)_i
           is the expected value of the sub-ensemble i (i=1...N).
*/

#include <nekrs.hpp>
#include <nekInterfaceAdapter.hpp>
#include "avg.hpp"

// private members
namespace {
  static ogs_t *ogs;
  static ins_t *ins;

  static occa::memory o_Uavg;
  static occa::memory o_Uav2;
  static occa::memory o_Pavg;
  static occa::memory o_Savg;
  
  static occa::kernel avgX_XXKernel;
  static occa::kernel avgXYKernel;

  static bool buildKernelCalled = 0;
  static bool setupCalled = 0;

  static dfloat atime = 0;
  static dfloat timel;
}

void avg::buildKernel(ins_t *ins)
{
  mesh_t *mesh = ins->mesh; 

  string fileName;
  int rank = mesh->rank;
  fileName.assign(getenv("NEKRS_INSTALL_DIR"));
  fileName += "/okl/avg/avg.okl";
  occa::properties& kernelInfo = *ins->kernelInfo;
  for (int r=0;r<2;r++){
    if ((r==0 && rank==0) || (r==1 && rank>0)) {
       avgX_XXKernel     = mesh->device.buildKernel(fileName.c_str(), "avgX_XX"    , kernelInfo);
       avgXYKernel = mesh->device.buildKernel(fileName.c_str(), "avgXY", kernelInfo);
    }
    MPI_Barrier(mesh->comm);
  }
} 

void avg::run(dfloat time)
{
  mesh_t *mesh = ins->mesh;

  const dfloat dtime = time - timel;
  atime += dtime; 

  if(atime == 0 || dtime == 0) return;
 
  const dfloat b = dtime/atime;
  const dfloat a = 1-b;
  const dlong N  = ins->fieldOffset; 

  avgX_XXKernel(N, ins->fieldOffset, ins->fieldOffset*ins->NVfields, 
                ins->NVfields, a, b, ins->o_U, o_Uavg);
  avgX_XXKernel(N, ins->fieldOffset, ins->fieldOffset, 
                1, a, b, ins->o_P, o_Pavg);

  if(ins->Nscalar) {
    cds_t *cds = ins->cds; 
    avgX_XXKernel(N, ins->fieldOffset, cds->fieldOffset*cds->NSfields , 
                  cds->NSfields, a, b, cds->o_S, o_Savg);
  }

  avgXYKernel(N, ins->fieldOffset, a, b, ins->o_U, o_Uav2);

  timel = time;
}

void avg::setup(ins_t *ins_)
{
  ins = ins_;
  mesh_t *mesh = ins->mesh;

  if(setupCalled) return;

  o_Uavg = mesh->device.malloc(2*ins->fieldOffset*ins->NVfields*sizeof(dfloat));
  ins->setScalarKernel(2*ins->fieldOffset*ins->NVfields, 0.0, o_Uavg);

  o_Uav2 = mesh->device.malloc(ins->fieldOffset*ins->NVfields*sizeof(dfloat));
  ins->setScalarKernel(ins->fieldOffset*ins->NVfields, 0.0, o_Uav2);

  o_Pavg = mesh->device.malloc(2*ins->fieldOffset*sizeof(dfloat));
  ins->setScalarKernel(2*ins->fieldOffset, 0.0, o_Pavg);
 
  if(ins->Nscalar) {
    cds_t *cds = ins->cds; 
    o_Savg = mesh->device.malloc(2*cds->fieldOffset*cds->NSfields*sizeof(dfloat));
    ins->setScalarKernel(2*cds->fieldOffset*cds->NSfields, 0.0, o_Savg);
  }

  setupCalled = 1;
  timel = ins->startTime;
}

void avg::outfld()
{
  cds_t *cds = ins->cds; 
  mesh_t *mesh = ins->mesh;
  const int FP64 = 1;

  occa::memory o_null;
  occa::memory o_S; 

  const int Nscalar = ins->Nscalar;
  if(ins->Nscalar)
    occa::memory o_S = o_Savg + cds->NSfields*cds->fieldOffset*sizeof(dfloat); 

  nek_outfld("avg", atime, o_null, 
             o_Uavg, 
             o_Pavg, 
             o_S, 
             Nscalar, 
             FP64); 

  nek_outfld("rms", atime, o_null, 
             o_Uavg + ins->NVfields*ins->fieldOffset*sizeof(dfloat), 
             o_Pavg + ins->fieldOffset*sizeof(dfloat), 
             o_S, 
             Nscalar, 
             FP64); 

  nek_outfld("rm2", atime, o_null, 
             o_Uav2, 
             o_null, 
             o_null,
             0, 
             FP64); 

  atime = 0;
}

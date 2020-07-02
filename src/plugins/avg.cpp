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

  static occa::memory o_Uavg, o_Urms;
  static occa::memory o_Urm2;
  static occa::memory o_Pavg, o_Prms;
  static occa::memory o_Savg, o_Srms;
  
  static occa::kernel EXKernel;
  static occa::kernel EXXKernel;
  static occa::kernel EXYKernel;

  static bool buildKernelCalled = 0;
  static bool setupCalled = 0;

  static int counter = 0;

  static dfloat atime;
  static dfloat timel;
}

void avg::buildKernel(ins_t *ins)
{
  mesh_t *mesh = ins->mesh; 

  string fileName;
  int rank = mesh->rank;
  fileName.assign(getenv("NEKRS_INSTALL_DIR"));
  fileName += "/okl/plugins/avg.okl";
  occa::properties& kernelInfo = *ins->kernelInfo;
  for (int r=0;r<2;r++){
    if ((r==0 && rank==0) || (r==1 && rank>0)) {
       EXKernel  = mesh->device.buildKernel(fileName.c_str(), "EX"  , kernelInfo);
       EXXKernel = mesh->device.buildKernel(fileName.c_str(), "EXX" , kernelInfo);
       EXYKernel = mesh->device.buildKernel(fileName.c_str(), "EXY" , kernelInfo);
    }
    MPI_Barrier(mesh->comm);
  }
  buildKernelCalled = 1;
} 

void avg::reset()
{
  counter = 0;
  atime   = 0;
}

void avg::EX (dlong N, dfloat a, dfloat b, int nflds, occa::memory o_x, occa::memory o_EX)
{
  EXKernel(N, ins->fieldOffset, nflds, a, b, o_x, o_EX); 
}

void avg::EXX(dlong N, dfloat a, dfloat b, int nflds, occa::memory o_x, occa::memory o_EXX)
{
  EXXKernel(N, ins->fieldOffset, nflds, a, b, o_x, o_EXX); 
}

void avg::EXY(dlong N, dfloat a, dfloat b, int nflds, occa::memory o_x, occa::memory o_y, occa::memory o_EXY)
{
  EXYKernel(N, ins->fieldOffset, nflds, a, b, o_x, o_y, o_EXY); 
}

void avg::run(dfloat time)
{
  if(!setupCalled || !buildKernelCalled) { 
    cout << "avg::run() was called prior to avg::setup()!\n";
    ABORT(1);
  }

  if(!counter) {
    atime = 0;
    timel = time;
  }
  counter++;

  const dfloat dtime = time - timel;
  atime += dtime; 

  if(atime == 0 || dtime == 0) return;
 
  const dfloat b = dtime/atime;
  const dfloat a = 1-b;

  mesh_t *mesh = ins->mesh;
  const dlong N = mesh->Nelements*mesh->Np; 

  // velocity
  EX (N, a, b, ins->NVfields, ins->o_U, o_Uavg);
  EXX(N, a, b, ins->NVfields, ins->o_U, o_Urms);

  const dlong offsetByte = ins->fieldOffset*sizeof(dfloat);
  occa::memory o_vx = ins->o_U + 0*offsetByte;
  occa::memory o_vy = ins->o_U + 1*offsetByte; 
  occa::memory o_vz = ins->o_U + 2*offsetByte; 

  EXY(N, a, b, 1, o_vx, o_vy, o_Urm2 + 0*offsetByte);
  EXY(N, a, b, 1, o_vy, o_vz, o_Urm2 + 1*offsetByte);
  EXY(N, a, b, 1, o_vz, o_vx, o_Urm2 + 2*offsetByte);

  // pressure
  EX (N, a, b, 1, ins->o_P, o_Pavg);
  EXX(N, a, b, 1, ins->o_P, o_Prms);

  // scalars
  if(ins->Nscalar) {
    cds_t *cds = ins->cds; 
    const dlong N = cds->mesh->Nelements*cds->mesh->Np; 
    EX (N, a, b, cds->NSfields, cds->o_S, o_Savg);
    EXX(N, a, b, cds->NSfields, cds->o_S, o_Srms);
  }

  timel = time;
}

void avg::setup(ins_t *ins_)
{
  if(!buildKernelCalled) { 
    cout << "avg::setup() was called prior avg::buildKernel()!\n";
    ABORT(1);
  }

  ins = ins_;
  mesh_t *mesh = ins->mesh;

  if(setupCalled) return;

  o_Uavg = mesh->device.malloc(ins->fieldOffset*ins->NVfields*sizeof(dfloat));
  o_Urms = mesh->device.malloc(ins->fieldOffset*ins->NVfields*sizeof(dfloat));
  ins->setScalarKernel(ins->fieldOffset*ins->NVfields, 0.0, o_Uavg);
  ins->setScalarKernel(ins->fieldOffset*ins->NVfields, 0.0, o_Urms);

  o_Urm2 = mesh->device.malloc(ins->fieldOffset*ins->NVfields*sizeof(dfloat));
  ins->setScalarKernel(ins->fieldOffset*ins->NVfields, 0.0, o_Urm2);

  o_Pavg = mesh->device.malloc(ins->fieldOffset*sizeof(dfloat));
  o_Prms = mesh->device.malloc(ins->fieldOffset*sizeof(dfloat));
  ins->setScalarKernel(ins->fieldOffset, 0.0, o_Pavg);
  ins->setScalarKernel(ins->fieldOffset, 0.0, o_Prms);
 
  if(ins->Nscalar) {
    cds_t *cds = ins->cds; 
    o_Savg = mesh->device.malloc(cds->fieldOffset*cds->NSfields*sizeof(dfloat));
    o_Srms = mesh->device.malloc(cds->fieldOffset*cds->NSfields*sizeof(dfloat));
    ins->setScalarKernel(cds->fieldOffset*cds->NSfields, 0.0, o_Savg);
    ins->setScalarKernel(cds->fieldOffset*cds->NSfields, 0.0, o_Srms);
  }

  setupCalled = 1;
}

void avg::outfld()
{
  cds_t *cds = ins->cds; 
  mesh_t *mesh = ins->mesh;
  const int FP64 = 1;
  const int coords = 0;

  occa::memory o_null;
  occa::memory o_Tavg, o_Trms;

  const int Nscalar = ins->Nscalar;
  if(ins->Nscalar) {
    occa::memory o_Tavg  = o_Savg; 
    occa::memory o_Trms = o_Srms; 
  }

  nek_outfld("avg", atime, coords, 
             o_Uavg, 
             o_Pavg, 
             o_Tavg, 
             Nscalar, 
             FP64); 

  nek_outfld("rms", atime, coords, 
             o_Urms, 
             o_Prms, 
             o_Trms, 
             Nscalar, 
             FP64); 

  nek_outfld("rm2", atime, coords, 
             o_Urm2, 
             o_null, 
             o_null,
             0, 
             FP64); 

  atime = 0;
}

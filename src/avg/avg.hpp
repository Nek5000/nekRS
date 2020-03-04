/*
     compute averages E(X), E(X*X) and E(X*Y)

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

// private members
namespace {
  static ogs_t *ogs;
  static ins_t *ins;

  static dfloat *Uavg;
  static dfloat *Pavg;
  static dfloat *Savg;
 
  static occa::memory o_Uavg;
  static occa::memory o_Pavg;
  static occa::memory o_Savg;
  
  static occa::kernel avgX_XXKernel;
  static occa::kernel avgXY_YZ_XZKernel;

  static bool buildKernelCalled = 0;
  static bool setupCalled = 0;

  static dfloat atime = 0;
  static dfloat timel;
}

namespace avg {

void buildKernel(ins_t *ins)
{
  mesh_t *mesh = ins->mesh; 

  string fileName;
  int rank = mesh->rank;
  fileName.assign(getenv("NEKRS_INSTALL_DIR"));
  fileName += "/plugin/okl/avg.okl";
  occa::properties& kernelInfo = *ins->kernelInfo;
  for (int r=0;r<2;r++){
    if ((r==0 && rank==0) || (r==1 && rank>0)) {
       avgX_XXKernel     = mesh->device.buildKernel(fileName.c_str(), "avgX_XX"    , kernelInfo);
       avgXY_YZ_XZKernel = mesh->device.buildKernel(fileName.c_str(), "avgXY_YZ_XZ", kernelInfo);
    }
    MPI_Barrier(mesh->comm);
  }
} 

void run(dfloat time)
{
  mesh_t *mesh = ins->mesh;
  cds_t *cds = ins->cds; 

  const dfloat dtime = time - timel;
  atime += dtime; 

  if(atime == 0 || dtime == 0) return;
 
  const dfloat b = dtime/atime;
  const dfloat a = 1-b;

  const dlong N = ins->fieldOffset; 

  avgX_XXKernel(N, ins->fieldOffset, ins->NVfields, a, b, ins->o_U, o_Uavg);
  const dlong UavgFieldOffset = 2*ins->fieldOffset*ins->NVfields;
  avgXY_YZ_XZKernel(N, ins->fieldOffset, UavgFieldOffset, a, b, ins->o_U, o_Uavg);

  avgX_XXKernel(N, ins->fieldOffset, 1, a, b, ins->o_P, o_Pavg);
  avgX_XXKernel(N, ins->fieldOffset, cds->NSfields, a, b, cds->o_S, o_Savg);

  timel = time;
}


void setup(ins_t *ins_)
{
  ins = ins_;
  cds_t *cds = ins->cds; 
  mesh_t *mesh = ins->mesh;

  if(setupCalled) return;

  const dlong NtotalUavg = 2*ins->fieldOffset*ins->NVfields + 3*ins->fieldOffset;
  Uavg = (dfloat *)calloc(NtotalUavg, sizeof(dfloat));
  o_Uavg = mesh->device.malloc(NtotalUavg*sizeof(dfloat), Uavg);

  Pavg = (dfloat *)calloc(2*ins->fieldOffset, sizeof(dfloat));
  o_Pavg = mesh->device.malloc(2*ins->fieldOffset*sizeof(dfloat), Pavg);

  Savg = (dfloat *)calloc(2*cds->fieldOffset*cds->NSfields, sizeof(dfloat));
  o_Savg = mesh->device.malloc(2*cds->fieldOffset*cds->NSfields*sizeof(dfloat), Savg);

  setupCalled = 1;
  timel = ins->startTime;
}

void outfld()
{
  cds_t *cds = ins->cds; 
  mesh_t *mesh = ins->mesh;

  o_Uavg.copyTo(Uavg);
  o_Pavg.copyTo(Pavg);
  if(cds->NSfields) o_Savg.copyTo(Savg);

  *(nekData.time) = atime;
  const int p62 = nekData.param[62];
  nekData.param[62] = 1; // enforce in DP

  for(int fld=0; fld<2; fld++) {
    dlong Nlocal = mesh->Nelements * mesh->Np;
 
    dfloat *vx = Uavg + ins->fieldOffset*(2*ins->NVfields + 3 + 0);
    dfloat *vy = Uavg + ins->fieldOffset*(2*ins->NVfields + 3 + 1);
    dfloat *vz = Uavg + ins->fieldOffset*(2*ins->NVfields + 3 + 2);
    dfloat *p  = Pavg + fld*ins->fieldOffset; 
    memcpy(nekData.vx, vx, sizeof(dfloat)*Nlocal);
    memcpy(nekData.vy, vy, sizeof(dfloat)*Nlocal);
    memcpy(nekData.vz, vz, sizeof(dfloat)*Nlocal);
    memcpy(nekData.pr, p , sizeof(dfloat)*Nlocal);
    if(cds->NSfields) {
      const dlong nekFieldOffset = nekData.lelt*mesh->Np;
      for(int is=0; is<ins->Nscalar; is++) {
        mesh_t *mesh;
        (is) ? mesh = cds->meshV : mesh = cds->mesh;
        const dlong Nlocal = mesh->Nelements * mesh->Np;
        dfloat *Ti = nekData.t   + is*nekFieldOffset;
        dfloat *Si = Savg + fld*cds->NSfields*cds->fieldOffset + is*cds->fieldOffset;
        memcpy(Ti, Si, Nlocal*sizeof(dfloat));
      }
    } 
    if(fld == 0) nek_outfld("avg"); 
    if(fld == 1) nek_outfld("rms");
  }

  nekData.param[62] = p62;
  atime = 0;
}

} // namespace 

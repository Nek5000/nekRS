/*
   copy velocity data of a given slab (slabIdSrc) to another slab
   (slideIdDst) also known as recycling

   Note: This implementation relies on a special global element 
         numbering which is only true for extruded meshes in z from nek!
*/

#include <nekInterfaceAdapter.hpp>

namespace velRecycling {

static ogs_t *ogs;
static ins_t *ins;

static occa::memory o_wrk;

static dfloat *flux, *area;
static occa::memory o_flux, o_area;

static dfloat *tmp1, *tmp2;
static occa::memory o_tmp1, o_tmp2;

static occa::kernel setValueBCKernel;
static occa::kernel getBCFluxKernel;
static occa::kernel sumReductionKernel; 
static occa::kernel scalarMultiplyKernel; 

bool buildKernelCalled = 0;
bool setupCalled = 0;

void buildKernel(ins_t *ins)
{
  mesh_t *mesh = ins->mesh; 

  string fileName;
  int rank = mesh->rank;
  fileName.assign(getenv("NEKRS_INSTALL_DIR"));
  fileName += "/okl/recycling.okl";
  occa::properties& kernelInfo = *ins->kernelInfo;
  for (int r=0;r<2;r++){
    if ((r==0 && rank==0) || (r==1 && rank>0)) {
       setValueBCKernel     =  mesh->device.buildKernel(fileName.c_str(), "setBCVectorValue", kernelInfo);
       getBCFluxKernel      =  mesh->device.buildKernel(fileName.c_str(), "getBCFlux", kernelInfo);
       sumReductionKernel   =  mesh->device.buildKernel(fileName.c_str(), "sumReduction", kernelInfo);
       scalarMultiplyKernel =  mesh->device.buildKernel(fileName.c_str(), "scalarMultiply", kernelInfo);
    }
    MPI_Barrier(mesh->comm);
  }
} 

void copy()
{
  mesh_t *mesh = ins->mesh; 
  const int bc = 2;
  const dfloat wbar = 1.0; 

  // copy recycling plane in interior to inlet
  o_wrk.copyFrom(ins->o_U, ins->NVfields*ins->Ntotal*sizeof(dfloat));
  setValueBCKernel(mesh->Nelements, 0.0, bc, ins->fieldOffset,
                   o_wrk, mesh->o_vmapM, ins->o_EToB);

  //ogsGatherScatterMany(o_wrk, ins->NVfields, ins->fieldOffset,
  //                     ogsDfloat, ogsAdd, ogs);
  for(int k=0;k<ins->dim;++k)
    ogsGatherScatter(o_wrk+k*ins->fieldOffset*sizeof(dfloat),
                     ogsDfloat, ogsAdd, ogs);
  
  // rescale
  getBCFluxKernel(mesh->Nelements, bc, ins->fieldOffset, o_wrk,
                  mesh->o_vmapM, ins->o_EToB, mesh->o_sgeo, o_area, o_flux);

  const int NfpTotal = mesh->Nelements*mesh->Nfaces*mesh->Nfp; 
  sumReductionKernel(NfpTotal, o_area, o_flux, o_tmp1, o_tmp2); 

  o_tmp1.copyTo(tmp1);
  o_tmp2.copyTo(tmp2);
  dfloat sbuf[2] = {0,0};
  for(int n=0; n<blockSize; n++){
    sbuf[0] += tmp1[n]; 
    sbuf[1] += tmp2[n]; 
  }
  MPI_Allreduce(MPI_IN_PLACE, sbuf, 2, MPI_DFLOAT, MPI_SUM, mesh->comm); 

  const dfloat scale = -wbar*sbuf[0] / sbuf[1]; 
  //printf("here: %f %f %f\n", scale, sbuf[0], sbuf[1]);
  scalarMultiplyKernel(ins->NVfields*ins->Ntotal, scale, o_wrk);
}


void setup(ins_t *ins_, occa::memory o_wrk_, const hlong eOffset)
{
  ins = ins_;
  o_wrk = o_wrk_; 
  mesh_t *mesh = ins->mesh;

  const int Ntotal = mesh->Np * mesh->Nelements;
  hlong *ids = (hlong *) calloc(Ntotal, sizeof(hlong));

  for (int e=0; e < mesh->Nelements; e++){
    // establish a unique numbering
    const int eg = nek_lglel(e); // 0-based
    for (int n=0; n < mesh->Np; n++)  {
      ids[e*mesh->Np + n] = eg*mesh->Np + n+1; 
    }
    
    for (int n=0; n < mesh->Nfp*mesh->Nfaces; n++) {
      const int f = n/mesh->Nfp;
      const int idM = ins->mesh->vmapM[e*mesh->Nfp*mesh->Nfaces + n];
      const int id  = ins->EToB[f + e*mesh->Nfaces];
      if (id == 2) ids[idM] += eOffset*mesh->Np; 
    }
   
  }

  ogs = ogsSetup(Ntotal, ids, mesh->comm, 1, mesh->device);
  free(ids);

  const int NfpTotal = mesh->Nelements*mesh->Nfaces*mesh->Nfp;
  const int ntmp = (NfpTotal+blockSize-1)/blockSize;
  tmp1   = (dfloat *) calloc(ntmp, sizeof(dfloat));
  tmp2   = (dfloat *) calloc(ntmp, sizeof(dfloat));

  o_tmp1 = mesh->device.malloc(ntmp*sizeof(dfloat), tmp1);
  o_tmp2 = mesh->device.malloc(ntmp*sizeof(dfloat), tmp2);

  flux   = (dfloat *)calloc(NfpTotal, sizeof(dfloat));
  area   = (dfloat *)calloc(NfpTotal, sizeof(dfloat));
  o_flux = mesh->device.malloc(NfpTotal*sizeof(dfloat), flux);
  o_area = mesh->device.malloc(NfpTotal*sizeof(dfloat), area);
}


} // namespace 

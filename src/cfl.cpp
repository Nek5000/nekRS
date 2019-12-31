#include "nrs.hpp"

static int firstTime = 1;
static dfloat *tmp;
static occa::memory o_tmp;

void setup(ins_t *ins){

  mesh_t *mesh = ins->mesh; 
  
  dfloat *dH; 
  if(ins->elementType==QUADRILATERALS || ins->elementType==HEXAHEDRA){
    dH = (dfloat*) calloc((mesh->N+1),sizeof(dfloat));

    for(int n=0; n<(mesh->N+1); n++){
      if(n==0)
        dH[n] = mesh->gllz[n+1] - mesh->gllz[n];
      else if(n==mesh->N)
        dH[n] = mesh->gllz[n] - mesh->gllz[n-1];
      else
        dH[n] = 0.5*( mesh->gllz[n+1] - mesh->gllz[n-1]); 
    }
    for(int n=0; n< (mesh->N+1); n++)
      dH[n] = 1.0/dH[n]; 

    ins->o_idH = mesh->device.malloc((mesh->N+1)*sizeof(dfloat), dH); 
    free(dH); 
  }

  tmp = (dfloat *) calloc(ins->Nblock, sizeof(dfloat));
  o_tmp = mesh->device.malloc(ins->Nblock*sizeof(dfloat), tmp);

  firstTime = 0; 

}

dfloat computeCFL(ins_t *ins, dfloat time, int tstep){

  mesh_t *mesh = ins->mesh; 
  if(firstTime) setup(ins);

  // Compute cfl factors i.e. dt* U / h 
  ins->cflKernel(mesh->Nelements,
                 ins->dt, 
                 mesh->o_vgeo,
                 ins->o_idH,
                 ins->fieldOffset,
                 ins->o_U,
                 ins->o_wrk0);  
  
  // find the local maximum of CFL number
  ins->maxKernel(ins->Nlocal, ins->o_wrk0, o_tmp);
  o_tmp.copyTo(tmp);
  
  // finish reduction
  dfloat cfl = 0.f; 
  for(dlong n=0;n<ins->Nblock;++n){
    cfl  = mymax(cfl, tmp[n]);
  }

  dfloat gcfl = 0.f;
  MPI_Allreduce(&cfl, &gcfl, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
  
  return gcfl; 
}


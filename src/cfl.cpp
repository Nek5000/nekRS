#include "nekrs.hpp"

using namespace nekrs::mpi;

void getDh(ins_t *ins){

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

  ins->computedDh = 1; 

}

dfloat computeCFL(ins_t *ins, dfloat time, int tstep){

  mesh_t *mesh = ins->mesh; 
  // create dH i.e. nodal spacing in reference coordinates
  if(!ins->computedDh)
    getDh(ins);
  // Compute cfl factors i.e. dt* U / h 
  ins->cflKernel(mesh->Nelements,
                 ins->dt, 
                 mesh->o_vgeo,
                 ins->o_idH,
                 ins->fieldOffset,
                 ins->o_U,
                 ins->o_rhsU);  
  
  // find the local maximum of CFL number
  const dlong Ntotal = mesh->Np*mesh->Nelements; 
  dfloat *tmp        = (dfloat *) calloc(ins->Nblock, sizeof(dfloat));
  occa::memory o_tmp = mesh->device.malloc(ins->Nblock*sizeof(dfloat), tmp);
  ins->maxKernel(Ntotal, ins->o_rhsU, o_tmp);
  o_tmp.copyTo(tmp);
  
  // finish reduction
  dfloat cfl = 0.f; 
  for(dlong n=0;n<ins->Nblock;++n){
    cfl  = mymax(cfl, tmp[n]);
  }

  dfloat gcfl = 0.f;
  Allreduce(&cfl, &gcfl, sizeof(dfloat), MPI_MAX);
  
  free(tmp);
  return gcfl; 
}


#include "nrs.hpp"
#include "platform.hpp"
#include "linAlg.hpp"
#include "postProcessing.hpp"

namespace {
  dfloat *drag;
  occa::memory o_drag;
}

dfloat postProcessing::viscousDrag(nrs_t *nrs, int nbID, const occa::memory& o_bID, occa::memory& o_Sij)
{
  mesh_t *mesh = nrs->meshV;

  if(o_drag.size() == 0) { 
    drag = (dfloat *) calloc(mesh->Nelements, sizeof(dfloat));
    o_drag = platform->device.malloc(mesh->Nelements * sizeof(dfloat), drag);
  } 

  auto dragKernel = platform->kernels.get("drag");
  
  dragKernel(mesh->Nelements,
	     nrs->fieldOffset,
	     nbID,
	     o_bID,
	     mesh->o_sgeo,
	     mesh->o_vmapM,
	     mesh->o_EToB,
	     nrs->o_mue,
	     o_Sij,
	     o_drag);

  o_drag.copyTo(drag);

  dfloat sum = 0;
  for (dlong i = 0; i < mesh->Nelements; i++) {
    sum += drag[i];
  }

  MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);
    
  return sum;
}

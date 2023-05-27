#include "nrs.hpp"
#include "platform.hpp"
#include "linAlg.hpp"
#include "postProcessing.hpp"

void postProcessing::strainRotationRate(nrs_t *nrs, bool smooth, bool rotationRate, occa::memory& o_SO)
{
  mesh_t *mesh = nrs->meshV;

  const int nFields = (rotationRate) ? 2*nrs->NVfields + nrs->NVfields : 2*nrs->NVfields;

  nrsCheck(o_SO.size() < (nFields*nrs->fieldOffset*sizeof(dfloat)),
           MPI_COMM_SELF, EXIT_FAILURE, "o_SO too small to store %d fields!\n", nFields); 
 
  nrs->SijOijKernel(mesh->Nelements,
		    nrs->fieldOffset,
		    (int) rotationRate,
                    (int) smooth,
		    mesh->o_vgeo,
		    mesh->o_D,
		    nrs->o_U,
		    o_SO);

  if(smooth) { 
    oogs::startFinish(o_SO, 
                      nFields,
                      nrs->fieldOffset,
                      ogsDfloat,
                      ogsAdd,
                      nrs->gsh);

    platform->linAlg->axmyMany(mesh->Nlocal, 
                               nFields, 
                               nrs->fieldOffset, 
                               0, 
                               1.0, 
                               mesh->o_invLMM, 
                               o_SO);
  } 
}

void postProcessing::strainRate(nrs_t *nrs, bool smooth, occa::memory& o_S)
{
  strainRotationRate(nrs, smooth, false, o_S); 
}

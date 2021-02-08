#include <mesh.h>
#include <linAlg.hpp>
#include <nekInterfaceAdapter.hpp>
void mesh_t::computeInvLMM()
{
  o_invLMM.copyFrom(o_LMM, Nelements * Np * sizeof(dfloat));
  oogs::startFinish(o_invLMM, 1, 0, ogsDfloat, ogsAdd, oogs);
  linAlg->ady(Nelements * Np, 1.0, o_invLMM);
}
void mesh_t::computeBdivW()
{
  strongDivergenceKernel(
    Nelements,
    o_vgeo,
    o_Dmatrices,
    fieldOffset,
    o_U,
    o_BdivW
  );
}
void mesh_t::move(){
  // update o_x, o_y and o_z based on mesh->o_U using AB formula
  nStagesSumVectorKernel(
      Nelements * Np,
      fieldOffset,
      o_ABCoeff,
      o_U,
      o_x,
      o_y,
      o_z
  );
  update();
}
void mesh_t::update(){
    geometricFactorsKernel(
        Nelements,
        o_D,
        o_gllw,
        o_x,
        o_y,
        o_z,
        o_cubInterpT,
        o_cubw,
        o_LMM,
        o_vgeo,
        o_ggeo,
        o_cubvgeo
    );
    volume = linAlg->sum(Nelements * Np, o_LMM, comm);
    computeInvLMM();
    surfaceGeometricFactorsKernel(
        Nelements,
        o_D,
        o_gllw,
        o_faceNodes,
        o_x,
        o_y,
        o_z,
        o_sgeo
    );
}

#include <mesh.h>
#include <linAlg.hpp>
#include <nekInterfaceAdapter.hpp>
void mesh_t::computeInvMassMatrix()
{
    o_scratch.copyFrom(o_LMM, Nelements * Np * sizeof(dfloat));
    ogsGatherScatter(o_scratch, ogsDfloat, ogsAdd, ogs);
    // invLMM[n] = 1.0 / LMM[n]
    scalarDivideKernel(Nelements * Np, 1.0, o_scratch);
    o_invLMM.copyFrom(o_scratch, Nelements * Np * sizeof(dfloat));
}
void mesh_t::move(int tstep){

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
    volume = linAlg_t::getInstance()->sum(Nelements * Np, o_LMM, comm);
    computeInvMassMatrix();
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

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
void mesh_t::move(){

  linAlg_t* linAlg = linAlg_t::getInstance();
  // update o_x, o_y and o_z based on mesh->o_U using AB formula
  occa::memory o_x_increment = o_scratch + 0 * fieldOffset * sizeof(dfloat);
  occa::memory o_y_increment = o_scratch + 1 * fieldOffset * sizeof(dfloat);
  occa::memory o_z_increment = o_scratch + 2 * fieldOffset * sizeof(dfloat);
  nStagesSumVectorKernel(
      Nelements * Np,
      fieldOffset,
      o_ABCoeff,
      o_U,
      o_x_increment,
      o_y_increment,
      o_z_increment
  );

  // for testing
  std::cout << "sum(o_x_increment) = " << linAlg->sum(Nelements*Np, o_x_increment, comm) << "\n";
  std::cout << "sum(o_y_increment) = " << linAlg->sum(Nelements*Np, o_y_increment, comm) << "\n";
  std::cout << "sum(o_z_increment) = " << linAlg->sum(Nelements*Np, o_z_increment, comm) << "\n";


  // o_x += o_x_increment
  linAlg->axpby(Nelements*Np, 1.0, o_x_increment, 1.0, o_x);
  linAlg->axpby(Nelements*Np, 1.0, o_y_increment, 1.0, o_y);
  linAlg->axpby(Nelements*Np, 1.0, o_z_increment, 1.0, o_z);
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

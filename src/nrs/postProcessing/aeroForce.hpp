#ifndef AERO_FORCE_HPP_
#define AERO_FORCE_HPP_ 

#include "platform.hpp"

class AeroForce
{
public:
  AeroForce() {};
  AeroForce(std::tuple< std::array<dfloat, 3>, std::array<dfloat, 3> > f) 
  {
    forceV = std::get<0>(f);
    forceP = std::get<1>(f); 
  };

  std::array<dfloat, 3> forceViscous() const { return forceV; };
  void forceViscous(std::array<dfloat, 3> f) { forceV = f; }

  std::array<dfloat, 3> forcePressure() const { return forceP; };
  void forcePressure(std::array<dfloat, 3> f) { forceP = f; }

  std::array<dfloat, 3> forceEff() const 
  { 
     return {forceV[0] + forceP[0], forceV[1] + forceP[1], forceV[2] + forceP[2]}; 
  };

  void rho(const occa::memory& o_rhoIn) { o_rho = o_rhoIn; }
  occa::memory rho() const { return o_rho; }

private:
  occa::memory o_rho;

  std::array<dfloat, 3> forceV;
  std::array<dfloat, 3> forceP;
};

#endif

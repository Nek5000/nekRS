#ifndef ELLIPTIC_HPP
#define ELLIPTIC_HPP 1

#include "ellipticBcTypes.h"

struct elliptic_t;

class elliptic 
{
public:
  elliptic(const std::string& name, mesh_t *mesh, dlong fieldOffset, const std::vector<int>& EToBIn, const occa::memory& o_lambda0, const occa::memory& o_lambda1);

  ~elliptic() { };  

  void solve(const occa::memory& o_lambda0, const occa::memory& o_lambda1, const occa::memory& RHS, occa::memory x);

  std::string& name() const;
  setupAide& options();

  int Niter() const;
  void Niter(int val);

  dlong fieldOffset() const;

  bool nullSpace() const;

  dfloat initialResidual() const;
  void initialResidual(dfloat val);
  dfloat initialGuessResidual() const;
  void initialGuessResidual(dfloat val);
  dfloat finalResidual() const;
  void finalResidual(dfloat val);

  dlong Nmasked();
  occa::memory o_maskIds() const;
  occa::memory o_EToB() const;
  std::vector<int> EToB() const; 

  int Nfields() const; 

  void applyZeroNormalMask(const std::function<void(dlong Nelements, const occa::memory &o_elementList, occa::memory &o_x)>& f);
  void userPreconditioner(const std::function<void(const occa::memory &o_r, occa::memory &o_z)>& f);
  std::tuple<int, int> projectionCounters() const; 

private:
  elliptic_t *solver;
  std::vector<int> _EToB;
};

#endif

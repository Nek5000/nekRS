#include "elliptic.h"
#include "elliptic.hpp"

elliptic::elliptic(const std::string &name,
                   mesh_t *mesh,
                   dlong fieldOffset,
                   const std::vector<int> &EToBIn,
                   const occa::memory &o_lambda0,
                   const occa::memory &o_lambda1)
{
  solver = new elliptic_t();
  solver->name = name;
  solver->mesh = mesh;

  _EToB = EToBIn;
  solver->EToB = _EToB.data();

  solver->fieldOffset = fieldOffset;

  ellipticSolveSetup(solver, o_lambda0, o_lambda1);
};

void elliptic::solve(const occa::memory &o_lambda0,
                     const occa::memory &o_lambda1,
                     const occa::memory &RHS,
                     occa::memory x)
{
  ellipticSolve(solver, o_lambda0, o_lambda1, RHS, x);
};

std::string &elliptic::name() const
{
  return solver->name;
}

setupAide &elliptic::options()
{
  return solver->options;
}

int elliptic::Niter() const
{
  return solver->Niter;
};

void elliptic::Niter(int val)
{
  solver->Niter = val;
};

dlong elliptic::fieldOffset() const
{
  return solver->fieldOffset;
};

bool elliptic::nullSpace() const
{
  return solver->nullspace;
};

dfloat elliptic::initialResidual() const
{
  return solver->res00Norm;
};

void elliptic::initialResidual(dfloat val)
{
  solver->res00Norm = val;
};

dfloat elliptic::initialGuessResidual() const
{
  return solver->res0Norm;
};

void elliptic::initialGuessResidual(dfloat val)
{
  solver->res0Norm = val;
};

dfloat elliptic::finalResidual() const
{
  return solver->resNorm;
};

void elliptic::finalResidual(dfloat val)
{
  solver->resNorm = val;
};

dlong elliptic::Nmasked()
{
  return solver->Nmasked;
};

occa::memory elliptic::o_maskIds() const
{
  return solver->o_maskIds;
};

occa::memory elliptic::o_EToB() const
{
  return solver->o_EToB;
};

std::vector<int> elliptic::EToB() const
{
  return _EToB;
};

int elliptic::Nfields() const
{
  return solver->Nfields;
};

void elliptic::applyZeroNormalMask(
    const std::function<void(dlong Nelements, const occa::memory &o_elementList, occa::memory &o_x)> &f)
{
  solver->applyZeroNormalMask = f;
};

void elliptic::userPreconditioner(const std::function<void(const occa::memory &o_r, occa::memory &o_z)> &f)
{
  solver->userPreconditioner = f;
};

std::tuple<int, int> elliptic::projectionCounters() const
{
  if (solver->solutionProjection) {
    return {solver->solutionProjection->getPrevNumVecsProjection(),
            solver->solutionProjection->getMaxNumVecsProjection()};
  }
  return {0, 0};
};

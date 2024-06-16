#include <elliptic.h>
#include <ellipticApplyMask.hpp>

void ellipticApplyMask(elliptic_t *solver, occa::memory &o_x, std::string precision)
{
  auto mesh = solver->mesh;
  ellipticApplyMask(solver,
                    mesh->Nelements,
                    solver->Nmasked,
                    mesh->o_elementList,
                    solver->o_maskIds,
                    o_x,
                    precision);
}

void ellipticApplyMask(elliptic_t *solver,
                       dlong Nelements,
                       dlong Nmasked,
                       const occa::memory &o_elementList,
                       const occa::memory &o_maskIds,
                       occa::memory &o_x,
                       std::string precision)
{
  auto mesh = solver->mesh;

  if (solver->applyZeroNormalMask) {
    nekrsCheck(precision != dfloatString,
               MPI_COMM_SELF,
               EXIT_FAILURE,
               "Precision level (%s) not supported in applyZeroNormalMask\n",
               precision.c_str());
    solver->applyZeroNormalMask(Nelements, o_elementList, o_x);
  }

  if (precision != dfloatString) {
    platform->linAlg->pmask(Nmasked, o_maskIds, o_x);
  } else {
    platform->linAlg->mask(Nmasked, o_maskIds, o_x);
  }
}

#include <elliptic.h>
#include <ellipticApplyMask.hpp>
void ellipticApplyMask(elliptic_t *solver, occa::memory &o_x, std::string precision)
{
  mesh_t *mesh = solver->mesh;
  ellipticApplyMask(solver, mesh->Nelements, solver->Nmasked, mesh->o_elementList, solver->o_maskIds, o_x, precision);
}
void ellipticApplyMask(elliptic_t *solver,
                       dlong Nelements,
                       dlong Nmasked,
                       occa::memory &o_elementList,
                       occa::memory &o_maskIds,
                       occa::memory &o_x,
                       std::string precision)
{
  mesh_t *mesh = solver->mesh;
  occa::kernel &maskKernel = (precision != dfloatString) ? mesh->maskPfloatKernel : mesh->maskKernel;

  if (solver->applyZeroNormalMask) {
    nrsCheck(precision != dfloatString, MPI_COMM_SELF, EXIT_FAILURE,
             "Precision level (%s) not supported in applyZeroNormalMask\n", precision.c_str());
    solver->applyZeroNormalMask(Nelements, o_elementList, o_x);
  }
  if (Nmasked) {
    maskKernel(Nmasked, o_maskIds, o_x);
  }
}

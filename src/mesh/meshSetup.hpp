#if !defined(nekrs_meshsetup_hpp_)
#define nekrs_meshsetup_hpp_

#include "nrs.hpp"
mesh_t *createMesh(
                MPI_Comm comm,
                int N,
                int cubN,
                bool cht,
                occa::properties &kernelInfo);
#endif

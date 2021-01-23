#if !defined(nekrs_meshsetup_hpp_)
#define nekrs_meshsetup_hpp_

#include "nrs.hpp"
mesh_t* createMeshDummy(MPI_Comm comm,
                        int N,
                        int cubN,
                        occa::device device,
                        occa::properties &kernelInfo);

mesh_t* createMesh(MPI_Comm comm,
                   int N,
                   int cubN,
                   int isMeshT,
                   occa::device device,
                   occa::properties &kernelInfo);

mesh_t* createMeshV(MPI_Comm comm,
                    int N,
                    int cubN,
                    mesh_t* meshT,
                    occa::properties &kernelInfo);

#endif

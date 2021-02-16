#if !defined(nekrs_meshsetup_hpp_)
#define nekrs_meshsetup_hpp_

#include "nrs.hpp"
mesh_t* createMeshDummy(MPI_Comm comm,
                        int N,
                        int cubN,
                        setupAide &options,
                        occa::properties &kernelInfo);

mesh_t* createMesh(MPI_Comm comm,
                   int N,
                   int cubN,
                   int isMeshT,
                   setupAide &options,
                   occa::properties &kernelInfo);

mesh_t* createMeshV(MPI_Comm comm,
                    int N,
                    int cubN,
                    mesh_t* meshT,
                    setupAide &options,
                    occa::properties &kernelInfo);

#endif

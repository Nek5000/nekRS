#if !defined(nekrs_meshsetup_hpp_)
#define nekrs_meshsetup_hpp_

#include "nekrs.hpp"
mesh_t* createMeshDummy(MPI_Comm comm,
                        int N,
                        int cubN,
                        setupAide &options,
                        occa::device device,
                        occa::properties &kernelInfo);
mesh_t* createMeshT(MPI_Comm comm,
                    int N,
                    int cubN,
                    int isMeshT,
                    setupAide &options,
                    occa::device device,
                    occa::properties &kernelInfo);
mesh_t* createMeshV(MPI_Comm comm,
                    int N,
                    int cubN,
                    mesh_t* meshT,
                    setupAide &options,
                    occa::properties &kernelInfo);

#endif

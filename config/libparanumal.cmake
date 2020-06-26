# TODO: Need to get value from OCCA?
include(FetchContent)
set(FETCHCONTENT_QUIET OFF)

set(USE_OCCA_MEM_BYTE_ALIGN 64)

add_definitions( -DUSE_OCCA_MEM_BYTE_ALIGN=${USE_OCCA_MEM_BYTE_ALIGN} )

# ---------------------------------------------------------
# Download dependencies
# ---------------------------------------------------------

FetchContent_Declare(
    paranumal_content
    GIT_REPOSITORY https://gitlab.com/nekrs/libparanumal.git
    GIT_TAG next)
FetchContent_GetProperties(paranumal_content)
if(NOT paranumal_content_POPULATED)
  FetchContent_Populate(paranumal_content)
endif()

# FetchContent_MakeAvailable uses the CMakeLists in SOURCE_SUBDIR
FetchContent_Declare(
  hypre_content
  URL https://github.com/hypre-space/hypre/archive/v2.18.2.tar.gz )
#FetchContent_MakeAvailable(hypre_content)
FetchContent_GetProperties(hypre_content)
if (NOT hypre_content_POPULATED)
  FetchContent_Populate(hypre_content)
  add_subdirectory(${hypre_content_SOURCE_DIR}/src ${hypre_content_BINARY_DIR})
endif()

set(PARANUMAL_DIR ${paranumal_content_SOURCE_DIR})
set(OGS_DIR ${PARANUMAL_DIR}/libs/gatherScatter)
set(PARALMOND_DIR ${PARANUMAL_DIR}/libs/parAlmond)
set(ELLIPTIC_DIR ${PARANUMAL_DIR}/solvers/elliptic)
set(HYPRE_DIR ${hypre_content_SOURCE_DIR})

# ---------------------------------------------------------
# libogs
# ---------------------------------------------------------

set(OGS_SOURCES
        ${OGS_DIR}/src/ogsGather.cpp
        ${OGS_DIR}/src/ogsGatherMany.cpp
        ${OGS_DIR}/src/ogsGatherScatter.cpp
        ${OGS_DIR}/src/ogsGatherScatterMany.cpp
        ${OGS_DIR}/src/ogsGatherScatterVec.cpp
        ${OGS_DIR}/src/ogsGatherVec.cpp
        ${OGS_DIR}/src/ogsHostGather.c
        ${OGS_DIR}/src/ogsHostGatherMany.c
        ${OGS_DIR}/src/ogsHostGatherScatter.c
        ${OGS_DIR}/src/ogsHostGatherScatterMany.c
        ${OGS_DIR}/src/ogsHostGatherScatterVec.c
        ${OGS_DIR}/src/ogsHostGatherVec.c
        ${OGS_DIR}/src/ogsHostScatter.c
        ${OGS_DIR}/src/ogsHostScatterMany.c
        ${OGS_DIR}/src/ogsHostScatterVec.c
        ${OGS_DIR}/src/ogsHostSetup.c
        ${OGS_DIR}/src/ogsKernels.cpp
        ${OGS_DIR}/src/ogsMappedAlloc.cpp
        ${OGS_DIR}/src/ogsScatter.cpp
        ${OGS_DIR}/src/ogsScatterMany.cpp
        ${OGS_DIR}/src/ogsScatterVec.cpp
        ${OGS_DIR}/src/ogsSetup.cpp)

add_library(libogs ${OGS_SOURCES})
set_target_properties(libogs PROPERTIES OUTPUT_NAME ogs)
target_compile_definitions(libogs PUBLIC -DDOGS="${OGS_DIR}")
target_include_directories(libogs PUBLIC
        ${OGS_DIR}/include
        ${OGS_DIR}
        ${PARANUMAL_DIR}/include)
target_link_libraries(libogs PUBLIC libocca gs)

# ---------------------------------------------------------
# libparanumal
# ---------------------------------------------------------

set(PARANUMAL_SOURCES
        src/core/occaDeviceConfig.cpp
        ${PARANUMAL_DIR}/src/hash.c
        ${PARANUMAL_DIR}/src/matrixConditionNumber.c
        ${PARANUMAL_DIR}/src/matrixInverse.c
        ${PARANUMAL_DIR}/src/meshApplyElementMatrix.c
        ${PARANUMAL_DIR}/src/meshConnect.c
        ${PARANUMAL_DIR}/src/meshConnectBoundary.c
        ${PARANUMAL_DIR}/src/meshConnectFaceNodes2D.c
        ${PARANUMAL_DIR}/src/meshConnectFaceNodes3D.c
        ${PARANUMAL_DIR}/src/meshConnectPeriodicFaceNodes2D.c
        ${PARANUMAL_DIR}/src/meshConnectPeriodicFaceNodes3D.c
        ${PARANUMAL_DIR}/src/meshGeometricFactorsHex3D.c
        ${PARANUMAL_DIR}/src/meshGeometricFactorsQuad2D.c
        ${PARANUMAL_DIR}/src/meshGeometricFactorsQuad3D.c
        ${PARANUMAL_DIR}/src/meshGeometricFactorsTet3D.c
        ${PARANUMAL_DIR}/src/meshGeometricFactorsTri2D.c
        ${PARANUMAL_DIR}/src/meshGeometricFactorsTri3D.c
        ${PARANUMAL_DIR}/src/meshGeometricPartition2D.c
        ${PARANUMAL_DIR}/src/meshGeometricPartition3D.c
        ${PARANUMAL_DIR}/src/meshHaloExchange.c
        ${PARANUMAL_DIR}/src/meshHaloExtract.c
        ${PARANUMAL_DIR}/src/meshHaloSetup.c
        ${PARANUMAL_DIR}/src/meshLoadReferenceNodesHex3D.c
        ${PARANUMAL_DIR}/src/meshLoadReferenceNodesQuad2D.c
        ${PARANUMAL_DIR}/src/meshLoadReferenceNodesTet3D.c
        ${PARANUMAL_DIR}/src/meshLoadReferenceNodesTri2D.c
        ${PARANUMAL_DIR}/src/meshOccaSetup2D.c
        ${PARANUMAL_DIR}/src/meshOccaSetup3D.c
        ${PARANUMAL_DIR}/src/meshOccaSetupQuad3D.c
        ${PARANUMAL_DIR}/src/meshParallelConnectNodes.c
        ${PARANUMAL_DIR}/src/meshParallelConnectOpt.c
        ${PARANUMAL_DIR}/src/meshParallelConsecutiveGlobalNumbering.c
        ${PARANUMAL_DIR}/src/meshParallelGatherScatterSetup.c
        ${PARANUMAL_DIR}/src/meshParallelReaderHex3D.c
        ${PARANUMAL_DIR}/src/meshParallelReaderQuad2D.c
        ${PARANUMAL_DIR}/src/meshParallelReaderQuad3D.c
        ${PARANUMAL_DIR}/src/meshParallelReaderTet3D.c
        ${PARANUMAL_DIR}/src/meshParallelReaderTri2D.c
        ${PARANUMAL_DIR}/src/meshParallelReaderTri3D.c
        ${PARANUMAL_DIR}/src/meshPartitionStatistics.c
        ${PARANUMAL_DIR}/src/meshPhysicalNodesQuad2D.c
        ${PARANUMAL_DIR}/src/meshPhysicalNodesQuad3D.c
        ${PARANUMAL_DIR}/src/meshPhysicalNodesTet3D.c
        ${PARANUMAL_DIR}/src/meshPhysicalNodesTri2D.c
        ${PARANUMAL_DIR}/src/meshPhysicalNodesTri3D.c
        ${PARANUMAL_DIR}/src/meshPlotVTU2D.c
        ${PARANUMAL_DIR}/src/meshPlotVTU3D.c
        ${PARANUMAL_DIR}/src/meshPrint2D.c
        ${PARANUMAL_DIR}/src/meshPrint3D.c
        ${PARANUMAL_DIR}/src/meshSetup.c
        ${PARANUMAL_DIR}/src/meshSetupBoxHex3D.c
        ${PARANUMAL_DIR}/src/meshSetupBoxQuad2D.c
        ${PARANUMAL_DIR}/src/meshSetupHex3D.c
        ${PARANUMAL_DIR}/src/meshSetupQuad2D.c
        ${PARANUMAL_DIR}/src/meshSetupQuad3D.c
        ${PARANUMAL_DIR}/src/meshSetupTet3D.c
        ${PARANUMAL_DIR}/src/meshSetupTri2D.c
        ${PARANUMAL_DIR}/src/meshSetupTri3D.c
        ${PARANUMAL_DIR}/src/meshSurfaceGeometricFactorsHex3D.c
        ${PARANUMAL_DIR}/src/meshSurfaceGeometricFactorsQuad2D.c
        ${PARANUMAL_DIR}/src/meshSurfaceGeometricFactorsQuad3D.c
        ${PARANUMAL_DIR}/src/meshSurfaceGeometricFactorsTet3D.c
        ${PARANUMAL_DIR}/src/meshSurfaceGeometricFactorsTri2D.c
        ${PARANUMAL_DIR}/src/meshSurfaceGeometricFactorsTri3D.c
        ${PARANUMAL_DIR}/src/meshVTU2D.c
        ${PARANUMAL_DIR}/src/meshVTU3D.c
        ${PARANUMAL_DIR}/src/mysort.c
        ${PARANUMAL_DIR}/src/occaHostMallocPinned.c
        ${PARANUMAL_DIR}/src/parallelSort.c
        ${PARANUMAL_DIR}/src/readArray.c
        ${PARANUMAL_DIR}/src/setupAide.c
        ${PARANUMAL_DIR}/src/timer.c)

set_source_files_properties(${PARANUMAL_SOURCES} PROPERTIES LANGUAGE CXX)
add_library(libparanumal ${PARANUMAL_SOURCES})
set_target_properties(libparanumal PROPERTIES OUTPUT_NAME paranumal)
target_compile_definitions(libparanumal PUBLIC -DDHOLMES="${PARANUMAL_DIR}")
target_include_directories(libparanumal PUBLIC ${PARANUMAL_DIR}/include src/core/)
target_link_libraries(libparanumal PUBLIC libogs libocca blasLapack)

# ---------------------------------------------------------
# libparAlmond
# ---------------------------------------------------------

set(PARALMOND_SOURCES
        ${PARALMOND_DIR}/hypre/hypre.c
        ${PARALMOND_DIR}/src/SpMV.cpp
        ${PARALMOND_DIR}/src/agmgLevel.cpp
        ${PARALMOND_DIR}/src/agmgSetup/agmgSetup.cpp
        ${PARALMOND_DIR}/src/agmgSetup/constructProlongation.cpp
        ${PARALMOND_DIR}/src/agmgSetup/formAggregates.cpp
        ${PARALMOND_DIR}/src/agmgSetup/galerkinProd.cpp
        ${PARALMOND_DIR}/src/agmgSetup/strongGraph.cpp
        ${PARALMOND_DIR}/src/agmgSetup/transpose.cpp
        ${PARALMOND_DIR}/src/agmgSmoother.cpp
        ${PARALMOND_DIR}/src/coarseSolver.cpp
        ${PARALMOND_DIR}/src/kernels.cpp
        ${PARALMOND_DIR}/src/level.cpp
        ${PARALMOND_DIR}/src/matrix.cpp
        ${PARALMOND_DIR}/src/multigrid.cpp
        ${PARALMOND_DIR}/src/parAlmond.cpp
        ${PARALMOND_DIR}/src/pcg.cpp
        ${PARALMOND_DIR}/src/pgmres.cpp
        ${PARALMOND_DIR}/src/solver.cpp
        ${PARALMOND_DIR}/src/utils.cpp
        ${PARALMOND_DIR}/src/vector.cpp)

add_library(libparAlmond ${PARALMOND_SOURCES})
set_target_properties(libparAlmond PROPERTIES OUTPUT_NAME parAlmond)
target_compile_definitions(libparAlmond PUBLIC -DDPARALMOND="${PARALMOND_DIR}" PRIVATE -DHYPRE)
target_include_directories(libparAlmond 
        PUBLIC
        ${PARALMOND_DIR}/include
        ${PARALMOND_DIR}
        ${PARALMOND_DIR}/hypre
        ${PARANUMAL_DIR}/include
        PRIVATE
        ${HYPRE_DIR}/src
        ${HYPRE_DIR}/src/utilities
        ${HYPRE_DIR}/src/seq_mv
        ${HYPRE_DIR}/src/parcsr_mv
        ${HYPRE_DIR}/src/parcsr_ls
        ${HYPRE_DIR}/src/IJ_mv
        ${HYPRE_DIR}/src/multivector
        ${HYPRE_DIR}/src/krylov
        ${hypre_content_BINARY_DIR})
target_link_libraries(libparAlmond PUBLIC libogs libocca PRIVATE HYPRE)
# This conflicts with the stdlib "version" header...
file(REMOVE ${HYPRE_DIR}/src/utilities/version)

# ---------------------------------------------------------
# libelliptic
# ---------------------------------------------------------

set(ELLIPTIC_SOURCES
        ${ELLIPTIC_DIR}/src/NBFPCG.c
        ${ELLIPTIC_DIR}/src/NBPCG.c
        ${ELLIPTIC_DIR}/src/PCG.c
        ${ELLIPTIC_DIR}/src/ellipticBuildContinuous.c
        ${ELLIPTIC_DIR}/src/ellipticBuildContinuousGalerkin.c
        ${ELLIPTIC_DIR}/src/ellipticBuildIpdg.c
        ${ELLIPTIC_DIR}/src/ellipticBuildJacobi.c
        ${ELLIPTIC_DIR}/src/ellipticBuildLocalPatches.c
        ${ELLIPTIC_DIR}/src/ellipticBuildMultigridLevel.c
        ${ELLIPTIC_DIR}/src/ellipticHaloExchange.c
        ${ELLIPTIC_DIR}/src/ellipticKernelInfo.c
        ${ELLIPTIC_DIR}/src/ellipticMixedCopy.c
        ${ELLIPTIC_DIR}/src/ellipticMultiGridLevel.c
        ${ELLIPTIC_DIR}/src/ellipticMultiGridLevelSetup.c
        ${ELLIPTIC_DIR}/src/ellipticMultiGridSetup.c
        ${ELLIPTIC_DIR}/src/ellipticOperator.c
        ${ELLIPTIC_DIR}/src/ellipticPlotVTUHex3D.c
        ${ELLIPTIC_DIR}/src/ellipticPreconditioner.c
        ${ELLIPTIC_DIR}/src/ellipticPreconditionerSetup.c
        ${ELLIPTIC_DIR}/src/ellipticSEMFEMSetup.c
        ${ELLIPTIC_DIR}/src/ellipticScaledAdd.c
        ${ELLIPTIC_DIR}/src/ellipticSetScalar.c
        ${ELLIPTIC_DIR}/src/ellipticSolve.c
        ${ELLIPTIC_DIR}/src/ellipticSolveSetup.c
        ${ELLIPTIC_DIR}/src/ellipticThinOas.c
        ${ELLIPTIC_DIR}/src/ellipticThinOasSetup.c
        ${ELLIPTIC_DIR}/src/ellipticUpdateNBFPCG.c
        ${ELLIPTIC_DIR}/src/ellipticUpdateNBPCG.c
        ${ELLIPTIC_DIR}/src/ellipticUpdatePCG.c
        ${ELLIPTIC_DIR}/src/ellipticVectors.c
        ${ELLIPTIC_DIR}/src/ellipticWeightedInnerProduct.c
        ${ELLIPTIC_DIR}/src/ellipticWeightedNorm2.c
        ${ELLIPTIC_DIR}/src/ellipticZeroMean.c)

set_source_files_properties(${ELLIPTIC_SOURCES} PROPERTIES LANGUAGE CXX)
add_library(libelliptic ${ELLIPTIC_SOURCES})
set_target_properties(libelliptic PROPERTIES OUTPUT_NAME elliptic)
target_compile_definitions(libelliptic PUBLIC -DDELLIPTIC="${ELLIPTIC_DIR}")
target_include_directories(libelliptic PUBLIC ${ELLIPTIC_DIR})
target_link_libraries(libelliptic PUBLIC libparanumal libparAlmond libogs libocca blasLapack)


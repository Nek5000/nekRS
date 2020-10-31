set(LIBP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/libP)
set(PARALMOND_SOURCE_DIR ${LIBP_SOURCE_DIR}/parAlmond)
set(ELLIPTIC_SOURCE_DIR ${LIBP_SOURCE_DIR}/solvers/elliptic)

set(LIBP_SOURCES
        ${LIBP_SOURCE_DIR}/src/hash.c
        ${LIBP_SOURCE_DIR}/src/matrixConditionNumber.c
        ${LIBP_SOURCE_DIR}/src/matrixInverse.c
        ${LIBP_SOURCE_DIR}/src/matrixEig.cpp
        ${LIBP_SOURCE_DIR}/src/matrixTranspose.cpp
        ${LIBP_SOURCE_DIR}/src/matrixRightSolve.cpp
        ${LIBP_SOURCE_DIR}/src/meshBasis1D.cpp
        ${LIBP_SOURCE_DIR}/src/meshBasisHex3D.cpp
        ${LIBP_SOURCE_DIR}/src/meshApplyElementMatrix.c
        ${LIBP_SOURCE_DIR}/src/meshConnect.c
        ${LIBP_SOURCE_DIR}/src/meshConnectBoundary.c
        ${LIBP_SOURCE_DIR}/src/meshConnectFaceNodes2D.c
        ${LIBP_SOURCE_DIR}/src/meshConnectFaceNodes3D.c
        ${LIBP_SOURCE_DIR}/src/meshConnectPeriodicFaceNodes2D.c
        ${LIBP_SOURCE_DIR}/src/meshConnectPeriodicFaceNodes3D.c
        ${LIBP_SOURCE_DIR}/src/meshFree.c
        ${LIBP_SOURCE_DIR}/src/meshGeometricFactorsHex3D.c
        ${LIBP_SOURCE_DIR}/src/meshGeometricFactorsQuad2D.c
        ${LIBP_SOURCE_DIR}/src/meshGeometricFactorsQuad3D.c
        ${LIBP_SOURCE_DIR}/src/meshGeometricFactorsTet3D.c
        ${LIBP_SOURCE_DIR}/src/meshGeometricFactorsTri2D.c
        ${LIBP_SOURCE_DIR}/src/meshGeometricFactorsTri3D.c
        ${LIBP_SOURCE_DIR}/src/meshGeometricPartition2D.c
        ${LIBP_SOURCE_DIR}/src/meshGeometricPartition3D.c
        ${LIBP_SOURCE_DIR}/src/meshHaloExchange.c
        ${LIBP_SOURCE_DIR}/src/meshHaloExtract.c
        ${LIBP_SOURCE_DIR}/src/meshHaloSetup.c
        ${LIBP_SOURCE_DIR}/src/meshLoadReferenceNodesHex3D.c
        ${LIBP_SOURCE_DIR}/src/meshLoadReferenceNodesQuad2D.c
        ${LIBP_SOURCE_DIR}/src/meshLoadReferenceNodesTet3D.c
        ${LIBP_SOURCE_DIR}/src/meshLoadReferenceNodesTri2D.c
        ${LIBP_SOURCE_DIR}/src/meshOccaSetup2D.c
        ${LIBP_SOURCE_DIR}/src/meshOccaSetup3D.c
        ${LIBP_SOURCE_DIR}/src/meshOccaSetupQuad3D.c
        ${LIBP_SOURCE_DIR}/src/meshParallelConnectOpt.c
        ${LIBP_SOURCE_DIR}/src/meshParallelConsecutiveGlobalNumbering.c
        ${LIBP_SOURCE_DIR}/src/meshParallelGatherScatterSetup.c
        ${LIBP_SOURCE_DIR}/src/meshParallelReaderHex3D.c
        ${LIBP_SOURCE_DIR}/src/meshParallelReaderQuad2D.c
        ${LIBP_SOURCE_DIR}/src/meshParallelReaderQuad3D.c
        ${LIBP_SOURCE_DIR}/src/meshParallelReaderTet3D.c
        ${LIBP_SOURCE_DIR}/src/meshParallelReaderTri2D.c
        ${LIBP_SOURCE_DIR}/src/meshParallelReaderTri3D.c
        ${LIBP_SOURCE_DIR}/src/meshPartitionStatistics.c
        ${LIBP_SOURCE_DIR}/src/meshPhysicalNodesQuad2D.c
        ${LIBP_SOURCE_DIR}/src/meshPhysicalNodesQuad3D.c
        ${LIBP_SOURCE_DIR}/src/meshPhysicalNodesTet3D.c
        ${LIBP_SOURCE_DIR}/src/meshPhysicalNodesTri2D.c
        ${LIBP_SOURCE_DIR}/src/meshPhysicalNodesTri3D.c
        ${LIBP_SOURCE_DIR}/src/meshSurfaceGeometricFactorsHex3D.c
        ${LIBP_SOURCE_DIR}/src/meshSurfaceGeometricFactorsQuad2D.c
        ${LIBP_SOURCE_DIR}/src/meshSurfaceGeometricFactorsQuad3D.c
        ${LIBP_SOURCE_DIR}/src/meshSurfaceGeometricFactorsTet3D.c
        ${LIBP_SOURCE_DIR}/src/meshSurfaceGeometricFactorsTri2D.c
        ${LIBP_SOURCE_DIR}/src/meshSurfaceGeometricFactorsTri3D.c
        ${LIBP_SOURCE_DIR}/src/mysort.c
        ${LIBP_SOURCE_DIR}/src/occaHostMallocPinned.c
        ${LIBP_SOURCE_DIR}/src/parallelSort.c
        ${LIBP_SOURCE_DIR}/src/readArray.c
        ${LIBP_SOURCE_DIR}/src/setupAide.c
        ${LIBP_SOURCE_DIR}/src/timer.c)

set_source_files_properties(${LIBP_SOURCES} PROPERTIES LANGUAGE CXX)

set(PARALMOND_SOURCES
        ${PARALMOND_SOURCE_DIR}/hypre/hypre.c
        ${PARALMOND_SOURCE_DIR}/src/SpMV.cpp
        ${PARALMOND_SOURCE_DIR}/src/agmgLevel.cpp
        ${PARALMOND_SOURCE_DIR}/src/agmgSetup/agmgSetup.cpp
        ${PARALMOND_SOURCE_DIR}/src/agmgSetup/constructProlongation.cpp
        ${PARALMOND_SOURCE_DIR}/src/agmgSetup/formAggregates.cpp
        ${PARALMOND_SOURCE_DIR}/src/agmgSetup/galerkinProd.cpp
        ${PARALMOND_SOURCE_DIR}/src/agmgSetup/strongGraph.cpp
        ${PARALMOND_SOURCE_DIR}/src/agmgSetup/transpose.cpp
        ${PARALMOND_SOURCE_DIR}/src/agmgSmoother.cpp
        ${PARALMOND_SOURCE_DIR}/src/coarseSolver.cpp
        ${PARALMOND_SOURCE_DIR}/src/kernels.cpp
        ${PARALMOND_SOURCE_DIR}/src/level.cpp
        ${PARALMOND_SOURCE_DIR}/src/matrix.cpp
        ${PARALMOND_SOURCE_DIR}/src/multigrid.cpp
        ${PARALMOND_SOURCE_DIR}/src/parAlmond.cpp
        ${PARALMOND_SOURCE_DIR}/src/pcg.cpp
        ${PARALMOND_SOURCE_DIR}/src/pgmres.cpp
        ${PARALMOND_SOURCE_DIR}/src/solver.cpp
        ${PARALMOND_SOURCE_DIR}/src/timer.cpp
        ${PARALMOND_SOURCE_DIR}/src/utils.cpp
        ${PARALMOND_SOURCE_DIR}/src/vector.cpp)

# ---------------------------------------------------------
# libelliptic
# ---------------------------------------------------------

set(ELLIPTIC_SOURCES
        ${ELLIPTIC_SOURCE_DIR}/src/NBFPCG.c
        ${ELLIPTIC_SOURCE_DIR}/src/NBPCG.c
        ${ELLIPTIC_SOURCE_DIR}/src/PCG.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticBuildContinuous.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticBuildContinuousGalerkin.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticBuildIpdg.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticBuildJacobi.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticHaloExchange.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticKernelInfo.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticBuildMultigridLevelFine.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticBuildMultigridLevel.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticMultiGridLevel.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticMultiGridLevelSetup.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticMultiGridSchwarz.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticMultiGridSetup.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticOperator.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticPreconditioner.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticPreconditionerSetup.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticResidualProjection.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticScaledAdd.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticSolve.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticSolveSetup.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticUpdateNBFPCG.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticUpdateNBPCG.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticUpdatePCG.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticVectors.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticWeightedInnerProduct.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticWeightedNorm2.c
        ${ELLIPTIC_SOURCE_DIR}/src/ellipticZeroMean.c)

set_source_files_properties(${ELLIPTIC_SOURCES} PROPERTIES LANGUAGE CXX)

# ---------------------------------------------------------
# install
# ---------------------------------------------------------

set(file_pattern "\.cu$|\.hip$|\.okl$|\.c$|\.hpp$|\.tpp$|\.h$|hex.*\.dat$")

install(DIRECTORY 
  ${LIBP_SOURCE_DIR}/include 
  ${LIBP_SOURCE_DIR}/nodes
  ${LIBP_SOURCE_DIR}/okl 
  DESTINATION libparanumal
  FILES_MATCHING REGEX ${file_pattern})

install(DIRECTORY
  ${PARALMOND_SOURCE_DIR}/include
  ${PARALMOND_SOURCE_DIR}/okl
  DESTINATION parAlmond
  FILES_MATCHING REGEX ${file_pattern})
install(FILES ${PARALMOND_SOURCE_DIR}/parAlmond.hpp DESTINATION gatherScatter)

install(DIRECTORY
  ${ELLIPTIC_SOURCE_DIR}/data
  ${ELLIPTIC_SOURCE_DIR}/okl
  DESTINATION elliptic
  FILES_MATCHING REGEX ${file_pattern})
install(FILES 
  ${ELLIPTIC_SOURCE_DIR}/elliptic.h
  ${ELLIPTIC_SOURCE_DIR}/ellipticMultiGrid.h
  ${ELLIPTIC_SOURCE_DIR}/ellipticResidualProjection.h
  ${ELLIPTIC_SOURCE_DIR}/ellipticPrecon.h
  DESTINATION gatherScatter)

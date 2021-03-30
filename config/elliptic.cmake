set(ELLIPTIC_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/elliptic)

set(ELLIPTIC_SOURCES
        ${ELLIPTIC_SOURCE_DIR}/linearSolver/PCG.cpp
        ${ELLIPTIC_SOURCE_DIR}/linearSolver/PGMRES.cpp
	      ${ELLIPTIC_SOURCE_DIR}/ellipticBuildContinuous.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticBuildContinuousGalerkin.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticJacobi.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticKernelInfo.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticBuildMultigridLevelFine.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticBuildMultigridLevel.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticMultiGridLevel.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticMultiGridLevelSetup.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticMultiGridSchwarz.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticMultiGridSetup.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticOperator.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticPreconditioner.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticPreconditionerSetup.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticResidualProjection.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticSolve.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticSolveSetup.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticUpdatePCG.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticZeroMean.cpp)

set(PARALMOND_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/elliptic/parAlmond)

set(PARALMOND_SOURCES
        ${PARALMOND_SOURCE_DIR}/crs_hypre.cpp
        ${PARALMOND_SOURCE_DIR}/SpMV.cpp
        ${PARALMOND_SOURCE_DIR}/agmgLevel.cpp
        ${PARALMOND_SOURCE_DIR}/agmgSetup/agmgSetup.cpp
        ${PARALMOND_SOURCE_DIR}/agmgSetup/constructProlongation.cpp
        ${PARALMOND_SOURCE_DIR}/agmgSetup/formAggregates.cpp
        ${PARALMOND_SOURCE_DIR}/agmgSetup/galerkinProd.cpp
        ${PARALMOND_SOURCE_DIR}/agmgSetup/strongGraph.cpp
        ${PARALMOND_SOURCE_DIR}/agmgSetup/transpose.cpp
        ${PARALMOND_SOURCE_DIR}/agmgSmoother.cpp
        ${PARALMOND_SOURCE_DIR}/coarseSolver.cpp
        ${PARALMOND_SOURCE_DIR}/kernels.cpp
        ${PARALMOND_SOURCE_DIR}/level.cpp
        ${PARALMOND_SOURCE_DIR}/matrix.cpp
        ${PARALMOND_SOURCE_DIR}/multigrid.cpp
        ${PARALMOND_SOURCE_DIR}/parAlmond.cpp
        ${PARALMOND_SOURCE_DIR}/pcg.cpp
        ${PARALMOND_SOURCE_DIR}/pgmres.cpp
        ${PARALMOND_SOURCE_DIR}/solver.cpp
        ${PARALMOND_SOURCE_DIR}/timer.cpp
        ${PARALMOND_SOURCE_DIR}/utils.cpp
        ${PARALMOND_SOURCE_DIR}/vector.cpp)

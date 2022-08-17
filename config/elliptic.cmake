set(ELLIPTIC_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/elliptic)

set(ELLIPTIC_SOURCES
        ${ELLIPTIC_SOURCE_DIR}/linearSolver/PCG.cpp
        ${ELLIPTIC_SOURCE_DIR}/linearSolver/PGMRES.cpp
        ${ELLIPTIC_SOURCE_DIR}/amgSolver/amgx/amgx.c
        ${ELLIPTIC_SOURCE_DIR}/ellipticApplyMask.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticBuildSEMFEM.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticBuildContinuous.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticBuildContinuousGalerkin.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticUpdateJacobi.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticSEMFEM.cpp
        ${ELLIPTIC_SOURCE_DIR}/registerEllipticKernels.cpp
        ${ELLIPTIC_SOURCE_DIR}/registerEllipticPreconditionerKernels.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticBuildPreconditionerKernels.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticBuildMultigridLevelFine.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticBuildMultigridLevel.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticMultiGridUpdateLambda.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticMultiGridLevel.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticMultiGridLevelSetup.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticMultiGridSchwarz.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticMultiGridSetup.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticOperator.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticPreconditioner.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticPreconditionerSetup.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticSolutionProjection.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticSolve.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticOgs.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticSetup.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticUpdatePCG.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticZeroMean.cpp)

set(PARALMOND_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/elliptic/amgSolver/parAlmond)

set(PARALMOND_SOURCES
        ${PARALMOND_SOURCE_DIR}/coarseSolver.cpp
#        ${PARALMOND_SOURCE_DIR}/kernels.cpp
        ${PARALMOND_SOURCE_DIR}/level.cpp
        ${PARALMOND_SOURCE_DIR}/multigrid.cpp
        ${PARALMOND_SOURCE_DIR}/parAlmond.cpp
        ${PARALMOND_SOURCE_DIR}/solver.cpp
)

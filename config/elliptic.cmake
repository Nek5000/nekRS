set(ELLIPTIC_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/solvers/elliptic)

set(ELLIPTIC_SOURCES
        ${ELLIPTIC_SOURCE_DIR}/linearSolver/PCG.cpp
        ${ELLIPTIC_SOURCE_DIR}/linearSolver/PGMRES.cpp
        ${ELLIPTIC_SOURCE_DIR}/amgSolver/amgx/AMGX.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticApplyMask.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticUpdateJacobi.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticBuildPreconditionerKernels.cpp
        ${ELLIPTIC_SOURCE_DIR}/MG/ellipticBuildMultigridLevelFine.cpp
        ${ELLIPTIC_SOURCE_DIR}/MG/ellipticBuildMultigridLevel.cpp
        ${ELLIPTIC_SOURCE_DIR}/MG/ellipticMultiGridUpdateLambda.cpp
        ${ELLIPTIC_SOURCE_DIR}/MG/optimalCoeffs.cpp
        ${ELLIPTIC_SOURCE_DIR}/MG/determineMGLevels.cpp
        ${ELLIPTIC_SOURCE_DIR}/MG/parseMultigridSchedule.cpp
        ${ELLIPTIC_SOURCE_DIR}/MG/ellipticMultiGridLevel.cpp
        ${ELLIPTIC_SOURCE_DIR}/MG/ellipticMultiGridLevelSetup.cpp
        ${ELLIPTIC_SOURCE_DIR}/MG/ellipticMultiGridSchwarz.cpp
        ${ELLIPTIC_SOURCE_DIR}/MG/ellipticMultiGridSetup.cpp
        ${ELLIPTIC_SOURCE_DIR}/MG/ellipticBuildFEM.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticOperator.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticPreconditioner.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticPreconditionerSetup.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticSolutionProjection.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticSolve.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticOgs.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticSetup.cpp
        ${ELLIPTIC_SOURCE_DIR}/SEMFEMSolver.cpp
        ${ELLIPTIC_SOURCE_DIR}/SEMFEMSolverBuild.cpp
        ${ELLIPTIC_SOURCE_DIR}/MG/coarseLevel.cpp
        ${ELLIPTIC_SOURCE_DIR}/MG/level.cpp
        ${ELLIPTIC_SOURCE_DIR}/MG/MGSolver.cpp
        ${ELLIPTIC_SOURCE_DIR}/ellipticZeroMean.cpp)

set(PARALMOND_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/parAlmond)

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

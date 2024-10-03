set(MESH_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/mesh)

set(MESH_SOURCES
    ${MESH_SOURCE_DIR}/planarAvg.cpp
    ${MESH_SOURCE_DIR}/meshSetup.cpp
    ${MESH_SOURCE_DIR}/meshIntp.cpp
    ${MESH_SOURCE_DIR}/meshSurface.cpp
    ${MESH_SOURCE_DIR}/meshDistance.cpp
    ${MESH_SOURCE_DIR}/meshNekReader.cpp
    ${MESH_SOURCE_DIR}/meshPhysicalNodesHex3D.cpp
    ${MESH_SOURCE_DIR}/meshGlobalIds.cpp
    ${MESH_SOURCE_DIR}/meshBasis1D.cpp
    ${MESH_SOURCE_DIR}/meshBasisHex3D.cpp
    ${MESH_SOURCE_DIR}/meshApplyElementMatrix.cpp
    ${MESH_SOURCE_DIR}/meshMetrics.cpp
    ${MESH_SOURCE_DIR}/meshConnect.cpp
    ${MESH_SOURCE_DIR}/meshConnectFaceNodes3D.cpp
    ${MESH_SOURCE_DIR}/meshFree.cpp
    ${MESH_SOURCE_DIR}/meshMove.cpp
    ${MESH_SOURCE_DIR}/meshGeometricFactorsHex3D.cpp
    ${MESH_SOURCE_DIR}/meshLoadReferenceNodesHex3D.cpp
    ${MESH_SOURCE_DIR}/meshOccaSetup3D.cpp
    ${MESH_SOURCE_DIR}/meshParallelConsecutiveGlobalNumbering.cpp
    ${MESH_SOURCE_DIR}/meshParallelGatherScatterSetup.cpp
    ${MESH_SOURCE_DIR}/meshSurfaceGeometricFactorsHex3D.cpp
    ${MESH_SOURCE_DIR}/meshComputeInvLMM.cpp
    ${MESH_SOURCE_DIR}/meshParallelConnectOpt.cpp)

include(cmake/platform.cmake)
include(cmake/core.cmake)
include(cmake/solver.cmake)

set(APP_SOURCES_DIR "src/app")
set(APP_SOURCES 
    ${APP_SOURCES_DIR}/nrs/nrs.cpp
    ${APP_SOURCES_DIR}/nrs/bdry/bdry.cpp
    ${APP_SOURCES_DIR}/nrs/constantFlowRate.cpp
    ${APP_SOURCES_DIR}/nrs/plugins/RANSktau.cpp
    ${APP_SOURCES_DIR}/nrs/plugins/lowMach.cpp
    ${APP_SOURCES_DIR}/nrs/postProcessing/registerPostProcessingKernels.cpp
    ${APP_SOURCES_DIR}/nrs/postProcessing/strainRotationRate.cpp
    ${APP_SOURCES_DIR}/nrs/postProcessing/aeroForces.cpp
    ${APP_SOURCES_DIR}/nrs/postProcessing/Qcriterion.cpp
)

set(APP_INCLUDE
    ${APP_SOURCES_DIR}
    ${APP_SOURCES_DIR}/nrs
    ${APP_SOURCES_DIR}/nrs/plugins
    ${APP_SOURCES_DIR}/nrs/bdry
    ${APP_SOURCES_DIR}/nrs/postProcessing
)

add_library(nekrs-lib SHARED src/lib/nekrs.cpp ${PLATFORM_SOURCES} ${CORE_SOURCES} ${SOLVER_SOURCES} ${APP_SOURCES})
if (NEKRS_BUILD_FLOAT)
  add_library(nekrs-lib-fp32 SHARED src/lib/nekrs.cpp ${PLATFORM_SOURCES} ${CORE_SOURCES} ${SOLVER_SOURCES} ${APP_SOURCES})
endif()

set_target_properties(nekrs-lib PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs)
target_link_libraries(nekrs-lib PUBLIC MPI::MPI_CXX)
if(ENABLE_CPPTRACE)
  target_link_libraries(nekrs-lib PUBLIC cpptrace::cpptrace)
  target_compile_definitions(nekrs-lib PUBLIC CPPTRACE_ENABLED)
endif()

if (NEKRS_BUILD_FLOAT)
  target_link_libraries(nekrs-lib-fp32 PUBLIC MPI::MPI_CXX)
  set_target_properties(nekrs-lib-fp32 PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs-fp32)
  if(ENABLE_CPPTRACE)
    target_link_libraries(nekrs-lib-fp32 PUBLIC cpptrace::cpptrace)
    target_compile_definitions(nekrs-lib-fp32 PUBLIC CPPTRACE_ENABLED)
  endif()
endif()

target_include_directories(nekrs-lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR} src src/lib ${APP_INCLUDE} ${SOLVER_INCLUDE} ${CORE_INCLUDE} ${PLATFORM_INCLUDE}) 
if (NEKRS_BUILD_FLOAT)
  target_include_directories(nekrs-lib-fp32 PUBLIC ${CMAKE_CURRENT_BINARY_DIR} src src/lib ${APP_INCLUDE} ${SOLVER_INCLUDE} ${CORE_INCLUDE} ${PLATFORM_INCLUDE}) 
  target_compile_definitions(nekrs-lib-fp32 PUBLIC NEKRS_USE_DFLOAT_FLOAT)
endif()

add_executable(nekrs-bin src/bin/driver.cpp)
if (NEKRS_BUILD_FLOAT)
  add_executable(nekrs-bin-fp32 src/bin/driver.cpp)
endif()

target_include_directories(nekrs-bin PRIVATE src/lib src/platform/utils)
set_target_properties(nekrs-bin PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs)
target_compile_definitions(nekrs-bin PRIVATE $<$<BOOL:${MPI_HONOR_APPNUM}>:MPI_HONOR_APPNUM=1>)
if (NEKRS_BUILD_FLOAT)
  target_include_directories(nekrs-bin-fp32 PRIVATE src/lib src/platform/utils)
  set_target_properties(nekrs-bin-fp32 PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs-fp32)
  target_compile_definitions(nekrs-bin-fp32 PRIVATE $<$<BOOL:${MPI_HONOR_APPNUM}>:MPI_HONOR_APPNUM=1>)
endif()

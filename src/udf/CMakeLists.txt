cmake_minimum_required(VERSION 3.21)

set_property(GLOBAL PROPERTY RULE_MESSAGES OFF)

set (CMAKE_CXX_STANDARD 17)

set(NEKRS_USE_DFLOAT_FLOAT OFF CACHE BOOL "use dfloat = float")

set(MPI_CXX_SKIP_MPICXX TRUE)
find_package(MPI REQUIRED)

set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")

#set(CMAKE_VERBOSE_MAKEFILE on)

set(INCLUDE_DIRS
  ${CASE_DIR}
  ${GSDIR}/gslib/include
  ${NEKRS_INSTALL_DIR}/occa/include
  ${NEKRS_INSTALL_DIR}/include
  ${NEKRS_INSTALL_DIR}/include/mesh
  ${NEKRS_INSTALL_DIR}/include/core
  ${NEKRS_INSTALL_DIR}/include/core/io
  ${NEKRS_INSTALL_DIR}/include/core/bdry
  ${NEKRS_INSTALL_DIR}/include/core/ogs
  ${NEKRS_INSTALL_DIR}/include/core/linAlg
  ${NEKRS_INSTALL_DIR}/include/utils
  ${NEKRS_INSTALL_DIR}/include/udf
  ${NEKRS_INSTALL_DIR}/include/plugins
  ${NEKRS_INSTALL_DIR}/include/nrs
  ${NEKRS_INSTALL_DIR}/include/nrs/io
  ${NEKRS_INSTALL_DIR}/include/nrs/plugins
  ${NEKRS_INSTALL_DIR}/include/nrs/postProcessing
  ${NEKRS_INSTALL_DIR}/include/nekInterface
  ${NEKRS_INSTALL_DIR}/include/nrs/neknek
  ${NEKRS_INSTALL_DIR}/include/nrs/bdry
  ${NEKRS_INSTALL_DIR}/include/setup
  ${NEKRS_INSTALL_DIR}/include/elliptic
  ${NEKRS_INSTALL_DIR}/include/nrs/cds
  ${NEKRS_INSTALL_DIR}/include/nrs/cds/cvode
  ${NEKRS_INSTALL_DIR}/include/nrs/cds/regularization
  ${NEKRS_INSTALL_DIR}/include/pointInterpolation
  ${NEKRS_INSTALL_DIR}/include/pointInterpolation/findpts
)

# dummy
add_library(OKL EXCLUDE_FROM_ALL okl.cpp)
target_include_directories(OKL PRIVATE $ENV{NEKRS_UDF_INCLUDES} ${INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(OKL PRIVATE __okl__)

add_library(udf SHARED udf.cpp)
target_link_libraries(udf PRIVATE MPI::MPI_CXX)
target_include_directories(udf PRIVATE $ENV{NEKRS_UDF_INCLUDES} ${INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})

if(cpptrace_FOUND)
  find_package(cpptrace QUIET PATHS ${NEKRS_INSTALL_DIR})
  target_link_libraries(udf PRIVATE cpptrace::cpptrace)
endif()

if(NEKRS_USE_DFLOAT_FLOAT)
  target_compile_definitions(udf PUBLIC NEKRS_USE_DFLOAT_FLOAT)
endif()

if(APPLE)
  target_link_options(udf PRIVATE "LINKER:-undefined,dynamic_lookup")
endif()

foreach(LIB $ENV{NEKRS_UDF_LIBS})
 target_link_libraries(udf PRIVATE ${LIB}) 
endforeach()

#foreach(PATH $ENV{NEKRS_UDF_RPATH})
#  target_link_options(udf PRIVATE "LINKER:-rpath,${PATH}")
#endforeach()

if(DEFINED ENV{NEKRS_ASCENT_INSTALL_DIR})
  find_package(Ascent
               PATHS $ENV{NEKRS_ASCENT_INSTALL_DIR} NO_DEFAULT_PATH)
  if (ASCENT_FOUND)
    target_link_libraries(udf PRIVATE ascent::ascent_mpi)
    target_compile_definitions(udf PUBLIC NEKRS_ENABLE_ASCENT)
  endif()
endif()

if(EXISTS "${CASE_DIR}/udf.cmake")
  message("-- Found ${CASE_DIR}/udf.cmake")
  include(${CASE_DIR}/udf.cmake)
endif()

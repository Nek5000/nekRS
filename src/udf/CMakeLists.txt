project(libUDF LANGUAGES CXX)
cmake_minimum_required(VERSION 3.11)

set_property(GLOBAL PROPERTY RULE_MESSAGES OFF)

set (CMAKE_CXX_STANDARD 17)

set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")

#set(CMAKE_VERBOSE_MAKEFILE on)

set(OGSDIR ${NEKRS_INSTALL_DIR}/gatherScatter)
set(FINDPTSDIR ${NEKRS_INSTALL_DIR}/findpts)
set(GSDIR ${NEKRS_INSTALL_DIR}/gslib)
set(INCLUDE_DIRS
  ${CASE_DIR}
  ${OGSDIR}
  ${OGSDIR}/include
  ${FINDPTSDIR}
  ${GSDIR}
  ${NEKRS_INSTALL_DIR}
  ${NEKRS_INSTALL_DIR}/include
  ${NEKRS_INSTALL_DIR}/include/mesh
  ${NEKRS_INSTALL_DIR}/include/io
  ${NEKRS_INSTALL_DIR}/include/core
  ${NEKRS_INSTALL_DIR}/include/utils
  ${NEKRS_INSTALL_DIR}/include/setup
  ${NEKRS_INSTALL_DIR}/include/udf
  ${NEKRS_INSTALL_DIR}/include/bdry
  ${NEKRS_INSTALL_DIR}/include/setup
  ${NEKRS_INSTALL_DIR}/include/navierStokes
  ${NEKRS_INSTALL_DIR}/include/solvers/elliptic
  ${NEKRS_INSTALL_DIR}/include/solvers/cvode
  ${NEKRS_INSTALL_DIR}/include/nekInterface
  ${NEKRS_INSTALL_DIR}/include/cds
  ${NEKRS_INSTALL_DIR}/include/linAlg
  ${NEKRS_INSTALL_DIR}/include/neknek
  ${NEKRS_INSTALL_DIR}/include/pointInterpolation
  ${NEKRS_INSTALL_DIR}/include/findpts
  ${NEKRS_INSTALL_DIR}/include/lpm
  ${NEKRS_INSTALL_DIR}/include/postProcessing
  ${NEKRS_INSTALL_DIR}/occa/include
)

# dummy
add_library(OKL EXCLUDE_FROM_ALL okl.cpp)
target_include_directories(OKL PRIVATE $ENV{NEKRS_UDF_INCLUDES} ${INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(OKL PRIVATE -D__okl__)

add_library(UDF SHARED udf.cpp)
target_include_directories(UDF PRIVATE $ENV{NEKRS_UDF_INCLUDES} ${INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})

if(APPLE)
    target_link_libraries(UDF "-Wl,-undefined -Wl,dynamic_lookup")
endif()
target_link_libraries(UDF "$ENV{NEKRS_UDF_LDFLAGS}")



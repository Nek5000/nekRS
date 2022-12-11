project(libUDF LANGUAGES CXX)
cmake_minimum_required(VERSION 3.11)

set_property(GLOBAL PROPERTY RULE_MESSAGES OFF)

set (CMAKE_CXX_STANDARD 17)

set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")

#set(CMAKE_VERBOSE_MAKEFILE on)

set(NEKRS_INSTALL_DIR $ENV{NEKRS_INSTALL_DIR})
set(OGSDIR      ${NEKRS_INSTALL_DIR}/gatherScatter)
set(GSDIR       ${NEKRS_INSTALL_DIR}/gslib)
set(INCLUDE_DIRS
  ${CASE_DIR}
  ${OGSDIR}
  ${GSDIR}
  ${NEKRS_INSTALL_DIR}
  ${NEKRS_INSTALL_DIR}/include
  ${NEKRS_INSTALL_DIR}/include/mesh
  ${NEKRS_INSTALL_DIR}/include/io
  ${NEKRS_INSTALL_DIR}/include/core
  ${NEKRS_INSTALL_DIR}/include/utils
  ${NEKRS_INSTALL_DIR}/include/setup
  ${NEKRS_INSTALL_DIR}/include/udf
  ${NEKRS_INSTALL_DIR}/include/setup
  ${NEKRS_INSTALL_DIR}/include/navierStokes
  ${NEKRS_INSTALL_DIR}/include/elliptic
  ${NEKRS_INSTALL_DIR}/include/nekInterface
  ${NEKRS_INSTALL_DIR}/include/cds
  ${NEKRS_INSTALL_DIR}/include/linAlg
  ${NEKRS_INSTALL_DIR}/include/postProcessing
  ${NEKRS_INSTALL_DIR}/occa/include
)

add_library(UDF SHARED udf.cpp)
target_include_directories(UDF PRIVATE $ENV{NEKRS_UDF_INCLUDES} ${INCLUDE_DIRS})

if(APPLE)
    target_link_libraries(UDF "-Wl,-undefined -Wl,dynamic_lookup")
endif()
target_link_libraries(UDF "$ENV{NEKRS_UDF_LDFLAGS}")

set_target_properties(UDF PROPERTIES LINKER_LANGUAGE CXX)

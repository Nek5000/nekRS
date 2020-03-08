project(libUDF)
cmake_minimum_required(VERSION 2.8)
enable_language(C CXX)

set (CMAKE_CXX_STANDARD 14)

set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")

set(CMAKE_VERBOSE_MAKEFILE on)

set (CMAKE_CXX_STANDARD 14)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} $ENV{NEKRS_LIBP_DEFINES}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} $ENV{NEKRS_LIBP_DEFINES}")

set(NEKRS_INSTALL_DIR $ENV{NEKRS_INSTALL_DIR})
set(HDRDIR      ${NEKRS_INSTALL_DIR}/libparanumal/include)
set(OGSDIR      ${NEKRS_INSTALL_DIR}/gatherScatter)
set(GSDIR       ${NEKRS_INSTALL_DIR}/gslib)
set(ALMONDDIR   ${NEKRS_INSTALL_DIR}/parAlmond)
set(ELLIPTICDIR ${NEKRS_INSTALL_DIR}/elliptic)
set(CDSDIR      ${NEKRS_INSTALL_DIR}/cds)
set(INSDIR      ${NEKRS_INSTALL_DIR}/ins)
set(INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}/../..
  ${INSDIR} 
  ${CDSDIR}
  ${HDRDIR}
  ${OGSDIR}
  ${GSDIR}
  ${ELLIPTICDIR}
  ${ALMONDDIR}
  ${NEKRS_INSTALL_DIR}
  ${NEKRS_INSTALL_DIR}/include
  ${NEKRS_INSTALL_DIR}/include/core
  ${NEKRS_INSTALL_DIR}/include/udf
  ${NEKRS_INSTALL_DIR}/include/nekInterface
  ${NEKRS_INSTALL_DIR}/occa/include
)

execute_process(COMMAND bash -c "cp -f ${FILENAME} ${CMAKE_CURRENT_SOURCE_DIR}/udf.cpp")
add_library(UDF SHARED udf.cpp)
target_include_directories(UDF PRIVATE ${INCLUDE_DIRS} ${UDF_DIR})

if(APPLE)
    target_link_libraries(UDF "-Wl,-undefined -Wl,dynamic_lookup")
endif()
set_target_properties(UDF PROPERTIES LINKER_LANGUAGE CXX)

message("UDF CXXFLAGS: ${CMAKE_CXX_FLAGS}")

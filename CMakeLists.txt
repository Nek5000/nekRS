cmake_minimum_required(VERSION 3.1)
project(NekRS)
enable_language(C CXX Fortran)
include(FortranCInterface)

set(CMAKE_VERBOSE_MAKEFILE on)

#Link:https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_PREFIX}/occa/lib")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_PREFIX}/elliptic")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_PREFIX}/parAlmond")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_PREFIX}/libparanumal")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_PREFIX}/gatherScatter")

###############################################################################
# Check compiler/MPI vendor, version and Fortran/C compatibility              #
###############################################################################
macro(check_compiler_id compiler)
  if(NOT "${compiler}" STREQUAL "GNU")
    message(FATAL_ERROR "NekRS only supports GNU compilers")
  endif()
endmacro()

macro(check_compiler_version compiler_version)
  if(${compiler_version} LESS 6.2)
    message(FATAL_ERROR "NekRS only supports GNU compiler versions >= 6.2")
  endif()
endmacro()

check_compiler_id("${CMAKE_CXX_COMPILER_ID}")
check_compiler_id("${CMAKE_C_COMPILER_ID}")
check_compiler_id("${CMAKE_Fortran_COMPILER_ID}")

check_compiler_version("${CMAKE_CXX_COMPILER_VERSION}")
check_compiler_version("${CMAKE_C_COMPILER_VERSION}")
check_compiler_version("${CMAKE_Fortran_COMPILER_VERSION}")

## Sanity checks for MPI
set(MPI_C_COMPILER ${CMAKE_C_COMPILER})
set(MPI_CXX_COMPILER ${CMAKE_CXX_COMPILER})
set(MPI_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})
find_package(MPI REQUIRED)
if(NOT MPI_FOUND)
  message(FATAL_ERROR "MPI is not found for one of the languages: C/CXX/Fortran")
endif()
#if (MPI_VERSION VERSION_LESS "3.1")
#  message(FATAL_ERROR "NekRS requires an MPI version >= 3.1")
#endif()

## Check Fortran/C interface
FortranCInterface_VERIFY(CXX QUIET)
if (NOT FortranCInterface_VERIFIED_CXX)
    message(WARNING "C and Fortran may not work properly")
endif()

## Use c++14
set (CMAKE_CXX_STANDARD 14)

###############################################################################
# Setup cmake and default env variables                                       #
###############################################################################
set(BIN_THIRD_PARTY_DIR "${CMAKE_CURRENT_BINARY_DIR}/3rd_party" CACHE PATH "")
set(SRC_THIRD_PARTY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/3rd_party" CACHE PATH "")

set(LIBPDIR "${BIN_THIRD_PARTY_DIR}/libparanumal" CACHE PATH "")
set(ELLIPTICDIR "${LIBPDIR}/solvers/elliptic" CACHE PATH "")
set(PARALMONDDIR "${LIBPDIR}/libs/parAlmond" CACHE PATH "")
set(OGSDIR "${LIBPDIR}/libs/gatherScatter" CACHE PATH "")
set(NEKDIR "${BIN_THIRD_PARTY_DIR}/nek5000" CACHE PATH "")
set(OCCADIR "${BIN_THIRD_PARTY_DIR}/occa" CACHE PATH "")

set(NEKINTERFACEDIR "${CMAKE_SOURCE_DIR}/src/nekInterface/" CACHE PATH "")
set(LIBP_DEFINES "$ENV{NEKRS_LIBP_DEFINES}" CACHE STRING "")
set(NEK5000_PPLIST "$ENV{NEKRS_NEK5000_PPLIST}" CACHE STRING "")

if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  set(NEKRS_DEBUG 1 CACHE BOOL "")
else()
  set(NEKRS_DEBUG 0 CACHE BOOL "")
endif()

if (DEFINED ENV{OCCA_CXX})
  set(CMAKE_OCCA_CXX_COMPILER "$ENV{OCCA_CXX}" CACHE STRING "")
else()
  set(CMAKE_OCCA_CXX_COMPILER "${CMAKE_CXX_COMPILER}" CACHE STRING "")
endif()

if (DEFINED ENV{OCCA_CXXFLAGS})
  set(CMAKE_OCCA_CXX_FLAGS "$ENV{OCCA_CXXFLAGS}" CACHE STRING "")
else()
  set(CMAKE_OCCA_CXX_FLAGS "${CMAKE_CXX_FLAGS}" CACHE STRING "")
endif()

if (DEFINED ENV{OCCA_CUDA_ENABLED})
  set(OCCA_CUDA_ENABLED "$ENV{OCCA_CUDA_ENABLED}" CACHE BOOL "")
else()
  set(OCCA_CUDA_ENABLED 1 CACHE BOOL "")
endif()

if (DEFINED ENV{OCCA_HIP_ENABLED})
  set(OCCA_HIP_ENABLED "$ENV{OCCA_HIP_ENABLED}" CACHE BOOL "")
else()
  set(OCCA_HIP_ENABLED 0 CACHE BOOL "")
endif()

if (DEFINED ENV{OCCA_OPENCL_ENABLED})
  set(OCCA_OPENCL_ENABLED "$ENV{OCCA_OPENCL_ENABLED}" CACHE BOOL "")
else()
  set(OCCA_OPENCL_ENABLED 0 CACHE BOOL "")
endif()

if (DEFINED ENV{OCCA_METAL_ENABLED})
  set(OCCA_METAL_ENABLED "$ENV{OCCA_METAL_ENABLED}" CACHE BOOL "")
else()
  set(OCCA_METAL_ENABLED 0 CACHE BOOL "")
endif()

###############################################################################
# Check if installed already                                                  #
###############################################################################
if (EXISTS "${CMAKE_INSTALL_PREFIX}/occa")
  message(WARNING "Found OCCA already installed at ${CMAKE_INSTALL_PREFIX}/.")
endif()
if (EXISTS "${CMAKE_INSTALL_PREFIX}/libparanumal")
  message(WARNING "Found libparanumal already installed at ${CMAKE_INSTALL_PREFIX}/.")
endif()

###############################################################################
# Set compiler flags                                                          #
###############################################################################
find_package(OpenMP)
if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")

if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
elseif ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}")
endif()

if (DEFINED ENV{NEKRS_LIBP_DEFINES})
  set(LIBP_DEFINES "$ENV{NEKRS_LIBP_DEFINES}" CACHE STRING "")
else()
  set(LIBP_DEFINES "-DUSE_NULL_PROJECTION=1" CACHE STRING "")
endif()

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${LIBP_DEFINES}")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LIBP_DEFINES}")

###############################################################################
# Generate nekrs.bashrc                                                       #
###############################################################################

configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/config/nekrs.conf.in
  ${CMAKE_INSTALL_PREFIX}/nekrs.conf
  @ONLY)

###############################################################################
# Start building nekRS                                                        #
###############################################################################
file(MAKE_DIRECTORY ${BIN_THIRD_PARTY_DIR})

set(NEKRS_INCLUDES)
set(THIRD_PARTY_INCLUDES)

macro(set_static_lib target_name builddir)
   set(${target_name}_lib
     "${builddir}/${CMAKE_STATIC_LIBRARY_PREFIX}${target_name}${CMAKE_STATIC_LIBRARY_SUFFIX}")
endmacro()

macro(set_shared_lib target_name builddir)
   set(${target_name}_lib
     "${builddir}/${CMAKE_SHARED_LIBRARY_PREFIX}${target_name}${CMAKE_SHARED_LIBRARY_SUFFIX}")
endmacro()

###############################################################################
# Build OCCA                                                                  #
###############################################################################

set_shared_lib(occa "${OCCADIR}/lib")
file(REMOVE ${occa_lib})
add_custom_target(occa_build DEPENDS ${occa_lib})
add_custom_command(OUTPUT ${occa_lib}
  WORKING_DIRECTORY ${BIN_THIRD_PARTY_DIR}
  COMMAND NEKRS_CC=${CMAKE_C_COMPILER}
    NEKRS_CFLAGS=${CMAKE_C_FLAGS}
    NEKRS_CXX=${CMAKE_CXX_COMPILER}
    NEKRS_CXXFLAGS=${CMAKE_CXX_FLAGS}
    NEKRS_SRC_THIRD_PARTY_DIR=${SRC_THIRD_PARTY_DIR}
    NEKRS_INSTALL_DIR=${CMAKE_INSTALL_PREFIX}
    OCCA_DIR=${OCCADIR}
    OCCA_CUDA_ENABLED="${OCCA_CUDA_ENABLED}"
    OCCA_HIP_ENABLED="${OCCA_HIP_ENABLED}"
    OCCA_OPENCL_ENABLED="${OCCA_OPENCL_ENABLED}"
    OCCA_METAL_ENABLED="${OCCA_METAL_ENABLED}"
    ${SRC_THIRD_PARTY_DIR}/occa.build
)
add_library(occa SHARED IMPORTED)
set_target_properties(occa PROPERTIES IMPORTED_LOCATION ${occa_lib})
set(THIRD_PARTY_INCLUDES ${OCCADIR}/include ${OCCADIR})

###############################################################################
# Build Nek5000                                                               #
###############################################################################
set(NEK_BUILD_LIBRARIES)
set(NEK_INCLUDES)
macro(FIND_NEK_DEPS target_name builddir include_)
  set_static_lib(${target_name} ${builddir})
  set(NEK_BUILD_LIBRARIES ${NEK_BUILD_LIBRARIES} ${${target_name}_lib})
  set(NEK_INCLUDES ${NEK_INCLUDES} ${include_})
endmacro()

FIND_NEK_DEPS(gs ${NEKDIR}/3rd_party/gslib/gslib/lib
  ${NEKDIR}/3rd_party/gslib/gslib/include)
FIND_NEK_DEPS(parRSB ${NEKDIR}/3rd_party/parRSB/parRSB/lib
  ${NEKDIR}/3rd_party/parRSB/parRSB/include)
FIND_NEK_DEPS(blasLapack ${NEKDIR}/3rd_party/blasLapack
  ${NEKDIR}/3rd_party/blasLapack)

file(REMOVE ${NEK_BUILD_LIBRARIES})
add_custom_target(nek_build DEPENDS ${NEK_BUILD_LIBRARIES})
add_custom_command(OUTPUT ${NEK_BUILD_LIBRARIES}
  WORKING_DIRECTORY ${BIN_THIRD_PARTY_DIR}
  COMMAND NEKRS_CC=${CMAKE_C_COMPILER}
    NEKRS_CFLAGS=${CMAKE_C_FLAGS}
    NEKRS_FC=${CMAKE_Fortran_COMPILER}
    NEKRS_FFLAGS=${CMAKE_Fortran_FLAGS}
    NEKRS_SRC_THIRD_PARTY_DIR=${SRC_THIRD_PARTY_DIR}
    NEKRS_NEK5000_PPLIST=${NEK5000_PPLIST}
    ${SRC_THIRD_PARTY_DIR}/nek5000.build
)

###############################################################################
# Set include dirs                                                            #
###############################################################################
set(THIRD_PARTY_INCLUDES ${THIRD_PARTY_INCLUDES}
  ${LIBPDIR}/include
  ${OGSDIR}
  ${PARALMONDDIR}
  ${ELLIPTICDIR}
  ${CDSDIR}
)

set(NEKRS_INCLUDES 
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/src/mesh
  ${CMAKE_CURRENT_SOURCE_DIR}/src/core
  ${CMAKE_CURRENT_SOURCE_DIR}/src/udf
  ${NEKINTERFACEDIR}
)

###############################################################################
# Build libP                                                                  #
###############################################################################
macro(FIND_DEPS target_name builddir)
  set_shared_lib(${target_name} ${builddir})
  add_library(${target_name} SHARED IMPORTED)
  set_target_properties(${target_name} PROPERTIES IMPORTED_LOCATION
    ${${target_name}_lib})
endmacro()

set(LIBP_BUILD_LIBRARIES)
macro(FIND_LIBP_DEPS target_name builddir)
  FIND_DEPS(${target_name} ${builddir})
  set(LIBP_BUILD_LIBRARIES ${LIBP_BUILD_LIBRARIES} ${${target_name}_lib})
endmacro()

FIND_LIBP_DEPS(P "${LIBPDIR}")
FIND_LIBP_DEPS(elliptic "${ELLIPTICDIR}")
FIND_LIBP_DEPS(parAlmond "${PARALMONDDIR}")
FIND_LIBP_DEPS(ogs "${OGSDIR}")

file(REMOVE ${LIBP_BUILD_LIBRARIES})
add_custom_target(libparanumal_build DEPENDS ${LIBP_BUILD_LIBRARIES})
add_custom_command(OUTPUT ${LIBP_BUILD_LIBRARIES}
  WORKING_DIRECTORY ${BIN_THIRD_PARTY_DIR}
  COMMAND NEKRS_CC=${CMAKE_C_COMPILER}
    NEKRS_CFLAGS=${CMAKE_C_FLAGS}
    NEKRS_CXX=${CMAKE_CXX_COMPILER}
    NEKRS_CXXFLAGS=${CMAKE_CXX_FLAGS}
    NEKRS_FC=${CMAKE_Fortran_COMPILER}
    NEKRS_FFLAGS=${CMAKE_Fortran_FLAGS}
    NEKRS_LIBP_DIR=${CMAKE_INSTALL_PREFIX}
    NEKRS_LIBP_DEFINES=${LIBP_DEFINES}
    NEKRS_SRC_THIRD_PARTY_DIR=${SRC_THIRD_PARTY_DIR}
    NEKRS_BIN_THIRD_PARTY_DIR=${BIN_THIRD_PARTY_DIR}
    NEKRS_OCCA_DIR=${BIN_THIRD_PARTY_DIR}/occa
    NEKRS_INSTALL_DIR=${CMAKE_INSTALL_PREFIX}
    ${SRC_THIRD_PARTY_DIR}/libparanumal/libparanumal.build
  DEPENDS occa_build
)

###############################################################################
# Add nekrs target                                                            #
###############################################################################
set(SRC 
    src/nekrs.cpp
    src/core/cfl.cpp
    src/core/filter.cpp
    src/core/bcMap.cpp
    src/core/insSetup.cpp
    src/core/runTime.cpp
    src/core/tombo.cpp
    src/core/cds.cpp
    src/udf/udf.cpp
    src/mesh/meshSetup.cpp
    src/mesh/meshNekReader.cpp
    src/mesh/meshPhysicalNodesHex3D.cpp
    src/mesh/meshParallelConnectNodes.cpp
    src/core/occaDeviceConfig.cpp
    src/nekInterface/nekInterfaceAdapter.cpp
    src/core/parReader.cpp
    src/core/configReader.cpp
    src/core/tinyexpr.c
    src/core/timer.cpp
    src/plugins/avg.cpp
    src/plugins/velRecycling.cpp
    src/plugins/RANSktau.cpp
    ## To get fortran flags
    src/dummy.f
)

add_library(nekrs-lib SHARED ${SRC})
set_target_properties(nekrs-lib PROPERTIES LINKER_LANGUAGE CXX)
add_dependencies(nekrs-lib occa_build libparanumal_build nek_build)
target_link_libraries(nekrs-lib ${CMAKE_DL_LIBS} elliptic parAlmond P ogs occa)
target_include_directories(nekrs-lib PRIVATE ${THIRD_PARTY_INCLUDES} ${NEKRS_INCLUDES})
set_target_properties(nekrs-lib PROPERTIES OUTPUT_NAME nekrs)

add_executable(nekrs-bin src/main.cpp)
set_target_properties(nekrs-bin PROPERTIES LINKER_LANGUAGE CXX)
add_dependencies(nekrs-bin nekrs-lib)
target_link_libraries(nekrs-bin nekrs-lib)
target_include_directories(nekrs-bin PRIVATE ${THIRD_PARTY_INCLUDES} ${NEKRS_INCLUDES})
set_target_properties(nekrs-bin PROPERTIES OUTPUT_NAME nekrs)


#################################################################################
### Install                                                                     #
#################################################################################
# Install dependencies and examples
configure_file(install_script.in.cmake install_script.cmake @ONLY)
install(SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/install_script.cmake)

#install nekRS
install(TARGETS nekrs-lib nekrs-bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include
  PRIVATE_HEADER DESTINATION include
)

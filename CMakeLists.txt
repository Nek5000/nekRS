cmake_minimum_required(VERSION 3.11)
project(NekRS LANGUAGES C CXX Fortran VERSION 20.1)
cmake_policy(SET CMP0048 NEW)   # For VERSION variables, since v3.0
include(FortranCInterface)

set(CMAKE_VERBOSE_MAKEFILE on)

#Link:https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_SKIP_BUILD_RPATH FALSE)
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
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
# Generate nekrs.conf                                                       #
###############################################################################

configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/config/nekrs.conf.in
  ${CMAKE_INSTALL_PREFIX}/nekrs.conf
  @ONLY)

###############################################################################
# Start building nekRS                                                        #
###############################################################################

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

include(config/occa.cmake)

###############################################################################
# Build Nek5000 dependencies                                                  #
###############################################################################

include(config/nek5000.cmake)

###############################################################################
# Build libP                                                                  #
###############################################################################

include(config/libparanumal.cmake)

###############################################################################
# Build NekRS
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
set_target_properties(nekrs-lib PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs)
target_link_libraries(nekrs-lib ${CMAKE_DL_LIBS} libelliptic)
target_include_directories(nekrs-lib PUBLIC 
  src
  src/mesh
  src/core
  src/udf
  ${NEKINTERFACEDIR}
  )

add_executable(nekrs-bin src/main.cpp)
set_target_properties(nekrs-bin PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs)
target_link_libraries(nekrs-bin nekrs-lib)

#################################################################################
### Install                                                                     #
#################################################################################
# Install dependencies and examples
#configure_file(install_script.in.cmake install_script.cmake @ONLY)
#install(SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/install_script.cmake)

#install nekRS
install(TARGETS nekrs-lib nekrs-bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
)

# Trailing slash prevents parent directory from being copied
install(DIRECTORY examples/ DESTINATION examples)
install(DIRECTORY scripts/ DESTINATION bin 
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
install(DIRECTORY okl/ DESTINATION okl FILES_MATCHING REGEX "\.okl$")
install(DIRECTORY src/ DESTINATION include FILES_MATCHING REGEX "\.hpp$|\.h$")
install(FILES src/udf/CMakeLists.txt DESTINATION udf)
install(DIRECTORY src/nekInterface/ DESTINATION nekInterface REGEX "\.hpp$|\.cpp$" EXCLUDE)


#################################################################################
# Make tarball
#################################################################################

set(TARNAME "nekRS-v${PROJECT_VERSION}")

add_custom_target(tar
  COMMAND ${CMAKE_COMMAND} -DTARNAME=${TARNAME} -DLIBP_SOURCE_DIR=${LIBP_SOURCE_DIR} -DHYPRE_SOURCE_DIR=${HYPRE_SOURCE_DIR} -P ${CMAKE_CURRENT_SOURCE_DIR}/config/maketar.cmake
  COMMAND tar -czf ${TARNAME}.tar.gz ${TARNAME}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})


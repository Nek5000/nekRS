cmake_minimum_required(VERSION 3.11)
project(NekRS LANGUAGES C CXX Fortran VERSION 20.1)

include(FortranCInterface)
include(ExternalProject)
include(FetchContent)

set(CMAKE_VERBOSE_MAKEFILE on)
set(FETCHCONTENT_QUIET off)
set(FETCHCONTENT_UPDATES_DISCONNECTED on)

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
# Configure Options
###############################################################################

set(NEK5000_PPLIST "PARRSB DPROCMAP" CACHE STRING 
  "Preprocessor macros for Nek5000")

set(NEKINTERFACEDIR "${CMAKE_CURRENT_SOURCE_DIR}/src/nekInterface/" CACHE PATH 
  "Directory for Nek5000 inteface source files")

set(LIBP_DEFINES "-DUSE_NULL_PROJECTION=1" CACHE STRING 
  "Preprocessor macros for libparanumal")

set(OCCA_CXX "${CMAKE_CXX_COMPILER}" CACHE STRING "C++ compiler for OCCA JIT compile")
set(OCCA_CXXFLAGS "-O2 -ftree-vectorize -funroll-loops -march=native -mtune=native"
        CACHE STRING "C++ flags for OCCA JIT compile")

set(OCCA_TAG "main" CACHE STRING "Git branch, tag, or hash for cloning OCCA")
set(NEK5000_TAG "master" CACHE STRING "Git branch, tag, or hash for cloning Nek5000")
set(HYPRE_VER "2.18.2" CACHE STRING "Version number of HYPRE to download")
set(USE_OCCA_MEM_BYTE_ALIGN "64" CACHE STRING "Memory allignment for OCCA kernels")

# From OCCA's CMakeLists. All default to ON, but are verified in config tests
#   - ENABLE_OPENMP
#   - ENABLE_CUDA
#   - ENABLE_OPENCL
#   - ENABLE_HIP
#   - ENABLE_METAL
#   - ENABLE_MPI

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
# Set compiler flags                                                          #
###############################################################################

set(CMAKE_POSITION_INDEPENDENT_CODE on)

add_definitions(-DUSE_OCCA_MEM_BYTE_ALIGN=${USE_OCCA_MEM_BYTE_ALIGN})

find_package(OpenMP)
if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# For CMake targets, CMAKE_<LANG>_FLAGS, CMAKE_<LANG>_FLAGS_<CONFIG>, and
# COMPILE_DEFINITIONS are included automatically, based on the current <CONFIG>.
# For # external projects, we explicitly collect them
function(set_external_flags lang)
  string(TOUPPER "${CMAKE_BUILD_TYPE}" mode)
  set(ext_flags "${CMAKE_${lang}_FLAGS} ${CMAKE_${lang}_FLAGS_${mode}}")

  get_directory_property(compile_defs
          DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" COMPILE_DEFINITIONS )
  foreach(def "${compile_defs}")
    set(ext_flags "${ext_flags} -D${def}")
  endforeach()

  set(EXTERNAL_${lang}_FLAGS "${ext_flags}" PARENT_SCOPE)
endfunction()

set_external_flags("C")
set_external_flags("CXX")
set_external_flags("Fortran")

###############################################################################
# Generate nekrs.conf                                                       #
###############################################################################

configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/config/nekrs.conf.in
  ${CMAKE_INSTALL_PREFIX}/nekrs.conf
  @ONLY)

###############################################################################
# Build OCCA                                                                  #
###############################################################################

include(config/occa.cmake)

###############################################################################
# Build HYPRE
###############################################################################

include(config/hypre.cmake)

###############################################################################
# Build Nek5000 dependencies                                                  #
###############################################################################

include(config/nek5000.cmake)

###############################################################################
# Definitions for libP, gslib, and blaslapack
###############################################################################

include(config/libp.cmake)
include(config/gslib.cmake)
include(config/blaslapack.cmake)

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
    src/linAlg/linAlg.cpp
    src/plugins/avg.cpp
    src/plugins/velRecycling.cpp
    src/plugins/RANSktau.cpp
    src/plugins/lowMach.cpp
    ## To get fortran flags
    src/dummy.f
    ${BLASLAPACK_SOURCES}
    ${OGS_SOURCES}
    ${LIBP_SOURCES}
    ${PARALMOND_SOURCES}
    ${ELLIPTIC_SOURCES}
)

add_library(nekrs-lib SHARED ${SRC})
set_target_properties(nekrs-lib PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs)
target_link_libraries(nekrs-lib PUBLIC libocca PRIVATE HYPRE ${GSLIB})

target_compile_definitions(nekrs-lib PUBLIC
    ${LIBP_DEFINES}
    -DDOGS="${CMAKE_INSTALL_PREFIX}/gatherScatter"
    -DDHOLMES="${CMAKE_INSTALL_PREFIX}/libparanumal" 
    -DDPARALMOND="${CMAKE_INSTALL_PREFIX}/parAlmond" 
    -DDELLIPTIC="${CMAKE_INSTALL_PREFIX}/elliptic" 
    -DHYPRE)

target_include_directories(nekrs-lib 
  PUBLIC 
  src
  src/mesh
  src/core
  src/udf
  src/linAlg
  ${NEKINTERFACEDIR}
  ${OGS_SOURCE_DIR}/include
  ${OGS_SOURCE_DIR}
  ${LIBP_SOURCE_DIR}/include
  ${ELLIPTIC_SOURCE_DIR}
  ${PARALMOND_SOURCE_DIR}/include
  ${PARALMOND_SOURCE_DIR}
  ${PARALMOND_SOURCE_DIR}/hypre
  PRIVATE
  ${HYPRE_SOURCE_DIR}/src
  ${HYPRE_SOURCE_DIR}/src/utilities
  ${HYPRE_SOURCE_DIR}/src/seq_mv
  ${HYPRE_SOURCE_DIR}/src/parcsr_mv
  ${HYPRE_SOURCE_DIR}/src/parcsr_ls
  ${HYPRE_SOURCE_DIR}/src/IJ_mv
  ${HYPRE_SOURCE_DIR}/src/multivector
  ${HYPRE_SOURCE_DIR}/src/krylov
  ${HYPRE_BINARY_DIR})

target_link_libraries(nekrs-lib PUBLIC libocca PRIVATE gs HYPRE)

add_executable(nekrs-bin src/main.cpp)
set_target_properties(nekrs-bin PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs)
target_link_libraries(nekrs-bin nekrs-lib)

set(ELLIPTIC_MAIN ${ELLIPTIC_SOURCE_DIR}/src/ellipticMain.c)
set_source_files_properties(${ELLIPTIC_MAIN} PROPERTIES LANGUAGE CXX)
add_executable(ellipticMain ${ELLIPTIC_MAIN})
set_target_properties(ellipticMain PROPERTIES LINKER_LANGUAGE CXX EXCLUDE_FROM_ALL 1)
target_link_libraries(ellipticMain nekrs-lib)

#################################################################################
### Install                                                                     #
#################################################################################

#install nekRS
install(TARGETS nekrs-lib nekrs-bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

# Trailing slash prevents parent directory from being copied
install(DIRECTORY scripts/ DESTINATION bin
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
install(DIRECTORY okl/ DESTINATION okl FILES_MATCHING REGEX "\.okl$")
install(DIRECTORY src/ DESTINATION include FILES_MATCHING REGEX "\.hpp$|\.h$")
install(FILES src/udf/CMakeLists.txt DESTINATION udf)
install(DIRECTORY src/nekInterface/ DESTINATION nekInterface REGEX "\.hpp$|\.cpp$" EXCLUDE)

# Examples will follow a different policy -- only copy if source tree is newer than install tree
file(GLOB_RECURSE EXAMPLE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/examples/*)
configure_file(${CMAKE_CURRENT_LIST_DIR}/config/install_examples.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/install_examples.cmake @ONLY)
install(SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/install_examples.cmake)

#################################################################################
# Uninstall
#################################################################################

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/config/uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/uninstall.cmake"
  IMMEDIATE @ONLY)

add_custom_target(uninstall
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/uninstall.cmake)

#################################################################################
# Make tarball
#################################################################################

set(DEST_DIR "nekRS-v${PROJECT_VERSION}")
set(TAR_NAME "nekRS-v${PROJECT_VERSION}.tar.gz")

add_custom_target(tar
  COMMAND ${CMAKE_COMMAND} 
      -DDEST_DIR=${DEST_DIR} 
      -DOCCA_SOURCE_DIR=${OCCA_SOURCE_DIR}
      -DLIBP_SOURCE_DIR=${LIBP_SOURCE_DIR}
      -DHYPRE_SOURCE_DIR=${HYPRE_SOURCE_DIR} 
      -DNEK5000_SOURCE_DIR=${NEK5000_SOURCE_DIR} 
      -P ${CMAKE_CURRENT_SOURCE_DIR}/config/maketar.cmake
  COMMAND tar -czf ${TAR_NAME} ${DEST_DIR}
  COMMAND rm -rf ${DEST_DIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})


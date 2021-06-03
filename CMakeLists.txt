cmake_minimum_required(VERSION 3.11)
project(NekRS LANGUAGES C CXX Fortran VERSION 21.1)

if(${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  MESSAGE(FATAL_ERROR "Error:  In-place builds are not supported. Please create a separate build directory")
endif(${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})

include(FortranCInterface)
include(ExternalProject)
include(FetchContent)

set_property(GLOBAL PROPERTY RULE_MESSAGES OFF)
set(CMAKE_INSTALL_MESSAGE NEVER)
set(CMAKE_VERBOSE_MAKEFILE OFF)
set(FETCHCONTENT_QUIET on)
set(FETCHCONTENT_UPDATES_DISCONNECTED on)

#Link:https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_SKIP_BUILD_RPATH FALSE)
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_PREFIX}/occa/lib")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_PREFIX}/elliptic")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_PREFIX}/amgSolver/parAlmond")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_PREFIX}/libparanumal")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_PREFIX}/gatherScatter")

###############################################################################
# Configure Options
###############################################################################

set(NEK5000_PPLIST "PARRSB DPROCMAP" CACHE STRING 
  "Preprocessor macros for Nek5000")

set(NEKINTERFACEDIR "${CMAKE_CURRENT_SOURCE_DIR}/src/nekInterface/" CACHE PATH 
  "Directory for Nek5000 inteface source files")

set(OCCA_CXX "${CMAKE_CXX_COMPILER}" CACHE STRING "C++ compiler for OCCA JIT compile")
set(OCCA_CXXFLAGS "-O2 -ftree-vectorize -funroll-loops -march=native -mtune=native"
        CACHE STRING "C++ flags for OCCA JIT compile")
set(OCCA_CUDA_COMPILER_FLAGS "-O3 --fmad=true" CACHE STRING "CUDA flags for OCCA JIT compile")
set(OCCA_HIP_COMPILER_FLAGS "-O3 -ffp-contract=fast" CACHE STRING "HIP flags for OCCA JIT compile")
set(OCCA_OPENCL_COMPILER_FLAGS "-cl-std=CL2.0 -cl-mad-enable -cl-no-signed-zeros" CACHE STRING "OPENCL flags for OCCA JIT compile")
set(GPU_MPI "0" CACHE STRING "GPU aware MPI")
set(BACKEND_DEFAULT "CUDA" CACHE STRING "Default occa mode")

# From OCCA's CMakeLists. All default to ON, but are verified in config tests
#   - ENABLE_OPENMP
#   - ENABLE_CUDA
#   - ENABLE_OPENCL
#   - ENABLE_HIP
#   - ENABLE_METAL
#   - ENABLE_MPI

set(ENABLE_AMGX OFF CACHE BOOL "Enable NVIDIA AMGX support")

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

set (CMAKE_C_STANDARD 99)
set (CMAKE_CXX_STANDARD 14)

###############################################################################
# Set compiler flags                                                          #
###############################################################################

set(CMAKE_POSITION_INDEPENDENT_CODE on)

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
  
    if ("${def}" STREQUAL "")
      set(ext_flags "${ext_flags}")
    else ()
      set(ext_flags "${ext_flags} -D${def}")
    endif()
  endforeach()

  set(EXTERNAL_${lang}_FLAGS "${ext_flags}" PARENT_SCOPE)
endfunction()

set_external_flags("C")
set_external_flags("CXX")
set_external_flags("Fortran")

execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if ("${GIT_COMMIT_HASH}" STREQUAL "")
  set(GIT_COMMIT_HASH "no sha")
endif()

###############################################################################
# Build OCCA                                                                  #
###############################################################################
include(config/occa.cmake)
set(ENABLE_CUDA ${OCCA_CUDA_ENABLED})
set(ENABLE_HIP ${OCCA_HIP_ENABLED})
set(ENABLE_OPENCL ${OCCA_OPENCL_ENABLED})

###############################################################################
# Build AMGX                                                                  #
###############################################################################

if(ENABLE_CUDA AND ENABLE_AMGX)
  add_subdirectory(3rd_party/AMGX)
endif()

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


include(config/mesh.cmake)
include(config/elliptic.cmake)
include(config/gslib.cmake)
include(config/blaslapack.cmake)

###############################################################################
# Build NekRS
###############################################################################

set(SRC 
    src/lib/nekrs.cpp
    src/io/writeFld.cpp
    src/io/utils.cpp
    src/core/utils/mysort.cpp
    src/core/utils/parallelSort.cpp
    src/core/utils/tinyexpr.c
    src/core/setupAide.cpp
    src/core/cfl.cpp
    src/core/filter.cpp
    src/core/bcMap.cpp
    src/core/setup.cpp
    src/timeStepper/timeStepper.cpp
    src/lns/tombo.cpp
    src/cds/cds.cpp
    src/core/parReader.cpp
    src/core/configReader.cpp
    src/core/timer.cpp
    src/core/platform.cpp
    src/linAlg/linAlg.cpp
    src/linAlg/matrixConditionNumber.cpp
    src/linAlg/matrixInverse.cpp
    src/linAlg/matrixEig.cpp
    src/linAlg/matrixTranspose.cpp
    src/linAlg/matrixRightSolve.cpp
    src/plugins/avg.cpp
    src/plugins/velRecycling.cpp
    src/plugins/RANSktau.cpp
    src/plugins/lowMach.cpp
    src/udf/udf.cpp
    src/nekInterface/nekInterfaceAdapter.cpp
    ${MESH_SOURCES}
    ${PARALMOND_SOURCES}
    ${ELLIPTIC_SOURCES}
    ${BLASLAPACK_SOURCES}
    ${OGS_SOURCES}
)

set_property(
   SOURCE src/lib/nekrs.cpp 
   APPEND PROPERTY COMPILE_DEFINITIONS
   GITCOMMITHASH="${GIT_COMMIT_HASH}"
   NEKRS_VERSION=${PROJECT_VERSION_MAJOR}
   NEKRS_SUBVERSION=${PROJECT_VERSION_MINOR}
)

add_library(nekrs-lib SHARED ${SRC})
set_target_properties(nekrs-lib PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs)
target_link_libraries(nekrs-lib PUBLIC libocca PRIVATE HYPRE ${GSLIB})

target_compile_definitions(nekrs-lib PUBLIC
    -DDOGS="${CMAKE_INSTALL_PREFIX}/gatherScatter"
    -DHYPRE)

target_include_directories(nekrs-lib 
  PUBLIC 
  ${CMAKE_CURRENT_BINARY_DIR}
  src/core
  src/core/utils
  src/lib
  src/io
  src/udf
  src/linAlg
  src/timeStepper
  src/lns
  src/cds
  ${MESH_SOURCE_DIR}
  ${NEKINTERFACEDIR}
  ${OGS_SOURCE_DIR}/include
  ${OGS_SOURCE_DIR}
  ${ELLIPTIC_SOURCE_DIR}
  PRIVATE
  ${ELLIPTIC_SOURCE_DIR}/amgSolver/hypre
  ${ELLIPTIC_SOURCE_DIR}/amgSolver/amgx
  ${PARALMOND_SOURCE_DIR}
  ${HYPRE_SOURCE_DIR}/src
  ${HYPRE_SOURCE_DIR}/src/utilities
  ${HYPRE_SOURCE_DIR}/src/seq_mv
  ${HYPRE_SOURCE_DIR}/src/parcsr_mv
  ${HYPRE_SOURCE_DIR}/src/parcsr_ls
  ${HYPRE_SOURCE_DIR}/src/IJ_mv
  ${HYPRE_SOURCE_DIR}/src/multivector
  ${HYPRE_SOURCE_DIR}/src/krylov
  ${HYPRE_BINARY_DIR})

if(ENABLE_AMGX AND ENABLE_CUDA)
  target_compile_definitions(nekrs-lib PUBLIC -DENABLE_AMGX)
  target_link_libraries(nekrs-lib PUBLIC amgxsh)
  target_include_directories(nekrs-lib PUBLIC 3rd_party/AMGX/base/include)
endif()
target_link_libraries(nekrs-lib PUBLIC libocca PRIVATE gs HYPRE)

add_executable(nekrs-bin src/main.cpp)
target_include_directories(nekrs-bin PRIVATE src/lib)
set_target_properties(nekrs-bin PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs)
target_link_libraries(nekrs-bin nekrs-lib)

message(NOTICE "")
message(NOTICE "------------ Summary ------------")
message(NOTICE "Installation directory: ${CMAKE_CURRENT_SOURCE_DIR}") 
message(NOTICE "C compiler: ${CMAKE_C_COMPILER}")
message(NOTICE "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(NOTICE "Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message(NOTICE "Default backend : ${BACKEND_DEFAULT}")
message(NOTICE "CPU backend compiler: ${OCCA_CXX} (flags: ${OCCA_CXXFLAGS})")
if(ENABLE_CUDA) 
  message(NOTICE "NVIDIA CUDA backend enabled (flags: ${OCCA_CUDA_COMPILER_FLAGS})")
endif()
if(ENABLE_HIP) 
  message(NOTICE "AMD HIP backend enabled (flags: ${OCCA_HIP_COMPILER_FLAGS})")
endif()
if(ENABLE_OPENCL) 
  message(NOTICE "OPENCL backend enabled (flags: ${OCCA_OPENCL_COMPILER_FLAGS})")
endif()
if(GPU_MPI)
  message(NOTICE "GPU aware MPI support enabled")
else(GPU_MPI)
  message(NOTICE "GPU aware MPI support disabled!")
endif()
message(NOTICE "")

###############################################################################
# Generate nekrs.conf                                                       #
###############################################################################

configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/config/nekrs.conf.in
  ${CMAKE_INSTALL_PREFIX}/nekrs.conf
  @ONLY)

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
install(DIRECTORY okl/ DESTINATION okl FILES_MATCHING REGEX "\.okl$|\.c$")
install(DIRECTORY src/ DESTINATION include FILES_MATCHING REGEX "\.hpp$|\.h$|\.tpp$")
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
      -DHYPRE_SOURCE_DIR=${HYPRE_SOURCE_DIR} 
      -DNEK5000_SOURCE_DIR=${NEK5000_SOURCE_DIR} 
      -P ${CMAKE_CURRENT_SOURCE_DIR}/config/maketar.cmake
  COMMAND tar -czf ${TAR_NAME} ${DEST_DIR}
  COMMAND rm -rf ${DEST_DIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})




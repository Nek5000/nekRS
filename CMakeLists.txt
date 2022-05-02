cmake_minimum_required(VERSION 3.13)
project(NekRS LANGUAGES C CXX Fortran VERSION 22.0.0)

if(${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  MESSAGE(FATAL_ERROR "Error:  In-place builds are not supported. Please create a separate build directory")
endif(${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})

include(CheckCCompilerFlag)
include(FortranCInterface)
include(ExternalProject)
include(FetchContent)

include(config/utils.cmake)

#set_property(GLOBAL PROPERTY RULE_MESSAGES OFF)
set(CMAKE_INSTALL_MESSAGE NEVER)
#set(CMAKE_VERBOSE_MAKEFILE ON)
set(FETCHCONTENT_QUIET on)
set(FETCHCONTENT_UPDATES_DISCONNECTED on)

#Link:https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_SKIP_BUILD_RPATH FALSE)
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

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

# From OCCA's CMakeLists. All default to ON, but are verified in config tests
#   - ENABLE_OPENMP
#   - ENABLE_CUDA
#   - ENABLE_OPENCL
#   - ENABLE_HIP
#   - ENABLE_METAL
#   - ENABLE_MPI

set(ENABLE_METAL OFF CACHE BOOL "Enable OCCA Metal support")
set(ENABLE_AMGX OFF CACHE BOOL "Enable NVIDIA AMGX support")
set(GPU_MPI "0" CACHE STRING "GPU aware MPI")

###############################################################################
# Check compiler/MPI vendor, version and Fortran/C compatibility              #
###############################################################################

macro(check_fcompiler_id compiler)
  if(NOT "${compiler}" STREQUAL "GNU")
  message(FATAL_ERROR "NekRS only supports gfortran")
  endif()
endmacro()

check_fcompiler_id("${CMAKE_Fortran_COMPILER_ID}")

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

find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  if(OpenMP_CXX_INCLUDE_DIR)
    include_directories(SYSTEM ${OpenMP_CXX_INCLUDE_DIRS})
  endif()
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

# avoid hijacked symbols with the same name
set(BSYMBOLIC_FLAG " -Wl,-Bsymbolic -Wl,-Bsymbolic-functions")
list(APPEND CMAKE_REQUIRED_LIBRARIES ${BSYMBOLIC_FLAG})
check_c_compiler_flag("" Allowed_LD_Flag_BSYMBOLIC)
if(NOT Allowed_LD_Flag_BSYMBOLIC OR APPLE)
  set(BSYMBOLIC_FLAG "")
else()
#  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${BSYMBOLIC_FLAG}")
endif()
list(REMOVE_ITEM CMAKE_REQUIRED_LIBRARIES ${BSYMBOLIC_FLAG})

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


include(config/bench.cmake)
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
    src/utils/inipp.cpp
    src/utils/mysort.cpp
    src/utils/parallelSort.cpp
    src/utils/tinyexpr.c
    src/utils/setupAide.cpp
    src/core/numberActiveFields.cpp
    src/core/printHeader.cpp
    src/timeStepper/cfl.cpp
    src/regularization/filter.cpp
    src/regularization/avm.cpp
    src/bdry/bcMap.cpp
    src/core/compileKernels.cpp
    src/setup/setup.cpp
    src/core/alignment.cpp
    src/core/registerNrsKernels.cpp
    src/bdry/createEToBV.cpp
    src/bdry/applyDirichlet.cpp
    src/timeStepper/timeStepper.cpp
    src/timeStepper/subCycling.cpp
    src/navierStokes/tombo.cpp
    src/navierStokes/constantFlowRate.cpp
    src/cds/cdsSolve.cpp
    src/cds/registerCdsKernels.cpp
    src/io/parReader.cpp
    src/io/re2Reader.cpp
    src/core/configReader.cpp
    src/core/timer.cpp
    src/core/platform.cpp
    src/core/flopCounter.cpp
    src/core/kernelRequestManager.cpp
    src/core/device.cpp
    src/linAlg/linAlg.cpp
    src/linAlg/matrixConditionNumber.cpp
    src/linAlg/matrixInverse.cpp
    src/linAlg/matrixEig.cpp
    src/linAlg/matrixTranspose.cpp
    src/linAlg/matrixRightSolve.cpp
    src/linAlg/registerLinAlgKernels.cpp
    src/plugins/tavg.cpp
    src/plugins/velRecycling.cpp
    src/plugins/RANSktau.cpp
    src/plugins/lowMach.cpp
    src/udf/udf.cpp
    src/udf/compileUDFKernels.cpp
    src/nekInterface/nekInterfaceAdapter.cpp
    src/postProcessing/planarAvg.cpp
    ${BENCH_SOURCES}
    ${MESH_SOURCES}
    ${PARALMOND_SOURCES}
    ${ELLIPTIC_SOURCES}
    ${BLASLAPACK_SOURCES}
    ${OGS_SOURCES}
)

set_property(
   SOURCE src/core/printHeader.cpp 
   APPEND PROPERTY COMPILE_DEFINITIONS
   GITCOMMITHASH="${GIT_COMMIT_HASH}"
   NEKRS_VERSION=${PROJECT_VERSION_MAJOR}
   NEKRS_SUBVERSION=${PROJECT_VERSION_MINOR}
   NEKRS_PATCHVERSION=${PROJECT_VERSION_PATCH}
)

add_library(nekrs-lib SHARED ${SRC})

add_dependencies(nekrs-lib HYPRE_BUILD)
#target_link_libraries(nekrs-lib PUBLIC libocca OpenMP::OpenMP_CXX PRIVATE gs ${GSLIB})
target_link_libraries(nekrs-lib PUBLIC libocca OpenMP::OpenMP_CXX PRIVATE gs ${GSLIB})

target_compile_definitions(nekrs-lib PUBLIC
    -DDOGS="${CMAKE_INSTALL_PREFIX}/gatherScatter"
    -DHYPRE)

target_include_directories(nekrs-lib 
  PUBLIC 
  ${CMAKE_CURRENT_BINARY_DIR}
  src
  src/setup
  src/bdry
  src/core
  src/utils
  src/lib
  src/io
  src/udf
  src/regularization
  src/linAlg
  src/timeStepper
  src/navierStokes
  src/cds
  src/postProcessing
  ${BENCH_SOURCE_DIR}
  ${BENCH_SOURCE_DIR}/core
  ${BENCH_SOURCE_DIR}/fdm
  ${BENCH_SOURCE_DIR}/axHelm
  ${BENCH_SOURCE_DIR}/advsub
  ${MESH_SOURCE_DIR}
  ${NEKINTERFACEDIR}
  ${OGS_SOURCE_DIR}/include
  ${OGS_SOURCE_DIR}
  ${ELLIPTIC_SOURCE_DIR}
  PRIVATE
  ${ELLIPTIC_SOURCE_DIR}/amgSolver/hypre
  ${ELLIPTIC_SOURCE_DIR}/amgSolver/amgx
  ${PARALMOND_SOURCE_DIR}
  ${HYPRE_INCLUDE_DIR}
)

if(ENABLE_AMGX AND ENABLE_CUDA)
  target_compile_definitions(nekrs-lib PUBLIC -DENABLE_AMGX)
  target_link_libraries(nekrs-lib PUBLIC amgxsh amgx)
  target_include_directories(nekrs-lib PUBLIC 3rd_party/AMGX/base/include 3rd_party/AMGX/core/include)
endif()

add_executable(nekrs-bin src/main.cpp)
target_include_directories(nekrs-bin PRIVATE src/lib)
set_target_properties(nekrs-bin PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs)
target_link_libraries(nekrs-bin PRIVATE nekrs-lib)

add_executable(axhelm-bin src/bench/axHelm/main.cpp)
set_target_properties(axhelm-bin PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs-bench-axhelm)
target_link_libraries(axhelm-bin PRIVATE nekrs-lib)

add_executable(advsub-bin src/bench/advsub/main.cpp)
set_target_properties(advsub-bin PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs-bench-advsub)
target_link_libraries(advsub-bin PRIVATE nekrs-lib)

add_executable(fdm-bin src/bench/fdm/main.cpp)
set_target_properties(fdm-bin PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs-bench-fdm)
target_link_libraries(fdm-bin PRIVATE nekrs-lib)

if(ENABLE_CUDA)
  set(BACKEND_DEFAULT "CUDA" CACHE STRING "Default occa mode")
elseif(ENABLE_HIP)
  set(BACKEND_DEFAULT "HIP" CACHE STRING "Default occa mode")
else()
  set(BACKEND_DEFAULT "SERIAL" CACHE STRING "Default occa mode")
endif() 

#################################################################################
### Install                                                                     #
#################################################################################

# Generate nekrs.conf
__MPI_underlying_compiler("CXX" MPI_UNDERLYING_COMPILER)

configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/config/nekrs.conf.in
  ${CMAKE_INSTALL_PREFIX}/nekrs.conf
  @ONLY)

#install nekRS
install(TARGETS nekrs-lib nekrs-bin axhelm-bin advsub-bin fdm-bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

# Trailing slash prevents parent directory from being copied
install(DIRECTORY scripts/ DESTINATION bin
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
install(DIRECTORY okl/ DESTINATION okl FILES_MATCHING REGEX "\.okl$|\.c|\.h$")
install(DIRECTORY src/ DESTINATION include FILES_MATCHING REGEX "\.hpp$|\.h$|\.tpp$")
install(FILES src/udf/CMakeLists.txt DESTINATION udf)
install(DIRECTORY src/nekInterface/ DESTINATION nekInterface REGEX "\.hpp$|\.cpp$" EXCLUDE)

# Examples will follow a different policy -- only copy if source tree is newer than install tree
file(GLOB_RECURSE EXAMPLE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/examples/*)
configure_file(${CMAKE_CURRENT_LIST_DIR}/config/install_examples.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/install_examples.cmake @ONLY)
install(SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/install_examples.cmake)

install(FILES src/io/parHelp.txt DESTINATION include)

message(NOTICE "")
message(NOTICE "----------------- Summary -----------------")
message(NOTICE "Installation directory: ${CMAKE_INSTALL_PREFIX}") 
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
message(NOTICE "-------------------------------------------")

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

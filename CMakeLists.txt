cmake_minimum_required(VERSION 3.18)
project(NekRS LANGUAGES C CXX Fortran VERSION 23.0.1)

if(${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  MESSAGE(FATAL_ERROR "In-place builds are not supported. Please create a separate build directory")
endif(${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})

include(CheckCCompilerFlag)
include(CheckLinkerFlag)
include(FortranCInterface)
include(ExternalProject)
include(FetchContent)
include(config/utils.cmake)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#set(CMAKE_INSTALL_MESSAGE NEVER)
set(CMAKE_INSTALL_MESSAGE ALWAYS)
set(FETCHCONTENT_QUIET on)
set(FETCHCONTENT_UPDATES_DISCONNECTED on)
set(CMAKE_ENABLE_EXPORTS TRUE)

if (NOT CMAKE_CONFIGURATION_TYPES AND 
    NOT CMAKE_NO_BUILD_TYPE AND
    NOT CMAKE_BUILD_TYPE AND
    CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

#Link:https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_SKIP_BUILD_RPATH FALSE)
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_PREFIX}/occa/lib")

###############################################################################
# Check compiler/MPI vendor, version and Fortran/C compatibility              #
###############################################################################

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
set (CMAKE_CXX_STANDARD 17)

__MPI_underlying_compiler("CXX" MPI_UNDERLYING_COMPILER)


###############################################################################
# Configure Options
###############################################################################

set(NEK5000_PPLIST "PARRSB DPROCMAP" CACHE STRING 
  "Preprocessor macros for Nek5000")
set(NEKINTERFACEDIR "${CMAKE_CURRENT_SOURCE_DIR}/src/nekInterface/" CACHE PATH 
"Directory for Nek5000 inteface source files")

string(COMPARE EQUAL "${CMAKE_CXX_COMPILER_ID}" "GNU" USING_GNU)
string(COMPARE EQUAL "${CMAKE_CXX_COMPILER_ID}" "IntelLLVM" USING_INTEL_LLVM)
string(COMPARE EQUAL "${CMAKE_CXX_COMPILER_ID}" "Clang" USING_LLVM)
string(COMPARE EQUAL "${CMAKE_CXX_COMPILER_ID}" "ARMClang" USING_LLVM)
string(COMPARE EQUAL "${CMAKE_CXX_COMPILER_ID}" "AppleClang" USING_APPLE_LLVM)

if(USING_GNU)
  set(NEKRS_COMPILER_FLAGS "-O2 -g -march=native -mtune=native -ftree-vectorize" CACHE STRING "compiler flags")
else()
  set(NEKRS_COMPILER_FLAGS "-O2 -g" CACHE STRING "compiler flags")
endif()

set(CMAKE_CXX_FLAGS_RELWITHDEBINFO ${NEKRS_COMPILER_FLAGS})
set(CMAKE_C_FLAGS_RELWITHDEBINFO ${NEKRS_COMPILER_FLAGS})
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO ${NEKRS_COMPILER_FLAGS})

set(OCCA_CXX "${MPI_UNDERLYING_COMPILER}" CACHE STRING "C++ compiler for OCCA JIT compile")
if(USING_INTEL_LLVM)
  set(OCCA_CXXFLAGS "-w -O3 -g -xHost" CACHE STRING "C++ flags for OCCA JIT compile")
elseif(USING_GNU)
  set(OCCA_CXXFLAGS "-w -O3 -g -march=native -mtune=native" CACHE STRING "C++ flags for OCCA JIT compile")
elseif(USING_LLVM)
  set(OCCA_CXXFLAGS "-w -O3 -g -march=native -mtune=native" CACHE STRING "C++ flags for OCCA JIT compile")
else()
  set(OCCA_CXXFLAGS "-w -O3 -g" CACHE STRING "C++ flags for OCCA JIT compile")
endif()

set(OCCA_CUDA_COMPILER_FLAGS "-w -O3 -lineinfo --use_fast_math" CACHE STRING "CUDA flags for OCCA JIT compile")
set(OCCA_HIP_COMPILER_FLAGS "-w -O3 -g -ffp-contract=fast -funsafe-math-optimizations -ffast-math -munsafe-fp-atomics" CACHE STRING "HIP flags for OCCA JIT compile")
set(OCCA_DPCPP_COMPILER_FLAGS "-w -O3 -fsycl -fdebug-info-for-profiling -gline-tables-only -ffp-contract=fast" CACHE STRING "DPC++ flags for OCCA JIT compile")
set(OCCA_OPENCL_COMPILER_FLAGS "-w -cl-std=CL2.0 -cl-unsafe-math-optimizations -cl-fast-relaxed-math -cl-mad-enable -cl-no-signed-zeros" CACHE STRING "OPENCL flags for OCCA JIT compile")

set(OCCA_FOUND ON)

if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(powerpc|ppc)64le")
  string(REPLACE "-march" "-mcpu" OCCA_CXXFLAGS ${OCCA_CXXFLAGS})
  string(REPLACE "-march" "-mcpu" CMAKE_Fortran_FLAGS_RELWITHDEBINFO ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO})
  string(REPLACE "-march" "-mcpu" CMAKE_C_FLAGS_RELWITHDEBINFO ${CMAKE_C_FLAGS_RELWITHDEBINFO})
  string(REPLACE "-march" "-mcpu" CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
endif()

set(ENABLE_CUDA ON CACHE BOOL "Enable OCCA CUDA support")
set(ENABLE_HIP ON CACHE BOOL "Enable OCCA HIP support")
set(ENABLE_DPCPP ON CACHE BOOL "Enable OCCA SYCL support")
set(ENABLE_OPENCL OFF CACHE BOOL "Enable OCCA OpenCL support")
set(ENABLE_METAL OFF CACHE BOOL "Enable OCCA Metal support")

set(ENABLE_AMGX OFF CACHE BOOL "Enable NVIDIA AMGX support")
set(NEKRS_GPU_MPI ON CACHE BOOL "Enable GPU aware MPI")
set(ENABLE_CVODE OFF CACHE BOOL "Enable CVODE support")

###############################################################################
# Set compiler flags                                                          #
###############################################################################

set(CMAKE_POSITION_INDEPENDENT_CODE on)

if(NOT USING_APPLE_LLVM)
  find_package(OpenMP REQUIRED)
endif()

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

set(BSYMBOLIC_FLAG "LINKER:-Bsymbolic LINKER:-Bsymbolic-functions")
check_linker_flag(C ${BSYMBOLIC_FLAG} Allowed_LD_Flag_BSYMBOLIC)
if(NOT Allowed_LD_Flag_BSYMBOLIC OR APPLE)
  set(BSYMBOLIC_FLAG "")
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1)
    message(FATAL_ERROR "GNU compiler version must be at least 9.1!")
  endif()
endif()

###############################################################################
# Build OCCA                                                                  #
###############################################################################
include(config/occa.cmake)
message("-- OCCA_CUDA_ENABLED: ${OCCA_CUDA_ENABLED}")
message("-- OCCA_HIP_ENABLED: ${OCCA_HIP_ENABLED}")
message("-- OCCA_DPCPP_ENABLED: ${OCCA_DPCPP_ENABLED}")
message("-- OCCA_OPENCL_ENABLED: ${OCCA_OPENCL_ENABLED}")

set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})

###############################################################################
# Build AMGX                                                                  #
###############################################################################
if(OCCA_CUDA_ENABLED AND ENABLE_AMGX)
  add_subdirectory(3rd_party/AMGX)
endif()

###############################################################################
# CVODE
###############################################################################
if(ENABLE_CVODE)
  # workaround to resolve name conflict (used be OCCA and CVODE) 
  set(ENABLE_CUDA ${OCCA_CUDA_ENABLED})
  set(ENABLE_HIP ${OCCA_HIP_ENABLED})
  set(ENABLE_OPENCL ${OCCA_OPENCL_ENABLED})

  if(OCCA_HIP_ENABLED)
    message(FATAL_ERROR "CVODE + HIP is not currently supported!")
  endif()
  set(ENABLE_MPI ON)
  set(BUILD_CVODE ON)
  set(BUILD_ARKODE OFF)
  set(BUILD_CVODES OFF)
  set(BUILD_IDA OFF)
  set(BUILD_IDAS OFF)
  set(BUILD_KINSOL OFF)
  set(BUILD_EXAMPLES OFF)
  set(BUILD_TESTING OFF)
  set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
  add_subdirectory(3rd_party/cvode)
endif()

###############################################################################
# Build HYPRE
###############################################################################
if(OCCA_CUDA_ENABLED)
  set(ENABLE_HYPRE_GPU OFF CACHE BOOL "Enable HYPRE GPU support")
endif()
include(config/hypre.cmake)

###############################################################################
# Build BLAS/LAPACK
###############################################################################
add_subdirectory(3rd_party/lapack)

###############################################################################
# Build Nek5000 dependencies                                                  #
###############################################################################
include(config/nek5000.cmake)

###############################################################################
# Build NekRS
###############################################################################
include(config/bench.cmake)
include(config/mesh.cmake)
include(config/elliptic.cmake)
include(config/gslib.cmake)

set(SRC 
    src/lib/nekrs.cpp
    src/io/writeFld.cpp
    src/io/fileUtils.cpp
    src/utils/inipp.cpp
    src/utils/unifdef.c
    src/utils/mysort.cpp
    src/utils/parallelSort.cpp
    src/utils/tinyexpr.c
    src/utils/setupAide.cpp
    src/core/numberActiveFields.cpp
    src/core/printHeader.cpp
    src/navierStokes/cfl.cpp
    src/regularization/filter.cpp
    src/regularization/avm.cpp
    src/bdry/bcMap.cpp
    src/core/compileKernels.cpp
    src/setup/setup.cpp
    src/bdry/alignment.cpp
    src/core/registerNrsKernels.cpp
    src/core/registerNekNekKernels.cpp
    src/core/registerPostProcessingKernels.cpp
    src/core/registerEllipticKernels.cpp
    src/core/registerEllipticPreconditionerKernels.cpp
    src/core/registerCdsKernels.cpp
    src/core/registerLinAlgKernels.cpp
    src/core/registerMeshKernels.cpp
    src/bdry/createEToBV.cpp
    src/navierStokes/applyDirichlet.cpp
    src/navierStokes/timeStepper.cpp
    src/navierStokes/subCycling.cpp
    src/navierStokes/tombo.cpp
    src/navierStokes/constantFlowRate.cpp
    src/navierStokes/Urst.cpp
    src/cds/cdsSolve.cpp
    src/cds/subCycling.cpp
    src/setup/parReader.cpp
    src/io/re2Reader.cpp
    src/setup/configReader.cpp
    src/core/timer.cpp
    src/core/platform.cpp
    src/core/comm.cpp
    src/core/flopCounter.cpp
    src/core/kernelRequestManager.cpp
    src/core/device.cpp
    src/linAlg/linAlg.cpp
    src/linAlg/matrixConditionNumber.cpp
    src/linAlg/matrixInverse.cpp
    src/linAlg/matrixEig.cpp
    src/linAlg/matrixTranspose.cpp
    src/linAlg/matrixRightSolve.cpp
    src/plugins/tavg.cpp
    src/plugins/velRecycling.cpp
    src/plugins/RANSktau.cpp
    src/plugins/lowMach.cpp
    src/plugins/lpm.cpp
    src/findpts/findpts.cpp
    src/pointInterpolation/pointInterpolation.cpp
    src/neknek/neknek.cpp
    src/udf/udf.cpp
    src/udf/compileUDFKernels.cpp
    src/nekInterface/nekInterfaceAdapter.cpp
    src/postProcessing/planarAvg.cpp
    src/postProcessing/strainRotationRate.cpp
    src/postProcessing/viscousDrag.cpp
    src/postProcessing/Qcriterion.cpp
    src/core/registerCvodeKernels.cpp
    src/solvers/cvode/cvode.cpp
    ${BENCH_SOURCES}
    ${MESH_SOURCES}
    ${ELLIPTIC_SOURCES}
    ${OGS_SOURCES}
    ${FINDPTS_SOURCES}
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

set_target_properties(nekrs-lib PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs)
target_link_libraries(nekrs-lib PUBLIC libocca PRIVATE gs ${GSLIB} blas lapack)
if(OpenMP_FOUND)
target_link_libraries(nekrs-lib PUBLIC OpenMP::OpenMP_CXX)
endif()

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
  src/navierStokes
  src/neknek
  src/cds
  src/findpts
  src/postProcessing
  src/pointInterpolation
  src/solvers/cvode
  src/lns
  ${BENCH_SOURCE_DIR}
  ${BENCH_SOURCE_DIR}/core
  ${BENCH_SOURCE_DIR}/fdm
  ${BENCH_SOURCE_DIR}/axHelm
  ${BENCH_SOURCE_DIR}/advsub
  ${MESH_SOURCE_DIR}
  ${NEKINTERFACEDIR}
  ${OGS_SOURCE_DIR}/include
  ${OGS_SOURCE_DIR}
  ${FINDPTS_SOURCE_DIR}
  ${ELLIPTIC_SOURCE_DIR}
  PRIVATE
  ${ELLIPTIC_SOURCE_DIR}/amgSolver/hypre
  ${ELLIPTIC_SOURCE_DIR}/amgSolver/amgx
  ${ELLIPTIC_SOURCE_DIR}/MG
)

if(ENABLE_AMGX AND OCCA_CUDA_ENABLED)
  target_compile_definitions(nekrs-lib PUBLIC -DENABLE_AMGX)
  target_link_libraries(nekrs-lib PUBLIC amgxsh amgx)
  target_include_directories(nekrs-lib PUBLIC 3rd_party/AMGX/base/include 3rd_party/AMGX/core/include)
endif()

target_link_libraries(nekrs-lib PRIVATE nekrs-hypre)
target_link_libraries(nekrs-lib PRIVATE nekrs-hypre-device)

if (ENABLE_CVODE)
  target_link_libraries (nekrs-lib PUBLIC SUNDIALS::cvode_shared SUNDIALS::nvecmpiplusx_shared SUNDIALS::nvecserial_shared)
  target_include_directories(nekrs-lib PUBLIC 3rd_party/cvode/include)
  target_compile_definitions(nekrs-lib PUBLIC -DENABLE_CVODE)
  if (OCCA_CUDA_ENABLED)
    find_package(CUDA REQUIRED)
    target_include_directories(nekrs-lib PRIVATE ${CUDA_INCLUDE_DIRS})
    target_link_libraries (nekrs-lib PUBLIC SUNDIALS::nveccuda_shared)
    target_compile_definitions(nekrs-lib PUBLIC -DENABLE_CUDA)
  endif()
endif()

add_executable(nekrs-bin src/main.cpp)
target_include_directories(nekrs-bin PRIVATE src/lib)
set_target_properties(nekrs-bin PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs)
target_link_libraries(nekrs-bin PRIVATE nekrs-lib ${CMAKE_DL_LIBS})

###############################################################################
# Optional nekRS plugins
###############################################################################
file(GLOB plugins LIST_DIRECTORIES true ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/plugins/*)
set(PLUGINS_RPATH "")
foreach(DIR ${plugins})
  get_filename_component(PLUGIN_NAME ${DIR} NAME)
  IF(IS_DIRECTORY ${DIR})
      message(NOTICE "-- Found plugin: ${PLUGIN_NAME}")
      add_subdirectory(${DIR})
  ELSE()
      CONTINUE()
  ENDIF()
  set(UDF_LIBS "${UDF_LIBS} ${CMAKE_INSTALL_PREFIX}/3rdParty/${PLUGIN_NAME}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}${PLUGIN_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}")
  set(PLUGINS_RPATH "${PLUGINS_RPATH};${CMAKE_INSTALL_PREFIX}/3rdParty/${PLUGIN_NAME}/lib")
  set(UDF_INCLUDES "${UDF_INCLUDES};${CMAKE_INSTALL_PREFIX}/3rdParty/${PLUGIN_NAME}/include")
endforeach()

set_target_properties(nekrs-bin PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${PLUGINS_RPATH}")

add_executable(axhelm-bin src/bench/axHelm/main.cpp)
set_target_properties(axhelm-bin PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs-bench-axhelm)
target_link_libraries(axhelm-bin PRIVATE nekrs-lib)

add_executable(advsub-bin src/bench/advsub/main.cpp)
set_target_properties(advsub-bin PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs-bench-advsub)
target_link_libraries(advsub-bin PRIVATE nekrs-lib)

add_executable(fdm-bin src/bench/fdm/main.cpp)
set_target_properties(fdm-bin PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME nekrs-bench-fdm)
target_link_libraries(fdm-bin PRIVATE nekrs-lib)

if(OCCA_CUDA_ENABLED)
  set(BACKEND_DEFAULT "CUDA" CACHE STRING "Default occa mode")
elseif(OCCA_HIP_ENABLED)
  set(BACKEND_DEFAULT "HIP" CACHE STRING "Default occa mode")
elseif(OCCA_DPCPP_ENABLED)
  set(BACKEND_DEFAULT "DPCPP" CACHE STRING "Default occa mode")
elseif(OCCA_OPENCL_ENABLED)
  set(BACKEND_DEFAULT "OPENCL" CACHE STRING "Default occa mode")
else()
  set(BACKEND_DEFAULT "SERIAL" CACHE STRING "Default occa mode")
  set(NEKRS_GPU_MPI OFF CACHE BOOL "Enable GPU aware MPI" FORCE)
endif() 

#################################################################################
### Install                                                                     #
#################################################################################

if(NEKRS_GPU_MPI)
  set(GPU_MPI "1")
else(NEKRS_GPU_MPI)
  set(GPU_MPI "0")
endif()

configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/config/nekrs.conf.in
  ${CMAKE_INSTALL_PREFIX}/nekrs.conf
  @ONLY)

install(
  TARGETS nekrs-lib nekrs-hypre nekrs-hypre-device nekrs-bin axhelm-bin advsub-bin fdm-bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

# Trailing slash prevents parent directory from being copied
install(DIRECTORY scripts/ DESTINATION bin
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
install(DIRECTORY kernels/ DESTINATION kernels FILES_MATCHING REGEX "\.okl$|\.c|\.h$|\.hpp$|")
install(DIRECTORY src/ DESTINATION include FILES_MATCHING REGEX "\.hpp$|\.h$|\.tpp$")
install(FILES src/udf/CMakeLists.txt DESTINATION udf)
install(DIRECTORY src/nekInterface/ DESTINATION nekInterface REGEX "\.hpp$|\.cpp$" EXCLUDE)

# Examples will follow a different policy -- only copy if source tree is newer than install tree
file(GLOB_RECURSE EXAMPLE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/examples/*)
configure_file(${CMAKE_CURRENT_LIST_DIR}/config/install_examples.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/install_examples.cmake @ONLY)
install(SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/install_examples.cmake)

install(FILES doc/parHelp.txt DESTINATION doc)

message(NOTICE "")
message(NOTICE "----------------- Summary -----------------")
message(NOTICE "Installation directory: ${CMAKE_INSTALL_PREFIX}") 
message(NOTICE "C compiler: ${CMAKE_C_COMPILER}")
message(NOTICE "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(NOTICE "Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message(NOTICE "Default backend : ${BACKEND_DEFAULT}")
message(NOTICE "CPU backend compiler: ${OCCA_CXX} (flags: ${OCCA_CXXFLAGS})")
if(OCCA_CUDA_ENABLED) 
  message(NOTICE "NVIDIA CUDA backend enabled (flags: ${OCCA_CUDA_COMPILER_FLAGS})")
endif()
if(OCCA_HIP_ENABLED) 
  message(NOTICE "AMD HIP backend enabled (flags: ${OCCA_HIP_COMPILER_FLAGS})")
endif()
if(OCCA_DPCPP_ENABLED)
  message(NOTICE "DPCPP backend enabled (flags: ${OCCA_DPCPP_COMPILER_FLAGS})")
endif()
if(OCCA_OPENCL_ENABLED) 
  message(NOTICE "OPENCL backend enabled (flags: ${OCCA_OPENCL_COMPILER_FLAGS})")
endif()
message(NOTICE "GPU aware MPI support: ${NEKRS_GPU_MPI}")
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

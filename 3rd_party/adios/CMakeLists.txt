#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.14)

# Fail immediately if not using an out-of-source build
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
  message(FATAL_ERROR
    "In-source builds are not supported.  Please create a build directory "
    "separate from the source directory")
endif()

include(${CMAKE_CURRENT_LIST_DIR}/cmake/ADIOSFunctions.cmake)

setup_version(2.10.1)

project(ADIOS2 VERSION ${ADIOS2_VERSION})

#------------------------------------------------------------------------------#
# Some boilerplate to setup nice output directories
#------------------------------------------------------------------------------#
include(GNUInstallDirs)

if (ADIOS2_USE_PIP)
  # Only UNIX is supported currently
  if (UNIX)
    # Linux bundles what libraries we have when they're put beside the modules.
    set(CMAKE_INSTALL_LIBDIR "adios2")
    set(CMAKE_INSTALL_INCLUDEDIR "adios2/include")
    # ELF loader settings.
    set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)
    list(APPEND CMAKE_INSTALL_RPATH "$ORIGIN")
  endif ()
endif ()

set(CMAKE_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/adios2
  CACHE STRING "Installation CMake subdirectory")
mark_as_advanced(CMAKE_INSTALL_CMAKEDIR)

list(INSERT CMAKE_MODULE_PATH 0 "${PROJECT_SOURCE_DIR}/cmake")
if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
endif()
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
endif()
if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
endif()

# Let windows builds auto-export dll symbols
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#------------------------------------------------------------------------------#
# Silence MSVC warnings
#------------------------------------------------------------------------------#

if(MSVC)
  add_definitions(
    -D_CRT_SECURE_NO_DEPRECATE
    -D_CRT_SECURE_NO_WARNINGS
    -D_SCL_SECURE_NO_DEPRECATE
    -D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING
  )
endif()

function(append_str_list_once var x)
  foreach(item IN ITEMS ${${var}})
    if(item STREQUAL x)
      return()
    endif()
  endforeach()
  set(out ${${var}} ${x})
  string(REPLACE ";" " " out "${out}")
  set(${var} "${out}" PARENT_SCOPE)
endfunction()

if(NOT MSVC)
  set(CMAKE_C_COMPILER_FLAGS_DEBUG "${CMAKE_C_COMPILER_FLAGS_DEBUG}"
    CACHE STRING "" FORCE
  )
  set(CMAKE_CXX_COMPILER_FLAGS_DEBUG "${CMAKE_CXX_COMPILER_FLAGS_DEBUG}"
    CACHE STRING "" FORCE
  )
endif()

#------------------------------------------------------------------------------#
# Deal with any pre-installation cleanup tasks
#------------------------------------------------------------------------------#
add_subdirectory(cmake/install/pre)

#------------------------------------------------------------------------------#
# Top level options
#------------------------------------------------------------------------------#

# Default to a debug build if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

# Force C++11 (No extensions) and C99
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT MSVC)
  set(CMAKE_C_STANDARD 99)
  set(CMAKE_C_STANDARD_REQUIRED True)
endif()

# Option for packagers to add a suffix after "adios" in library names (to
# install multiple builds in parallel)
set(ADIOS2_LIBRARY_SUFFIX "" CACHE STRING
  "Suffix to put after \"adios\" in library names")
set(ADIOS2_EXECUTABLE_SUFFIX "${ADIOS2_LIBRARY_SUFFIX}" CACHE STRING
  "Suffix to append to executable names")
mark_as_advanced(ADIOS2_LIBRARY_SUFFIX ADIOS2_EXECUTABLE_SUFFIX)

# Use meta-compile features if available, otherwise use specific language
# features
if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Intel|Clang|AppleClang|MSVC)$")
  set(ADIOS2_CXX11_FEATURES cxx_auto_type cxx_nullptr)
else()
  set(ADIOS2_CXX11_FEATURES cxx_std_11)
endif()
if(CMAKE_C_COMPILER_ID MATCHES "^(GNU|Intel|Clang|AppleClang|MSVC)$")
  set(ADIOS2_C99_FEATURES c_restrict)
else()
  set(ADIOS2_C99_FEATURES c_std_99)
endif()

include(CMakeDependentOption)

# Setup shared library defaults.  If explicitly specified somehow, then default
# to that.  Otherwise base the default on whether or not shared libs are even
# supported.
get_property(SHARED_LIBS_SUPPORTED GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS)
cmake_dependent_option(BUILD_SHARED_LIBS
  "Build shared libraries (so/dylib/dll)." ${SHARED_LIBS_SUPPORTED}
  "SHARED_LIBS_SUPPORTED" OFF
)
mark_as_advanced(BUILD_SHARED_LIBS)
if((NOT BUILD_SHARED_LIBS) AND (NOT DEFINED CMAKE_POSITION_INDEPENDENT_CODE))
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Ctest creates BUILD_TESTING option and sets it to true by default
# Here we disable BUILD_TESTING option by default
option(BUILD_TESTING "Build testing" OFF)
include(CTest)
mark_as_advanced(CLEAR BUILD_TESTING)

adios_option(Blosc2     "Enable support for c-blosc-2 transforms" AUTO)
adios_option(BZip2      "Enable support for BZip2 transforms" AUTO)
adios_option(ZFP        "Enable support for ZFP transforms" AUTO)
adios_option(SZ         "Enable support for SZ transforms" AUTO)
adios_option(LIBPRESSIO "Enable support for LIBPRESSIO transforms" AUTO)
adios_option(MGARD      "Enable support for MGARD transforms" AUTO)
adios_option(PNG        "Enable support for PNG transforms" AUTO)
adios_option(CUDA       "Enable support for Cuda" OFF)
adios_option(Kokkos     "Enable support for Kokkos" OFF)
adios_option(MPI        "Enable support for MPI" AUTO)
adios_option(DAOS       "Enable support for DAOS" AUTO)
adios_option(DataMan    "Enable support for DataMan" AUTO)
adios_option(DataSpaces "Enable support for DATASPACES" AUTO)
adios_option(MHS        "Enable support for MHS" AUTO)
adios_option(SST        "Enable support for SST" AUTO)
adios_option(ZeroMQ     "Enable support for ZeroMQ" AUTO)
adios_option(HDF5       "Enable support for the HDF5 engine" AUTO)
adios_option(HDF5_VOL   "Enable support for HDF5 ADIOS2 VOL" AUTO)
adios_option(IME        "Enable support for DDN IME transport" AUTO)
adios_option(Python     "Enable support for Python bindings" AUTO)
adios_option(Fortran    "Enable support for Fortran bindings" AUTO)
adios_option(SysVShMem  "Enable support for SysV Shared Memory IPC on *NIX" AUTO)
adios_option(UCX        "Enable support for UCX DataPlane in SST" AUTO)
adios_option(Profiling  "Enable support for profiling" AUTO)
adios_option(Endian_Reverse "Enable support for Little/Big Endian Interoperability" AUTO)
adios_option(Sodium     "Enable support for Sodium for encryption" AUTO)
adios_option(Catalyst   "Enable support for in situ visualization plugin using ParaView Catalyst" AUTO)
adios_option(Campaign   "Enable support for Campaigns (requires SQLite3 and ZLIB)" AUTO)
adios_option(AWSSDK     "Enable support for S3 compatible storage using AWS SDK's S3 module" OFF)
adios_option(Derived_Variable    "Enable support for derived variables" OFF)
adios_option(PIP        "Enable support for pip packaging" OFF)

option(ADIOS2_Blosc2_PREFER_SHARED "Prefer shared Blosc2 libraries" ON)
mark_as_advanced(ADIOS2_Blosc2_PREFER_SHARED)
mark_as_advanced(ADIOS2_USE_PIP)
include(${PROJECT_SOURCE_DIR}/cmake/DetectOptions.cmake)

if(ADIOS2_HAVE_CUDA OR ADIOS2_HAVE_Kokkos_CUDA)
    set(CMAKE_CUDA_STANDARD 11)
    set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        if(DEFINED Kokkos_CUDA_ARCHITECTURES)
            set(CMAKE_CUDA_ARCHITECTURES ${Kokkos_CUDA_ARCHITECTURES})
        else()
            # Minimum common non-deprecated architecture
            set(CMAKE_CUDA_ARCHITECTURES 52)
        endif()
    endif()
endif()

if(NOT DEFINED CMAKE_HIP_ARCHITECTURES AND DEFINED Kokkos_HIP_ARCHITECTURES)
    set(CMAKE_HIP_ARCHITECTURES ${Kokkos_HIP_ARCHITECTURES})
endif()

if(ADIOS2_HAVE_MPI)
  if(MPIEXEC_MAX_NUMPROCS LESS 4 AND "$ENV{OMPI_MCA_rmaps_base_oversubscribe}")
    message(STATUS "OpenMPI oversubscribe detected: raising MPIEXEC_MAX_NUMPROCS to 4 for testing")
    set(MPIEXEC_MAX_NUMPROCS 4 CACHE STRING "" FORCE)
  endif()

  # Workaround for MPI forcing the link of C++ bindings
  add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)
endif()

#------------------------------------------------------------------------------#
# POSIX O_DIRECT is only working for Unix in adios for now
#------------------------------------------------------------------------------#
if(CYGWIN)
  #-tb O_DIRECT messes up cygwin
  set(ADIOS2_HAVE_O_DIRECT 0)
elseif(MSVC)
  # Windows has other things but we are not using them
  set(ADIOS2_HAVE_O_DIRECT 0)
elseif(APPLE)
  # Mac has other things but we are not using them
  set(ADIOS2_HAVE_O_DIRECT 0)
else()

  message(STATUS "Checking for O_DIRECT")
  include(CheckCXXSourceCompiles)
  check_cxx_source_compiles("
#include <unistd.h>
#include <fcntl.h>
int main(int argc, char * argv[]) { argc = O_DIRECT; }
" O_DIRECT_WORKS)

  if (O_DIRECT_WORKS)
    set(ADIOS2_HAVE_O_DIRECT 1)
  else()
    set(ADIOS2_HAVE_O_DIRECT 0)
  endif()

endif()

#if(NOT HAVE_O_DIRECT)
#  message(WARNING " -----  The open() flag O_DIRECT is not available! ---- ")
#else()
#  message(STATUS " -----  The open() flag O_DIRECT is available! ---- ")
#endif()


set(ADIOS2_CONFIG_OPTS
    DataMan DataSpaces HDF5 HDF5_VOL MHS SST Fortran MPI Python PIP Blosc2 BZip2
    LIBPRESSIO MGARD MGARD_MDR PNG SZ ZFP DAOS IME O_DIRECT Sodium Catalyst SysVShMem UCX
    ZeroMQ Profiling Endian_Reverse Derived_Variable AWSSDK GPU_Support CUDA Kokkos
    Kokkos_CUDA Kokkos_HIP Kokkos_SYCL Campaign
)

GenerateADIOSHeaderConfig(${ADIOS2_CONFIG_OPTS})
configure_file(
  ${PROJECT_SOURCE_DIR}/CTestCustom.cmake.in
  ${PROJECT_BINARY_DIR}/CTestCustom.cmake
  @ONLY
)

install(FILES ${PROJECT_BINARY_DIR}/source/adios2/common/ADIOSConfig.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/adios2/common
  COMPONENT adios2_core-development
)

#------------------------------------------------------------------------------#
# Some specialized shared library and RPATH handling
#------------------------------------------------------------------------------#
option(ADIOS2_RUN_INSTALL_TEST "Enable / disable the install test" TRUE)
mark_as_advanced(ADIOS2_RUN_INSTALL_TEST)
if(BUILD_SHARED_LIBS AND ADIOS2_RUN_INSTALL_TEST)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  string(REGEX REPLACE "[^/]+" ".." relative_base "${CMAKE_INSTALL_LIBDIR}")

  if(APPLE)
    # Various homebrew MPI fortran installations break rpath usage on OSX
    if(APPLE AND ADIOS2_HAVE_Fortran AND ADIOS2_HAVE_MPI)
      set(CMAKE_MACOSX_RPATH OFF)
      if(NOT CMAKE_VERSION VERSION_LESS 3.17)
        set(CMAKE_INSTALL_NAME_DIR
          "$<INSTALL_PREFIX>/${CMAKE_INSTALL_LIBDIR}")
      else()
        set(CMAKE_INSTALL_NAME_DIR
          "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
        message(STATUS "Skipping install tests due to platform specific rpath oddities")
        set(ADIOS2_RUN_INSTALL_TEST FALSE)
      endif()
    else()
      list(APPEND CMAKE_INSTALL_RPATH "@loader_path/${relative_base}/${CMAKE_INSTALL_LIBDIR}")
    endif()

  elseif(CMAKE_SYSTEM_NAME MATCHES "Linux")
    # Linux needs some specialized RPATH handling
    list(APPEND CMAKE_INSTALL_RPATH "$ORIGIN/${relative_base}/${CMAKE_INSTALL_LIBDIR}")
  endif()
endif()

if(MSVC AND BUILD_SHARED_LIBS AND ADIOS2_HAVE_HDF5)
  # See note about building with visual studio and shared libs, here:
  #     https://github.com/HDFGroup/hdf5/blob/develop/release_docs/USING_HDF5_VS.txt
  add_definitions(-DH5_BUILT_AS_DYNAMIC_LIB=1)
endif()

#------------------------------------------------------------------------------#
# Third party libraries
#------------------------------------------------------------------------------#
add_subdirectory(thirdparty)

#------------------------------------------------------------------------------#
# Main library source
#------------------------------------------------------------------------------#
add_subdirectory(source)

#------------------------------------------------------------------------------#
# Language bindings
#------------------------------------------------------------------------------#
add_subdirectory(bindings)

#------------------------------------------------------------------------------#
# Language Libraries
#------------------------------------------------------------------------------#
if(ADIOS2_HAVE_Python)
  add_subdirectory(python)
endif()

#------------------------------------------------------------------------------#
# Plugins
#------------------------------------------------------------------------------#
add_subdirectory(plugins)

#------------------------------------------------------------------------------#
# Examples
#------------------------------------------------------------------------------#
option(ADIOS2_BUILD_EXAMPLES "Build examples" OFF)
if(ADIOS2_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()


#------------------------------------------------------------------------------#
# Testing
#------------------------------------------------------------------------------#
cmake_dependent_option(ADIOS2_RUN_MPI_MPMD_TESTS
  "Enable the parallel MPMD tests" ON
  "BUILD_TESTING;ADIOS2_HAVE_MPI" OFF)
mark_as_advanced(ADIOS2_RUN_MPI_MPMD_TESTS)

if(BUILD_TESTING)
  find_program(DIFF_EXECUTABLE diff)
  if(DIFF_EXECUTABLE)
    set(DIFF_COMMAND ${DIFF_EXECUTABLE})
  elseif(GIT_COMMAND)
    set(DIFF_COMMAND ${GIT_COMMAND} diff)
  elseif(GITCOMMAND)
    set(DIFF_COMMAND ${GITCOMMAND} diff)
  endif()
  if(NOT DIFF_COMMAND)
    message(FATAL_ERROR "diff command is required for testing")
  endif()

  enable_testing()
  add_subdirectory(testing)
endif()

#------------------------------------------------------------------------------#
# Generating package configs
#------------------------------------------------------------------------------#
add_subdirectory(cmake/install/packaging)

#------------------------------------------------------------------------------#
# Generating static adios-config
#------------------------------------------------------------------------------#
option(ADIOS2_INSTALL_GENERATE_CONFIG "Install a generated adios2-config shell script for non-cmake projects" ON)
mark_as_advanced(ADIOS2_INSTALL_GENERATE_CONFIG)
if(ADIOS2_INSTALL_GENERATE_CONFIG)
  add_subdirectory(cmake/install/post)
endif()

#------------------------------------------------------------------------------#
# Configuration summary
#------------------------------------------------------------------------------#
message("")
message("ADIOS2 build configuration:")
message("  ADIOS Version: ${ADIOS2_VERSION}")
message("  C++ Compiler : ${CMAKE_CXX_COMPILER_ID} "
                         "${CMAKE_CXX_COMPILER_VERSION} "
                         "${CMAKE_CXX_COMPILER_WRAPPER}")
message("    ${CMAKE_CXX_COMPILER}")
message("")
if(ADIOS2_HAVE_CUDA)
	message("  Cuda Compiler : ${CMAKE_CUDA_COMPILER} ")
	message("    Cuda Architecture : ${CMAKE_CUDA_ARCHITECTURES} ")
endif()
if(ADIOS2_HAVE_Fortran)
  message("  Fortran Compiler : ${CMAKE_Fortran_COMPILER_ID} "
                               "${CMAKE_Fortran_COMPILER_VERSION} "
                               "${CMAKE_Fortran_COMPILER_WRAPPER}")
  message("    ${CMAKE_Fortran_COMPILER}")
  message("")
endif()
message("  Installation prefix: ${CMAKE_INSTALL_PREFIX}")
message("        bin: ${CMAKE_INSTALL_BINDIR}")
message("        lib: ${CMAKE_INSTALL_LIBDIR}")
message("    include: ${CMAKE_INSTALL_INCLUDEDIR}")
message("      cmake: ${CMAKE_INSTALL_CMAKEDIR}")
if(DEFINED CMAKE_INSTALL_PYTHONDIR)
  message("     python: ${CMAKE_INSTALL_PYTHONDIR}")
endif()
message("")
message("  Features:")
if(BUILD_SHARED_LIBS)
  set(msg_lib_type "shared")
elseif(CMAKE_POSITION_INDEPENDENT_CODE)
  set(msg_lib_type "static (with PIC)")
else()
  set(msg_lib_type "static (without PIC)")
endif()
message("    Library Type: ${msg_lib_type}")
message("    Build Type:   ${CMAKE_BUILD_TYPE}")
message("    Testing: ${BUILD_TESTING}")
message("    Examples: ${ADIOS2_BUILD_EXAMPLES}")
message("    Build Options:")

foreach(opt IN LISTS ADIOS2_CONFIG_OPTS)
  message_pad("      ${opt}" 25 label)
  if(${ADIOS2_HAVE_${opt}})
    message("${label}: ON")
  else()
    message("${label}: OFF")
  endif()
endforeach()

set (HPCDataPlaneList "")
if (ADIOS2_HAVE_SST)
  if (ADIOS2_SST_HAVE_LIBFABRIC)
     set (HPCDataPlaneList "${HPCDataPlaneList} fabric")
  endif()
  if (ADIOS2_SST_HAVE_UCX)
     set (HPCDataPlaneList "${HPCDataPlaneList} UCX")
  endif()
  if (ADIOS2_SST_HAVE_MPI_DP)
     set (HPCDataPlaneList "${HPCDataPlaneList} MPI")
  endif()
endif()
if ({HPCDataPlaneList} STREQUAL  "")
  message("    Possible RDMA DataPlanes for SST: <none>")
else()
  message("    Possible RDMA DataPlanes for SST: ${HPCDataPlaneList}")
endif()


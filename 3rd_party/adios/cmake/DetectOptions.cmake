#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

# This file contains the option and dependency logic.  The configuration
# options are designed to be tertiary: ON, OFF, or AUTO.  If AUTO, we try to
# determine if dependencies are available and enable the option if we find
# them, otherwise we disable it.  If explicitly ON then a failure to find
# dependencies is an error,

# Helper function to extract a common prefix from two strings
function(string_get_prefix in0 in1 outVar)
  string(LENGTH "${in0}" len0)
  string(LENGTH "${in1}" len1)
  set(lenMax ${len0})
  if(len1 LESS len0)
    set(lenMax ${len1})
  endif()
  set(lenPfx 0)
  if(lenMax GREATER 0)
    foreach(len RANGE 1 ${lenMax})
      string(SUBSTRING "${in0}" 0 ${len} sub0)
      string(SUBSTRING "${in1}" 0 ${len} sub1)
      if(NOT (sub0 STREQUAL sub1))
        break()
      endif()
      set(lenPfx ${len})
    endforeach()
  endif()
  string(SUBSTRING "${in0}" 0 ${lenPfx} outTmp)
  set(${outVar} "${outTmp}" PARENT_SCOPE)
endfunction()

# Helper function to strip a common prefix off an input string
function(string_strip_prefix pfx in outVar)
  string(LENGTH "${pfx}" lenPfx)
  string(LENGTH "${in}" lenIn)
  if(lenPfx GREATER lenIn)
    set(${outVar} "" PARENT_SCOPE)
    return()
  endif()
  string(SUBSTRING "${in}" 0 ${lenPfx} inPfx)
  if(NOT pfx STREQUAL inPfx)
    set(${outVar} "" PARENT_SCOPE)
    return()
  endif()
  string(SUBSTRING "${in}" ${lenPfx} -1 outTmp)
  set(${outVar} "${outTmp}" PARENT_SCOPE)
endfunction()

# Extract the common prefix from a collection of variables
function(lists_get_prefix listVars outVar)
  foreach(l IN LISTS listVars)
    foreach(d IN LISTS ${l})
      if(NOT prefix)
        set(prefix "${d}")
        continue()
      endif()
      string_get_prefix("${prefix}" "${d}" prefix)
      if(NOT prefix)
        set(${outVar} "" PARENT_SCOPE)
        return()
      endif()
    endforeach()
  endforeach()
  set(${outVar} "${prefix}" PARENT_SCOPE)
endfunction()

# Multithreading
find_package(Threads REQUIRED)

# Blosc2
if(ADIOS2_USE_Blosc2 STREQUAL AUTO)
  find_package(Blosc2 2.10.1 QUIET)
elseif(ADIOS2_USE_Blosc2)
  find_package(Blosc2 2.10.1)
endif()
if(Blosc2_FOUND)
  set(ADIOS2_HAVE_Blosc2 TRUE)
  if(TARGET Blosc2::blosc2_shared)
    set(blosc2_shlib_available ON)
  endif()

  if(TARGET Blosc2::blosc2_static)
    set(blosc2_slib_available ON)
  endif()

  if (blosc2_shlib_available AND (NOT blosc2_slib_available OR ADIOS2_Blosc2_PREFER_SHARED))
    set(adios2_blosc2_tgt Blosc2::blosc2_shared)
  elseif(blosc2_slib_available)
    set(adios2_blosc2_tgt Blosc2::blosc2_static)
  else()
    message(FATAL_ERROR "Blosc2 cmake package found but no targets exists inside it.")
  endif()

  add_library(adios2_blosc2 ALIAS ${adios2_blosc2_tgt})
endif()

# BZip2
if(ADIOS2_USE_BZip2 STREQUAL AUTO)
  find_package(BZip2)
elseif(ADIOS2_USE_BZip2)
  find_package(BZip2 REQUIRED)
endif()
if(BZIP2_FOUND)
  set(ADIOS2_HAVE_BZip2 TRUE)
endif()

# ZFP
if(ADIOS2_USE_ZFP)
  find_package(ZFP 1.0.0 CONFIG QUIET)
  if(NOT ZFP_FOUND)
    if(ADIOS2_USE_ZFP STREQUAL AUTO)
      find_package(ZFP 0.5.3 CONFIG)
    else()
      find_package(ZFP 0.5.3 REQUIRED CONFIG)
    endif()
  endif()
endif()
if(ZFP_FOUND)
  set(ADIOS2_HAVE_ZFP TRUE)
  set(ADIOS2_HAVE_ZFP_CUDA ${ZFP_WITH_CUDA})

  # Older versions of ZFP
  if(NOT ADIOS2_HAVE_ZFP_CUDA)
    get_target_property(ZFP_INCLUDE_DIRECTORIES zfp::zfp INTERFACE_INCLUDE_DIRECTORIES)
    set(CMAKE_REQUIRED_INCLUDES ${ZFP_INCLUDE_DIRECTORIES})
    set(CMAKE_REQUIRED_LIBRARIES zfp::zfp)
    include(CheckCSourceRuns)
    check_c_source_runs("
    #include <zfp.h>

    int main()
    {
      zfp_stream* stream = zfp_stream_open(NULL);
      return !zfp_stream_set_execution(stream, zfp_exec_cuda);
    }"
    ADIOS2_HAVE_ZFP_CUDA)
    unset(CMAKE_REQUIRED_INCLUDES)
    unset(CMAKE_REQUIRED_LIBRARIES)
  endif()

  if(ADIOS2_HAVE_ZFP_CUDA)
    add_compile_definitions(ADIOS2_HAVE_ZFP_CUDA)
  endif()
endif()

# SZ
if(ADIOS2_USE_SZ STREQUAL AUTO)
  find_package(SZ)
elseif(ADIOS2_USE_SZ)
  find_package(SZ REQUIRED)
endif()
if(SZ_FOUND)
  set(ADIOS2_HAVE_SZ TRUE)
endif()

# LIBPRESSIO
if(ADIOS2_USE_LIBPRESSIO STREQUAL AUTO)
 find_package(LibPressio QUIET CONFIG)
 find_package(std_compat QUIET CONFIG)
elseif(ADIOS2_USE_LIBPRESSIO)
  find_package(LibPressio REQUIRED CONFIG)
  find_package(std_compat REQUIRED CONFIG)
endif()
if(TARGET LibPressio::libpressio)
  set(ADIOS2_HAVE_LIBPRESSIO TRUE)
endif()

# MGARD
set(ADIOS2_HAVE_MGARD_MDR FALSE)
if(ADIOS2_USE_MGARD STREQUAL AUTO)
  find_package(mgard CONFIG)
elseif(ADIOS2_USE_MGARD)
  find_package(mgard REQUIRED CONFIG)
endif()
if(mgard_FOUND)
  set(ADIOS2_HAVE_MGARD TRUE)
  if(MGARD_ENABLE_MDR)
     set(ADIOS2_HAVE_MGARD_MDR TRUE)
 endif()
endif()

# PNG
if(ADIOS2_USE_PNG STREQUAL AUTO)
  find_package(PNG 1.6.0)
elseif(ADIOS2_USE_PNG)
  find_package(PNG 1.6.0 REQUIRED)
endif()
if(PNG_FOUND)
  set(ADIOS2_HAVE_PNG TRUE)
endif()

set(mpi_find_components C)

if(ADIOS2_USE_Derived_Variable)
    set(ADIOS2_HAVE_Derived_Variable TRUE)
endif()

if(ADIOS2_USE_Kokkos AND ADIOS2_USE_CUDA)
  message(FATAL_ERROR "ADIOS2_USE_Kokkos is incompatible with ADIOS2_USE_CUDA")
endif()

# Kokkos
if(ADIOS2_USE_Kokkos)
  if(ADIOS2_USE_Kokkos STREQUAL AUTO)
    find_package(Kokkos 3.7 QUIET)
  else()
    find_package(Kokkos 3.7 REQUIRED)
  endif()
  if(Kokkos_FOUND)
    set(ADIOS2_HAVE_Kokkos TRUE)
    if(Kokkos_ENABLE_CUDA OR Kokkos_ENABLE_HIP OR Kokkos_ENABLE_SYCL)
      if(Kokkos_ENABLE_CUDA)
        set(ADIOS2_HAVE_Kokkos_CUDA TRUE)
        enable_language(CUDA)
      endif()
      if(Kokkos_ENABLE_HIP)
        set(ADIOS2_HAVE_Kokkos_HIP TRUE)
        enable_language(HIP)
      endif()
      if(Kokkos_ENABLE_SYCL)
        set(ADIOS2_HAVE_Kokkos_SYCL TRUE)
      endif()
      set(ADIOS2_HAVE_GPU_Support TRUE)
    endif()
  endif()
endif()

# CUDA
if(ADIOS2_USE_CUDA)
  include(CheckLanguage)
  check_language(CUDA)
  if(ADIOS2_USE_CUDA STREQUAL AUTO)
    find_package(CUDAToolkit QUIET)
  else()
    find_package(CUDAToolkit REQUIRED)
  endif()
  if(CMAKE_CUDA_COMPILER AND CUDAToolkit_FOUND)
    enable_language(CUDA)
    set(ADIOS2_HAVE_CUDA TRUE)
    set(ADIOS2_HAVE_GPU_Support TRUE)
  endif()
endif()

if(ADIOS2_HAVE_Kokkos AND ADIOS2_HAVE_CUDA)
    message(FATAL_ERROR "The Kokkos and CUDA backends cannot be active concurrently")
endif()

# Fortran
if(ADIOS2_USE_Fortran STREQUAL AUTO)
  include(CheckLanguage)
  check_language(Fortran)
  if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
  endif()
elseif(ADIOS2_USE_Fortran)
  enable_language(Fortran)
  if(CMAKE_Fortran_COMPILER_LOADED)
  set(ADIOS2_HAVE_Fortran TRUE)
  list(APPEND mpi_find_components Fortran)

  include(CheckFortranSourceCompiles)
  check_fortran_source_compiles("
    program testargs
      integer :: n
      character(len=256) :: v
      n = command_argument_count()
      call get_command_argument(0, v)
    end program testargs"
    ADIOS2_HAVE_FORTRAN_F03_ARGS
    SRC_EXT F90
  )
  check_fortran_source_compiles("
    program testargs
      integer :: n
      character(len=256) :: v
      n = iargc()
      call getarg(0, v)
    end program testargs"
    ADIOS2_HAVE_FORTRAN_GNU_ARGS
    SRC_EXT F90
  )
 endif()
endif()

# MPI
list(APPEND mpi_find_components OPTIONAL_COMPONENTS CXX)
if(ADIOS2_USE_MPI STREQUAL AUTO)
  find_package(MPI COMPONENTS ${mpi_find_components})
elseif(ADIOS2_USE_MPI)
  find_package(MPI COMPONENTS ${mpi_find_components} REQUIRED)
endif()
if(MPI_FOUND)
  set(ADIOS2_HAVE_MPI TRUE)
endif()

# ZeroMQ
if(ADIOS2_USE_ZeroMQ STREQUAL AUTO)
    find_package(ZeroMQ 4.1)
elseif(ADIOS2_USE_ZeroMQ)
    find_package(ZeroMQ 4.1 REQUIRED)
endif()
if(ZeroMQ_FOUND)
    set(ADIOS2_HAVE_ZeroMQ TRUE)
endif()

# DataMan
if(ZeroMQ_FOUND)
    if(ADIOS2_USE_DataMan STREQUAL AUTO)
        set(ADIOS2_HAVE_DataMan TRUE)
    elseif(ADIOS2_USE_DataMan)
        set(ADIOS2_HAVE_DataMan TRUE)
    endif()
endif()

# MHS
if(ADIOS2_USE_MHS STREQUAL AUTO)
    set(ADIOS2_HAVE_MHS TRUE)
elseif(ADIOS2_USE_MHS)
    set(ADIOS2_HAVE_MHS TRUE)
endif()

# DataSpaces
if(MPI_FOUND)
    if(ADIOS2_USE_DataSpaces STREQUAL AUTO)
        find_package(DataSpaces 2.1.1)
    elseif(ADIOS2_USE_DataSpaces)
        find_package(DataSpaces 2.1.1 REQUIRED)
    endif()
    if(DATASPACES_FOUND)
        set(ADIOS2_HAVE_DataSpaces TRUE)
    endif()
endif()

# HDF5
if(ADIOS2_USE_HDF5 STREQUAL AUTO)
  find_package(HDF5 COMPONENTS C)
  if(HDF5_FOUND AND (NOT HDF5_IS_PARALLEL OR ADIOS2_HAVE_MPI))
    set(ADIOS2_HAVE_HDF5 TRUE)
  endif()
elseif(ADIOS2_USE_HDF5)
  find_package(HDF5 REQUIRED COMPONENTS C)
  if(HDF5_IS_PARALLEL AND NOT ADIOS2_HAVE_MPI)
    message(FATAL_ERROR "MPI is disabled but parallel HDF5 is detected.")
  endif()
  set(ADIOS2_HAVE_HDF5 TRUE)

  # HDF5 VOL requires 1.13+
  if(ADIOS2_USE_HDF5_VOL AND ADIOS2_HAVE_HDF5)
    if(HDF5_VERSION GREATER_EQUAL 1.14)
      set(ADIOS2_HAVE_HDF5_VOL TRUE)
    else()
      set(ADIOS2_HAVE_HDF5_VOL FALSE)
      if(ADIOS2_USE_HDF5_VOL STREQUAL ON)
        message(FATAL_ERROR "ADIOS HDF5_VOL requires HDF5>=1.14+")
      endif()
    endif()
  endif()
endif()

# IME
if(ADIOS2_USE_IME STREQUAL AUTO)
  find_package(IME)
elseif(ADIOS2_USE_IME)
  find_package(IME REQUIRED)
endif()
if(IME_FOUND)
  set(ADIOS2_HAVE_IME TRUE)
endif()


# Python

# Not supported on PGI
if(CMAKE_CXX_COMPILER_ID STREQUAL PGI)
  if(ADIOS2_USE_Python STREQUAL ON)
    message(FATAL_ERROR "Python bindings are not supported with the PGI compiler")
  elseif(ADIOS2_USE_Python STREQUAL AUTO)
    message(WARNING "Disabling python bindings as they are not supported with the PGI compiler")
    set(ADIOS2_USE_Python OFF)
  endif()
endif()

# Not supported without shared libs
if(NOT SHARED_LIBS_SUPPORTED)
  if(ADIOS2_USE_Python STREQUAL ON)
    message(FATAL_ERROR "Python bindings are not supported without shared library support")
  elseif(ADIOS2_USE_Python STREQUAL AUTO)
    message(WARNING "Disabling python bindings since no shared library support was detected.")
    set(ADIOS2_USE_Python OFF)
  endif()
endif()

if(ADIOS2_USE_PIP)
  find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development.Module)
  set(ADIOS2_HAVE_PIP TRUE)
elseif(ADIOS2_USE_Python STREQUAL AUTO)
  find_package(Python 3.8 COMPONENTS Interpreter Development)
  if(Python_FOUND AND ADIOS2_HAVE_MPI)
    find_package(PythonModule COMPONENTS mpi4py mpi4py/mpi4py.h)
  endif()
elseif(ADIOS2_USE_Python)
  find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development)
  if(ADIOS2_HAVE_MPI)
    find_package(PythonModule REQUIRED COMPONENTS mpi4py mpi4py/mpi4py.h)
  endif()
endif()
if(Python_FOUND)
  if(ADIOS2_HAVE_MPI)
    if(PythonModule_mpi4py_FOUND)
      set(ADIOS2_HAVE_Python ON)
    endif()
  else()
    set(ADIOS2_HAVE_Python ON)
  endif()
endif()

# Even if no python support, we still want the interpreter for tests
if(BUILD_TESTING AND NOT Python_Interpreter_FOUND)
  find_package(Python REQUIRED COMPONENTS Interpreter)
endif()

if(Python_Interpreter_FOUND)
  # Setup output directories
  if(Python_Development_FOUND)
    lists_get_prefix("Python_INCLUDE_DIRS;Python_LIBRARIES;Python_SITEARCH" _Python_DEVPREFIX)
  else()
    lists_get_prefix("Python_EXECUTABLE;Python_SITEARCH" _Python_DEVPREFIX)
  endif()
  string_strip_prefix(
    "${_Python_DEVPREFIX}" "${Python_SITEARCH}" CMAKE_INSTALL_PYTHONDIR_DEFAULT
  )
  set(CMAKE_INSTALL_PYTHONDIR "${CMAKE_INSTALL_PYTHONDIR_DEFAULT}"
    CACHE PATH "Install directory for python modules"
  )
  mark_as_advanced(CMAKE_INSTALL_PYTHONDIR)
  set(CMAKE_PYTHON_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_PYTHONDIR}
  )
endif()

# Sst
if(ADIOS2_USE_SST AND NOT WIN32)
  set(ADIOS2_HAVE_SST TRUE)
  find_package(LIBFABRIC 1.6)
  if(LIBFABRIC_FOUND)
    set(ADIOS2_SST_HAVE_LIBFABRIC TRUE)
    find_package(CrayDRC)
    if(CrayDRC_FOUND)
      set(ADIOS2_SST_HAVE_CRAY_DRC TRUE)
    endif()

    try_compile(ADIOS2_SST_HAVE_CRAY_CXI
      ${ADIOS2_BINARY_DIR}/check_libfabric_cxi
      ${ADIOS2_SOURCE_DIR}/cmake/check_libfabric_cxi.c
      CMAKE_FLAGS
        "-DINCLUDE_DIRECTORIES=${LIBFABRIC_INCLUDE_DIRS}"
        "-DLINK_DIRECTORIES=${LIBFABRIC_LIBRARIES}")
    message(STATUS "Libfabric support for the HPE CXI provider: ${ADIOS2_SST_HAVE_CRAY_CXI}")
  endif()
  if(ADIOS2_HAVE_MPI AND NOT "${ADIOS2_SST_HAVE_MPI_DP}")
    set(CMAKE_REQUIRED_LIBRARIES "MPI::MPI_C;Threads::Threads")
    include(CheckCXXSourceRuns)
    check_cxx_source_runs([=[
        #include <chrono>
        #include <future>
        #include <mpi.h>
	#include <thread>
        #include <stdlib.h>

        #if !defined(MPICH)
        #error "MPICH is the only supported library"
        #endif

        int main()
        {
          // Timeout after 5 second
          auto task = std::async(std::launch::async, []() {
                                 std::this_thread::sleep_for(std::chrono::seconds(5));
                                 exit(EXIT_FAILURE);
          });

          char* port_name = new char[MPI_MAX_PORT_NAME];
          MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, NULL);
          MPI_Open_port(MPI_INFO_NULL, port_name);
          MPI_Close_port(port_name);
          MPI_Finalize();
          exit(EXIT_SUCCESS);
        }]=]
    ADIOS2_HAVE_MPI_CLIENT_SERVER)
    unset(CMAKE_REQUIRED_LIBRARIES)
    if (ADIOS2_HAVE_MPI_CLIENT_SERVER)
      set(ADIOS2_SST_HAVE_MPI_DP TRUE)
    else()
      if ("${ADIOS2_SST_EXPECT_MPI_DP}") 
          message(FATAL_ERROR "Expected MPI to support Client-server connection model, but test failed.")
      endif()
    endif()
  endif()
  # UCX
  if(ADIOS2_USE_UCX STREQUAL AUTO)
    find_package(UCX 1.9.0)
  elseif(ADIOS2_USE_UCX)
    find_package(UCX 1.9.0)
  endif()
  if(UCX_FOUND)
    set(ADIOS2_SST_HAVE_UCX TRUE)
    set(ADIOS2_HAVE_UCX TRUE)
  endif()
endif()

# DAOS
find_package(DAOS)
if(DAOS_FOUND)
  set(ADIOS2_HAVE_DAOS TRUE)
endif()

#SysV IPC
if(UNIX)
  include(CheckSymbolExists)
  CHECK_SYMBOL_EXISTS(shmget "sys/ipc.h;sys/shm.h" HAVE_shmget)
  if(HAVE_shmget)
    set(ADIOS2_HAVE_SysVShMem ON)
  else()
    set(ADIOS2_HAVE_SysVShMem OFF)
  endif()
else()
  set(ADIOS2_HAVE_SysVShMem OFF)
endif()

#Profiling
if(ADIOS2_USE_Profiling STREQUAL AUTO)
  if(BUILD_SHARED_LIBS)
    set(ADIOS2_HAVE_Profiling ON)
  else()
    set(ADIOS2_HAVE_Profiling OFF)
  endif()
elseif(ADIOS2_USE_Profiling)
  set(ADIOS2_HAVE_Profiling ON)
endif()

if(ADIOS2_USE_Endian_Reverse STREQUAL ON)
  set(ADIOS2_HAVE_Endian_Reverse TRUE)
endif()

# Sodium for EncryptionOperator
if(ADIOS2_USE_Sodium STREQUAL AUTO)
  find_package(Sodium)
elseif(ADIOS2_USE_Sodium)
  find_package(Sodium REQUIRED)
endif()
if(Sodium_FOUND)
  set(ADIOS2_HAVE_Sodium TRUE)
endif()

# Catalyst stub library for ParaViewFidesEngine plugin for in situ vis
if(ADIOS2_USE_Catalyst STREQUAL AUTO)
  find_package(catalyst 2.0 QUIET)
elseif(ADIOS2_USE_Catalyst)
  find_package(catalyst 2.0 REQUIRED)
endif()
if(catalyst_FOUND)
  set(ADIOS2_HAVE_Catalyst TRUE)
endif()

# AWS S3
if(ADIOS2_USE_AWSSDK STREQUAL AUTO)
    find_package(AWSSDK 1.10.15 COMPONENTS s3)
elseif(ADIOS2_USE_AWSSDK)
    find_package(AWSSDK 1.10.15 REQUIRED COMPONENTS s3)
endif()
if(AWSSDK_FOUND)
    set(ADIOS2_HAVE_AWSSDK TRUE)
endif()

# Campaign Management
if(ADIOS2_USE_Campaign STREQUAL AUTO)
  find_package(SQLite3)
  if (SQLite3_FOUND)
    find_package(ZLIB)
  endif()
  if (SQLite3_FOUND AND ZLIB_FOUND)
    message(STATUS "Sqlite3 and ZLIB found. Turn on Campaign Management")
    set(ADIOS2_HAVE_Campaign TRUE)
  endif()
elseif(ADIOS2_USE_Campaign)
  find_package(SQLite3 REQUIRED)
  if (SQLite3_FOUND)
    find_package(ZLIB REQUIRED)
  endif()
  if (SQLite3_FOUND AND ZLIB_FOUND)
    message(STATUS "Sqlite3 and ZLIB found. Turn on Campaign Management")
    set(ADIOS2_HAVE_Campaign TRUE)
  endif()
endif()  

# Multithreading
find_package(Threads REQUIRED)

# Floating point detection
include(CheckTypeRepresentation)

#check_float_type_representation(float FLOAT_TYPE_C)
#check_float_type_representation(double DOUBLE_TYPE_C)
#check_float_type_representation("long double" LONG_DOUBLE_TYPE_C)

if(ADIOS2_HAVE_Fortran)
  #check_float_type_representation(real REAL_TYPE_Fortran LANGUAGE Fortran)
  #check_float_type_representation("real(kind=4)" REAL4_TYPE_Fortran LANGUAGE Fortran)
  #check_float_type_representation("real(kind=8)" REAL8_TYPE_Fortran LANGUAGE Fortran)
  #check_float_type_representation("real(kind=16)" REAL16_TYPE_Fortran LANGUAGE Fortran)

  include(CheckFortranCompilerFlag)
  check_fortran_compiler_flag("-fallow-argument-mismatch" ADIOS2_USE_Fortran_flag_argument_mismatch)
endif()

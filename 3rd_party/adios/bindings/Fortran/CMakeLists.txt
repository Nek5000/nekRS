#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
include(FortranCInterface)
FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")
FortranCInterface_VERIFY(CXX QUIET)

# Check whether the compiler supports Fortran submodule constructs we need.
adios2_check_fortran_submodules(ADIOS2_HAVE_FORTRAN_SUBMODULES)

# Leaving this workaround in place but commented out just in case we need to
# re-enable it in the future
#
# Cray submodules have naming issues so just don't use them for now
#if(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
#  set(ADIOS2_HAVE_FORTRAN_SUBMODULES 0 CACHE INTERNAL "" FORCE)
#endif()

if(ADIOS2_USE_Fortran_flag_argument_mismatch)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
endif()

add_library(adios2_fortran
  f2c/adios2_f2c_adios.cpp
  f2c/adios2_f2c_attribute.cpp
  f2c/adios2_f2c_io.cpp
  f2c/adios2_f2c_variable.cpp
  f2c/adios2_f2c_engine.cpp
  f2c/adios2_f2c_operator.cpp

  modules/adios2_mod.f90
  modules/adios2_functions_mod.f90
  modules/adios2_functions_allocate_mod.f90
  modules/adios2_parameters_mod.f90
  modules/adios2_adios_mod.f90
  modules/adios2_adios_init_mod.F90
  modules/adios2_adios_init_serial_smod.F90
  modules/adios2_attribute_mod.f90
  modules/adios2_attribute_data_mod.f90
  modules/adios2_io_mod.f90
  modules/adios2_io_open_mod.F90
  modules/adios2_io_open_serial_smod.F90
  modules/adios2_io_define_variable_mod.f90
  modules/adios2_io_define_derived_variable_mod.f90
  modules/adios2_io_define_attribute_mod.f90
  modules/adios2_engine_mod.f90
  modules/adios2_engine_begin_step_mod.f90
  modules/adios2_engine_put_mod.f90
  modules/adios2_engine_get_mod.f90
  modules/adios2_variable_mod.f90
  modules/adios2_variable_min_mod.f90
  modules/adios2_variable_max_mod.f90
  modules/adios2_operator_mod.f90
)

set_property(TARGET adios2_fortran PROPERTY EXPORT_NAME fortran)
set_property(TARGET adios2_fortran PROPERTY OUTPUT_NAME adios2${ADIOS2_LIBRARY_SUFFIX}_fortran)

target_include_directories(adios2_fortran
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/adios2/fortran>
  PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(adios2_fortran
  PRIVATE
    adios2_c
    adios2_core
)

if(ADIOS2_HAVE_FORTRAN_SUBMODULES)
  target_compile_definitions(adios2_fortran PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:ADIOS2_HAVE_FORTRAN_SUBMODULES>")
endif()

# Set library version information
set_target_properties(adios2_fortran PROPERTIES
  VERSION ${ADIOS2_LIBRARY_VERSION}
  SOVERSION ${ADIOS2_LIBRARY_SOVERSION}
)

add_library(adios2::fortran ALIAS adios2_fortran)

if(ADIOS2_HAVE_MPI)
  target_compile_definitions(adios2_fortran PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:ADIOS2_HAVE_MPI_F>")

  set(adios2_fortran_mpi_srcs
    modules/adios2_adios_init_mpi_smod.F90
    modules/adios2_io_open_mpi_smod.F90
    f2c/adios2_f2c_adios_mpi.cpp
    f2c/adios2_f2c_io_mpi.cpp
    )

  if(ADIOS2_HAVE_FORTRAN_SUBMODULES)
    # With submodules we can put MPI implementations in a separate library.
    add_library(adios2_fortran_mpi ${adios2_fortran_mpi_srcs})
    target_include_directories(adios2_fortran_mpi PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
    target_link_libraries(adios2_fortran_mpi
      PUBLIC
        adios2_fortran
      PRIVATE
        adios2_c_mpi
        adios2_core_mpi
      INTERFACE
        MPI::MPI_Fortran
    )
    target_compile_definitions(adios2_fortran_mpi PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:ADIOS2_HAVE_FORTRAN_SUBMODULES;ADIOS2_HAVE_MPI_F>")
    set_target_properties(adios2_fortran_mpi PROPERTIES
      VERSION ${ADIOS2_LIBRARY_VERSION}
      SOVERSION ${ADIOS2_LIBRARY_SOVERSION}
    )
    set_property(TARGET adios2_fortran_mpi PROPERTY OUTPUT_NAME adios2${ADIOS2_LIBRARY_SUFFIX}_fortran_mpi)
  else()
    # Without submodules we need to put everything in one library.
    add_library(adios2_fortran_mpi INTERFACE)
    target_link_libraries(adios2_fortran_mpi INTERFACE adios2_fortran MPI::MPI_Fortran)
    target_sources(adios2_fortran PRIVATE ${adios2_fortran_mpi_srcs})
    target_link_libraries(adios2_fortran PRIVATE adios2_c_mpi adios2_core_mpi)
  endif()

  set_property(TARGET adios2_fortran_mpi PROPERTY EXPORT_NAME fortran_mpi)
  target_compile_definitions(adios2_fortran_mpi INTERFACE ADIOS2_USE_MPI)
  add_library(adios2::fortran_mpi ALIAS adios2_fortran_mpi)
  set(maybe_adios2_fortran_mpi adios2_fortran_mpi)
else()
  set(maybe_adios2_fortran_mpi)
endif()

add_library(adios2_f INTERFACE)
target_link_libraries(adios2_f INTERFACE ${maybe_adios2_fortran_mpi} adios2_fortran)

install(
  TARGETS adios2_fortran ${maybe_adios2_fortran_mpi} adios2_f EXPORT adios2FortranExports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT adios2_fortran-runtime
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT adios2_fortran-libraries NAMELINK_COMPONENT adios2_fortran-development
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT adios2_fortran-development
)
install(
  DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/adios2/fortran
  COMPONENT adios2_fortran-development
  FILES_MATCHING
    PATTERN "adios2*.mod"
    PATTERN "adios2*.smod"
    PATTERN "ADIOS2*.mod"
    PATTERN "ADIOS2*.smod"
    PATTERN "CMakeFiles" EXCLUDE
)

export(EXPORT adios2FortranExports
  FILE ${ADIOS2_BINARY_DIR}/adios2-fortran-targets.cmake
  NAMESPACE adios2::
)
install(EXPORT adios2FortranExports
  FILE adios2-fortran-targets.cmake
  NAMESPACE adios2::
  DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
  COMPONENT adios2_fortran-development
)

cmake_minimum_required(VERSION 3.12)

project(adios2-config-dummy C CXX)

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

find_package(adios2 REQUIRED)

if(ADIOS2_HAVE_Fortran)
  enable_language(Fortran)
endif()

add_executable(serial_without_C foo.c main.c)
add_executable(serial_with_C foo.c main.c)
target_compile_definitions(serial_with_C PRIVATE WITH_ADIOS2)
target_link_libraries(serial_with_C adios2::c)

add_executable(serial_without_CXX foo.cxx main.cxx)
add_executable(serial_with_CXX foo.cxx main.cxx)
target_compile_definitions(serial_with_CXX PRIVATE WITH_ADIOS2)
target_link_libraries(serial_with_CXX adios2::cxx11)

if(ADIOS2_HAVE_Fortran)
  add_executable(serial_without_Fortran foo.F90 main.f90)
  add_executable(serial_with_Fortran foo.F90 main.f90)
  target_compile_definitions(serial_with_Fortran PRIVATE WITH_ADIOS2)
  target_link_libraries(serial_with_Fortran adios2::fortran)
  set_target_properties(serial_without_Fortran serial_with_Fortran PROPERTIES
    LINKER_LANGUAGE Fortran
    )
endif()

if(ADIOS2_HAVE_MPI)
  find_package(MPI REQUIRED)

  add_executable(mpi_without_C foo.c main.c)
  target_link_libraries(mpi_without_C MPI::MPI_C)
  add_executable(mpi_with_C foo.c main.c)
  target_compile_definitions(mpi_with_C PRIVATE WITH_ADIOS2)
  target_link_libraries(mpi_with_C adios2::c_mpi adios2::c MPI::MPI_C)

  add_executable(mpi_without_CXX foo.cxx main.cxx)
  target_link_libraries(mpi_without_CXX MPI::MPI_C)
  add_executable(mpi_with_CXX foo.cxx main.cxx)
  target_compile_definitions(mpi_with_CXX PRIVATE WITH_ADIOS2)
  target_link_libraries(mpi_with_CXX adios2::cxx11_mpi adios2::cxx11 MPI::MPI_CXX)

  if(ADIOS2_HAVE_Fortran)
    add_executable(mpi_without_Fortran foo.F90 main.f90)
    target_link_libraries(mpi_without_Fortran MPI::MPI_Fortran)
    add_executable(mpi_with_Fortran foo.F90 main.f90)
    target_compile_definitions(mpi_with_Fortran PRIVATE WITH_ADIOS2)
    target_link_libraries(mpi_with_Fortran adios2::fortran_mpi adios2::fortran MPI::MPI_Fortran)
    set_target_properties(mpi_without_Fortran mpi_with_Fortran PROPERTIES
      LINKER_LANGUAGE Fortran
      )
  endif()
endif()

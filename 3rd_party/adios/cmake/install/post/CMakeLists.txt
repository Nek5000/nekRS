if(ADIOS2_HAVE_Fortran)
  set(ADIOS2_CONFIG_FORTRAN 1)
else()
  set(ADIOS2_CONFIG_FORTRAN 0)
endif()

if(ADIOS2_HAVE_MPI)
  set(ADIOS2_CONFIG_MPI 1)
else()
  set(ADIOS2_CONFIG_MPI 0)
endif()

if(ADIOS2_HAVE_MPI AND MPI_C_COMPILER)
  set(ADIOS2_CONFIG_MPICC 1)
else()
  set(ADIOS2_CONFIG_MPICC 0)
endif()

if(ADIOS2_HAVE_MPI AND MPI_CXX_COMPILER)
  set(ADIOS2_CONFIG_MPICXX 1)
else()
  set(ADIOS2_CONFIG_MPICXX 0)
endif()

if(ADIOS2_HAVE_MPI AND MPI_Fortran_COMPILER)
  set(ADIOS2_CONFIG_MPIFC 1)
else()
  set(ADIOS2_CONFIG_MPIFC 0)
endif()

if(NOT WIN32)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/generate-adios2-config.sh.in
    ${CMAKE_CURRENT_BINARY_DIR}/generate-adios2-config.sh
    @ONLY
  )
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/adios2-config.pre.sh.in
    ${CMAKE_CURRENT_BINARY_DIR}/adios2-config.pre.sh
    @ONLY
  )
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/adios2-config.post.sh.in
    ${CMAKE_CURRENT_BINARY_DIR}/adios2-config.post.sh
    @ONLY
  )
  install(CODE "
    execute_process(COMMAND ${CMAKE_CURRENT_BINARY_DIR}/generate-adios2-config.sh \"\${CMAKE_INSTALL_PREFIX}\" RESULT_VARIABLE result)
    if(result)
      message(FATAL_ERROR \"generate-adios2-config.sh returned \${result}, should be 0\")
    endif()
  "
    COMPONENT adios2_core-config
  )
endif()

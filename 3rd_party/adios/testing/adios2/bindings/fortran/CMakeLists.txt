#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

macro(fortran_add_test_helper testname mpi)
  set(test_targets "")
  if(NOT "${mpi}" MATCHES "^MPI_(ALLOW|ONLY|NONE)$")
    message(FATAL_ERROR "Invalid mpi argument value '${mpi}'.")
  endif()
  set(pfx Bindings.Fortran.)
  if(NOT "${mpi}" STREQUAL "MPI_ONLY")
    set(tgt Test.${pfx}${testname}.Serial)
    list(APPEND test_targets "${tgt}")
    add_executable(${tgt} Test${testname}.F90)
    set_target_properties(${tgt} PROPERTIES LINKER_LANGUAGE Fortran)
    target_link_libraries(${tgt} adios2::fortran)
    if (ADIOS2_HAVE_Derived_Variable)
        target_compile_definitions(${tgt} PRIVATE -DADIOS2_HAVE_Derived_Variable=1)
    endif()
    add_test(
      NAME ${pfx}${testname}.Serial
      COMMAND ${tgt}
    )
  endif()
  if(ADIOS2_HAVE_MPI AND NOT "${mpi}" STREQUAL "MPI_NONE")
    set(tgt Test.${pfx}${testname}.MPI)
    list(APPEND test_targets "${tgt}")
    add_executable(${tgt} Test${testname}.F90)
    set_target_properties(${tgt} PROPERTIES LINKER_LANGUAGE Fortran)
    target_link_libraries(${tgt} adios2::fortran_mpi MPI::MPI_Fortran)
    if (ADIOS2_HAVE_Derived_Variable)
        target_compile_definitions(${tgt} PRIVATE -DADIOS2_HAVE_Derived_Variable=1)
    endif()
    add_test(
      NAME ${pfx}${testname}.MPI
      COMMAND ${MPIEXEC_COMMAND} $<TARGET_FILE:${tgt}>
    )
    set_tests_properties(${pfx}${testname}.MPI PROPERTIES
      PROCESSORS "${MPIEXEC_MAX_NUMPROCS}"
    )
  endif()
  set("Test.${pfx}${testname}-TARGETS" "${test_targets}")
endmacro()

fortran_add_test_helper(BPWriteTypes MPI_ALLOW)
fortran_add_test_helper(Remove MPI_ALLOW)

fortran_add_test_helper(Adios2BindingsFortranIO MPI_ONLY)
fortran_add_test_helper(BPWriteReadAttributes MPI_ONLY)
fortran_add_test_helper(BPWriteVariableAttributes MPI_ONLY)
fortran_add_test_helper(BPWriteTypesByName MPI_ONLY)
fortran_add_test_helper(BPWriteTypesLocal MPI_ONLY)
fortran_add_test_helper(BPWriteReadHeatMap2D MPI_ONLY)
fortran_add_test_helper(BPWriteReadHeatMap3D MPI_ONLY)
fortran_add_test_helper(BPWriteReadHeatMap4D MPI_ONLY)
fortran_add_test_helper(BPWriteReadHeatMap5D MPI_ONLY)
fortran_add_test_helper(BPWriteReadHeatMap6D MPI_ONLY)
fortran_add_test_helper(BPReadGlobalsByName MPI_ONLY)
fortran_add_test_helper(BPWriteMemorySelectionRead2D MPI_ONLY)
fortran_add_test_helper(BPWriteMemorySelectionRead3D MPI_ONLY)
fortran_add_test_helper(NullEngine MPI_ONLY)
fortran_add_test_helper(BPMemorySpace MPI_NONE)
if(ADIOS2_HAVE_GPU_Support)
  fortran_add_test_helper(BPMemorySpaceGPU MPI_NONE)
endif()

if(ADIOS2_HAVE_MPI)
  add_subdirectory(operation) 
  
  # F2C 
  add_executable(Test.Bindings.Fortran.F2C.BPReadFBlocks
    TestF2C_BPReadFBlocks.cpp
  )
  target_link_libraries(Test.Bindings.Fortran.F2C.BPReadFBlocks
    adios2 adios2::thirdparty::gtest MPI::MPI_C
  )
  add_test(NAME Bindings.Fortran.F2C.BPReadFBlocks
    COMMAND
      ${MPIEXEC_COMMAND} $<TARGET_FILE:Test.Bindings.Fortran.F2C.BPReadFBlocks>
  )
  set_tests_properties(Bindings.Fortran.F2C.BPReadFBlocks PROPERTIES
    DEPENDS "Bindings.Fortran.BPWriteReadHeatMap2D;Bindings.Fortran.BPWriteReadHeatMap3D"
    PROCESSORS "${MPIEXEC_MAX_NUMPROCS}"
  )
endif()

add_library(SmallTestData_f OBJECT SmallTestData_mod.F90)
add_library(ISmallTestData_f INTERFACE)
target_sources(ISmallTestData_f INTERFACE $<TARGET_OBJECTS:SmallTestData_f>)

foreach(tgt
    ${Test.Bindings.Fortran.BPWriteTypes-TARGETS}
    ${Test.Bindings.Fortran.Remove-TARGETS}
    ${Test.Bindings.Fortran.BPWriteReadAttributes-TARGETS}
    ${Test.Bindings.Fortran.BPWriteVariableAttributes-TARGETS}
    ${Test.Bindings.Fortran.BPWriteTypesByName-TARGETS}
    ${Test.Bindings.Fortran.BPWriteTypesLocal-TARGETS}
    )
  target_link_libraries(${tgt} ISmallTestData_f)
endforeach()

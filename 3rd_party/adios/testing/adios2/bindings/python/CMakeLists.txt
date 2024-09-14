#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

function(add_python_mpi_test testname)
  python_add_test(NAME Bindings.Python.${testname}.MPI SCRIPT Test${testname}.py
    EXEC_WRAPPER ${MPIEXEC_COMMAND}
  )
  set_tests_properties(Bindings.Python.${testname}.MPI PROPERTIES
    PROCESSORS "${MPIEXEC_MAX_NUMPROCS}"
  )
endfunction()

python_add_test(NAME Bindings.Python.BPWriteReadTypes.Serial SCRIPT TestBPWriteReadTypes_nompi.py)
python_add_test(NAME Bindings.Python.TypeExceptionOnGet.Serial SCRIPT TestGetException_nompi.py)

if(ADIOS2_HAVE_MPI)
  add_python_mpi_test(BPWriteReadTypes)
  add_python_mpi_test(BPWriteReadString)
  add_python_mpi_test(BPReadMultisteps)
  add_python_mpi_test(BPWriteRead2D)
  add_python_mpi_test(BPBlocksInfo)
  add_python_mpi_test(NullEngine)
  add_python_mpi_test(Query)
  add_python_mpi_test(QueryLocalArray)
  # Currently hangs in H5Fclose for unknown reasons
  #if(ADIOS2_HAVE_HDF5)
  #  add_python_mpi_test(BPWriteTypesHighLevelAPI_HDF5)
  #endif()
endif()


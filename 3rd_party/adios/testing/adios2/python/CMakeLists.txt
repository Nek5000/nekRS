#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

function(add_python_mpi_test testname)
  python_add_test(NAME Api.Python.${testname}.MPI SCRIPT Test${testname}.py
    EXEC_WRAPPER ${MPIEXEC_COMMAND}
  )
  set_tests_properties(Api.Python.${testname}.MPI PROPERTIES
    PROCESSORS "${MPIEXEC_MAX_NUMPROCS}"
  )
endfunction()

python_add_test(NAME Api.Python.ADIOS SCRIPT TestADIOS.py)
python_add_test(NAME Api.Python.Engine SCRIPT TestEngine.py)
python_add_test(NAME Api.Python.IO SCRIPT TestIO.py)
python_add_test(NAME Api.Python.Operator SCRIPT TestOperator.py)
python_add_test(NAME Api.Python.Variable SCRIPT TestVariable.py)
python_add_test(NAME Api.Python.Attribute SCRIPT TestAttribute.py)
python_add_test(NAME Api.Python.Stream SCRIPT TestStream.py)
python_add_test(NAME Api.Python.FileReader SCRIPT TestFileReader.py)

if(ADIOS2_HAVE_MPI)
  add_python_mpi_test(BPWriteReadString)
  add_python_mpi_test(BPWriteTypesHighLevelAPI)
  add_python_mpi_test(BPWriteTypesHighLevelAPILocal)
  add_python_mpi_test(BPChangingShapeHighLevelAPI)
  if(ADIOS2_HAVE_PNG)
    add_python_mpi_test(BPPNGHighLevelAPI)
  endif()
  if(ADIOS2_HAVE_ZFP)
    add_python_mpi_test(BPZfpHighLevelAPI)
  endif()
  if(ADIOS2_HAVE_HDF5)
    python_add_test(NAME Api.Python.BPWriteTypesHighLevelAPI_HDF5 SCRIPT TestBPWriteTypesHighLevelAPI_HDF5.py)
    # MPI version currently fails in H5Fclose.
    # add_python_mpi_test(BPWriteTypesHighLevelAPI_HDF5)
  endif()
endif()

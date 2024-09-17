#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

include(ADIOSFunctions)

add_test(NAME HeatTransfer.BP3.zfp.Mx1.Write
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_EXTRA_FLAGS}
    ${MPIEXEC_NUMPROC_FLAG} 4
      $<TARGET_FILE:adios2_simulations_heatTransferWrite>
        ${PROJECT_SOURCE_DIR}/examples/simulations/heatTransfer/heat_bp3_zfp.xml
        Write.bp 2 2 10 10 10 10
)
set_tests_properties(HeatTransfer.BP3.zfp.Mx1.Write PROPERTIES PROCESSORS 4)

add_test(NAME HeatTransfer.BP3.zfp.Mx1.Read
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_EXTRA_FLAGS}
    ${MPIEXEC_NUMPROC_FLAG} 1
      $<TARGET_FILE:adios2_simulations_heatTransferRead>
        ${PROJECT_SOURCE_DIR}/examples/simulations/heatTransfer/heat_bp3.xml
        Write.bp Read.bp 1 1
)

add_test(NAME HeatTransfer.BP3.zfp.Mx1.Dump
  COMMAND ${CMAKE_COMMAND}
    -DARG1=-d 
    -DINPUT_FILE=Read.bp
    -DOUTPUT_FILE=Dump.txt
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

SetupTestPipeline(HeatTransfer.BP3.zfp.Mx1 "Write;Read;Dump" TRUE)

#############################################################################
add_test(NAME HeatTransfer.BP4.zfp.Mx1.Write
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_EXTRA_FLAGS}
    ${MPIEXEC_NUMPROC_FLAG} 4
      $<TARGET_FILE:adios2_simulations_heatTransferWrite>
        ${PROJECT_SOURCE_DIR}/examples/simulations/heatTransfer/heat_bp4_zfp.xml
        Write.bp 2 2 10 10 10 10
)
set_tests_properties(HeatTransfer.BP3.zfp.Mx1.Write PROPERTIES PROCESSORS 4)

add_test(NAME HeatTransfer.BP4.zfp.Mx1.Read
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_EXTRA_FLAGS}
    ${MPIEXEC_NUMPROC_FLAG} 1
      $<TARGET_FILE:adios2_simulations_heatTransferRead>
        ${PROJECT_SOURCE_DIR}/examples/simulations/heatTransfer/heat_bp4.xml
        Write.bp Read.bp 1 1
)

add_test(NAME HeatTransfer.BP4.zfp.Mx1.Dump
  COMMAND ${CMAKE_COMMAND}
    -DARG1=-d 
    -DINPUT_FILE=Read.bp
    -DOUTPUT_FILE=Dump.txt
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

SetupTestPipeline(HeatTransfer.BP4.zfp.Mx1 "Write;Read;Dump" TRUE)

#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

include(ADIOSFunctions)

add_executable(Test.Utils.ChangingShape TestUtilsChangingShape.cpp)

if(ADIOS2_HAVE_MPI)
  target_link_libraries(Test.Utils.ChangingShape adios2::cxx11_mpi MPI::MPI_C)
  set(cmd_executor ${MPIEXEC_COMMAND} ${MPIEXEC_NUMPROC_FLAG} 1)
else()
  target_link_libraries(Test.Utils.ChangingShape adios2::cxx11)
  set(cmd_executor)
endif()

add_test(NAME Utils.ChangingShape
  COMMAND ${cmd_executor} $<TARGET_FILE:Test.Utils.ChangingShape>
)

########################################
# bpls -la to screen for testing
########################################
add_test(NAME Utils.ChangingShape.Screen
  COMMAND ${CMAKE_COMMAND}
    -DARG1=-la
    -DINPUT_FILE=TestUtilsChangingShape.bp
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

add_test(NAME Utils.ChangingShape.Dump
  COMMAND ${CMAKE_COMMAND}
    -DARG1=-la
    -DINPUT_FILE=TestUtilsChangingShape.bp
    -DOUTPUT_FILE=TestUtilsChangingShape.bplsla.result.txt
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

if(ADIOS2_HAVE_MPI)
  add_test(NAME Utils.ChangingShape.Validate
    COMMAND ${DIFF_COMMAND} -u -w
      ${CMAKE_CURRENT_SOURCE_DIR}/TestUtilsChangingShape.bplsla.expected.txt
      TestUtilsChangingShape.bplsla.result.txt
  )
  SetupTestPipeline(Utils.ChangingShape ";Screen;Dump;Validate" TRUE)
else()
  SetupTestPipeline(Utils.ChangingShape ";Screen;Dump" TRUE)
endif()


########################################
# bpls -ld  AlternatingStepsAndChangingShapeVar
########################################
add_test(NAME Utils.ChangingShape.AlternatingStepsAndChangingShapeVar.Dump
  COMMAND ${CMAKE_COMMAND}
    -DARG1=-ld
    -DARG2=AlternatingStepsAndChangingShapeVar
    -DINPUT_FILE=TestUtilsChangingShape.bp
    -DOUTPUT_FILE=TestUtilsChangingShape.bplsldAlternatingStepsAndChangingShapeVar.result.txt
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

if(ADIOS2_HAVE_MPI)
  add_test(NAME Utils.ChangingShape.AlternatingStepsAndChangingShapeVar.Validate
    COMMAND ${DIFF_COMMAND} -u -w
      ${CMAKE_CURRENT_SOURCE_DIR}/TestUtilsChangingShape.bplsldAlternatingStepsAndChangingShapeVar.expected.txt
      TestUtilsChangingShape.bplsldAlternatingStepsAndChangingShapeVar.result.txt
  )
  SetupTestPipeline(Utils.ChangingShape
    ";AlternatingStepsAndChangingShapeVar.Dump;AlternatingStepsAndChangingShapeVar.Validate"
    FALSE
  )
else()
  SetupTestPipeline(Utils.ChangingShape
    ";AlternatingStepsAndChangingShapeVar.Dump" FALSE
  )
endif()

########################################
# bpls -ld AlternatingStepsVar -s "1,0,0" -c "3,-1,-1" -n 8
########################################
add_test(NAME Utils.ChangingShape.AlternatingStepsVarSelection.Dump
  COMMAND ${CMAKE_COMMAND}
    -DARG1=-ld
    -DARG2=AlternatingStepsVar
    -DARG3=-s
    -DARG4=1,0,0
    -DARG5=-c
    -DARG6=3,-1,-1
    -DARG7=-n
    -DARG8=8
    -DINPUT_FILE=TestUtilsChangingShape.bp
    -DOUTPUT_FILE=TestUtilsChangingShape.bplsldAlternatingStepsVarSelection.result.txt
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

if(ADIOS2_HAVE_MPI)
  add_test(NAME Utils.ChangingShape.AlternatingStepsVarSelection.Validate
    COMMAND ${DIFF_COMMAND} -u -w
      ${CMAKE_CURRENT_SOURCE_DIR}/TestUtilsChangingShape.bplsldAlternatingStepsVarSelection.expected.txt
      TestUtilsChangingShape.bplsldAlternatingStepsVarSelection.result.txt
  )
  SetupTestPipeline(Utils.ChangingShape
    ";AlternatingStepsVarSelection.Dump;AlternatingStepsVarSelection.Validate"
    FALSE
  )
else()
  SetupTestPipeline(Utils.ChangingShape
    ";AlternatingStepsVarSelection.Dump" FALSE
  )
endif()

########################################
# bpls -ld ChangingShapeVar -s "5,0,0" -c "1,-1,-1" -n 12
########################################
add_test(NAME Utils.ChangingShape.ChangingShapeVarOneStep.Dump
  COMMAND ${CMAKE_COMMAND}
    -DARG1=-ld
    -DARG2=ChangingShapeVar
    -DARG3=-s
    -DARG4=5,0,0
    -DARG5=-c
    -DARG6=1,-1,-1
    -DARG7=-n
    -DARG8=12
    -DINPUT_FILE=TestUtilsChangingShape.bp
    -DOUTPUT_FILE=TestUtilsChangingShape.bplsldChangingShapeVarOneStep.result.txt
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)
if(ADIOS2_HAVE_MPI)
  add_test(NAME Utils.ChangingShape.ChangingShapeVarOneStep.Validate
    COMMAND ${DIFF_COMMAND} -u -w
      ${CMAKE_CURRENT_SOURCE_DIR}/TestUtilsChangingShape.bplsldChangingShapeVarOneStep.expected.txt
      TestUtilsChangingShape.bplsldChangingShapeVarOneStep.result.txt
  )
  SetupTestPipeline(Utils.ChangingShape
    ";ChangingShapeVarOneStep.Dump;ChangingShapeVarOneStep.Validate"
    FALSE
  )
else()
  SetupTestPipeline(Utils.ChangingShape ";ChangingShapeVarOneStep.Dump" FALSE)
endif()

########################################
# bpls -ld FixedShapeVar -s "5,0,0" -c "10,-1,-1"
########################################
add_test(NAME Utils.ChangingShape.FixedShapeVarTooManySteps.Dump
  COMMAND ${CMAKE_COMMAND}
    -DARG1=-ld
    -DARG2=FixedShapeVar
    -DARG3=-s
    -DARG4=5,0,0
    -DARG5=-c
    -DARG6=10,-1,-1
    -DINPUT_FILE=TestUtilsChangingShape.bp
    -DOUTPUT_FILE=TestUtilsChangingShape.bplsldFixedShapeVarTooManySteps.result.txt
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

if(ADIOS2_HAVE_MPI)
  add_test(NAME Utils.ChangingShape.FixedShapeVarTooManySteps.Validate
    COMMAND ${DIFF_COMMAND} -u -w
      ${CMAKE_CURRENT_SOURCE_DIR}/TestUtilsChangingShape.bplsldFixedShapeVarTooManySteps.expected.txt
      TestUtilsChangingShape.bplsldFixedShapeVarTooManySteps.result.txt
  )
  SetupTestPipeline(Utils.ChangingShape
    ";FixedShapeVarTooManySteps.Dump;FixedShapeVarTooManySteps.Validate"
    FALSE
  )
else()
  SetupTestPipeline(Utils.ChangingShape ";FixedShapeVarTooManySteps.Dump" FALSE)
endif()



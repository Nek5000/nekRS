#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

include(ADIOSFunctions)

add_executable(Test.Utils.CWriter TestUtilsCWriter.c)

if(ADIOS2_HAVE_MPI)
  target_link_libraries(Test.Utils.CWriter adios2::c_mpi MPI::MPI_C)
  set(cmd_executor ${MPIEXEC_COMMAND} ${MPIEXEC_NUMPROC_FLAG} 1)
else()
  target_link_libraries(Test.Utils.CWriter adios2::c)
  set(cmd_executor)
endif()

add_test(NAME Utils.CWriter
  COMMAND ${cmd_executor} $<TARGET_FILE:Test.Utils.CWriter>
)
# This test produces TestUtilsCWriter.bp

########################################
# bpls -h
########################################
add_test(NAME Utils.CWriter.Bpls.h.Dump
  COMMAND ${CMAKE_COMMAND}
    -DARG1=-h
    -DINPUT_FILE=TestUtilsCWriter.bp
    -DOUTPUT_FILE=TestUtilsCWriter.bplsh.result.txt
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

if(ADIOS2_HAVE_MPI)
  add_test(NAME Utils.CWriter.Bpls.h.Validate
    COMMAND ${DIFF_COMMAND} -u -w
      ${CMAKE_CURRENT_SOURCE_DIR}/TestUtilsCWriter.bplsh.expected.txt
      TestUtilsCWriter.bplsh.result.txt
  )
  SetupTestPipeline(Utils.CWriter ";Bpls.h.Dump;Bpls.h.Validate" TRUE)
else()
  SetupTestPipeline(Utils.CWriter ";Bpls.h.Dump" TRUE)
endif()


########################################
# bpls -la
########################################
add_test(NAME Utils.CWriter.Bpls.la.Screen
  COMMAND ${CMAKE_COMMAND}
    -DARG1=-la
    -DINPUT_FILE=TestUtilsCWriter.bp
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

# bpls -la
add_test(NAME Utils.CWriter.Bpls.la.Dump
  COMMAND ${CMAKE_COMMAND}
    -DARG1=-la
    -DINPUT_FILE=TestUtilsCWriter.bp
    -DOUTPUT_FILE=TestUtilsCWriter.bplsla.result.txt
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

if(ADIOS2_HAVE_MPI)
  add_test(NAME Utils.CWriter.Bpls.la.Validate
    COMMAND ${DIFF_COMMAND} -u -w
      ${CMAKE_CURRENT_SOURCE_DIR}/TestUtilsCWriter.bplsla.expected.txt
      TestUtilsCWriter.bplsla.result.txt
  )
  SetupTestPipeline(Utils.CWriter
    ";Bpls.la.Screen;Bpls.la.Dump;Bpls.la.Validate" FALSE
  )
else()
  SetupTestPipeline(Utils.CWriter ";Bpls.la.Screen;Bpls.la.Dump" FALSE)
endif()


########################################
# bpls -Al
########################################
add_test(NAME Utils.CWriter.Bpls.Al.Dump
  COMMAND ${CMAKE_COMMAND}
    -DARG1=-Al
    -DINPUT_FILE=TestUtilsCWriter.bp
    -DOUTPUT_FILE=TestUtilsCWriter.bplsAl.result.txt
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

if(ADIOS2_HAVE_MPI)
  add_test(NAME Utils.CWriter.Bpls.Al.Validate
    COMMAND ${DIFF_COMMAND} -u -w
      ${CMAKE_CURRENT_SOURCE_DIR}/TestUtilsCWriter.bplsAl.expected.txt
      TestUtilsCWriter.bplsAl.result.txt
  )
  SetupTestPipeline(Utils.CWriter ";Bpls.Al.Dump;Bpls.Al.Validate" FALSE)
else()
  SetupTestPipeline(Utils.CWriter ";Bpls.Al.Dump" FALSE)
endif()

########################################
# bpls -ldDav
########################################
add_test(NAME Utils.CWriter.Bpls.ldDav.Dump
  COMMAND ${CMAKE_COMMAND}
    -DARG1='-ldDav'
    -DINPUT_FILE=TestUtilsCWriter.bp
    -DOUTPUT_FILE=TestUtilsCWriter.bplsldDav.result.txt
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

if(ADIOS2_HAVE_MPI)
  add_test(NAME Utils.CWriter.Bpls.ldDav.Validate
    COMMAND ${DIFF_COMMAND} -u -w
      ${CMAKE_CURRENT_SOURCE_DIR}/TestUtilsCWriter.bplsldDav.expected.txt
      TestUtilsCWriter.bplsldDav.result.txt
  )
  SetupTestPipeline(Utils.CWriter ";Bpls.ldDav.Dump;Bpls.ldDav.Validate" FALSE)
else()
  SetupTestPipeline(Utils.CWriter ";Bpls.ldDav.Dump" FALSE)
endif()


########################################
# bpls -ldDavvv
########################################
add_test(NAME Utils.CWriter.Bpls.ldDavvv.Dump
  COMMAND ${CMAKE_COMMAND}
    -DARG1='-ldDavvv'
    -DINPUT_FILE=TestUtilsCWriter.bp
    -DOUTPUT_FILE=TestUtilsCWriter.bplsldDavvv.result.txt
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

if(ADIOS2_HAVE_MPI)
  add_test(NAME Utils.CWriter.Bpls.ldDavvv.Validate
    COMMAND ${DIFF_COMMAND} -u -w
      ${CMAKE_CURRENT_SOURCE_DIR}/TestUtilsCWriter.bplsldDavvv.expected.txt
      TestUtilsCWriter.bplsldDavvv.result.txt
  )
  SetupTestPipeline(Utils.CWriter
    ";Bpls.ldDavvv.Dump;Bpls.ldDavvv.Validate" FALSE
  )
else()
  SetupTestPipeline(Utils.CWriter ";Bpls.ldDavvv.Dump" FALSE)
endif()


########################################
# bpls -ld varI16 -n 10
########################################
add_test(NAME Utils.CWriter.Bpls.ldvarI16.Dump
  COMMAND ${CMAKE_COMMAND}
    -DARG1=-ld
    -DARG2=varI16
    -DARG3=-n
    -DARG4=10
    -DINPUT_FILE=TestUtilsCWriter.bp
    -DOUTPUT_FILE=TestUtilsCWriter.bplsldvarI16.result.txt
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

if(ADIOS2_HAVE_MPI)
  add_test(NAME Utils.CWriter.Bpls.ldvarI16.Validate
    COMMAND ${DIFF_COMMAND} -u -w
      ${CMAKE_CURRENT_SOURCE_DIR}/TestUtilsCWriter.bplsldvarI16.expected.txt
      TestUtilsCWriter.bplsldvarI16.result.txt
  )
  SetupTestPipeline(Utils.CWriter
    ";Bpls.ldvarI16.Dump;Bpls.ldvarI16.Validate" FALSE
  )
else()
  SetupTestPipeline(Utils.CWriter ";Bpls.ldvarI16.Dump" FALSE)
endif()



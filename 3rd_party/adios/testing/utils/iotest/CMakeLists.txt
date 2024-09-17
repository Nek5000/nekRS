#---------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#---------------------------------------------------------------------#

# SetupTestPipeline
include(ADIOSFunctions)

#------------------------------------------
#  Pipe2 BP Write
#------------------------------------------
add_test(NAME Utils.IOTest.Pipe2.BP.Write
  COMMAND ${MPIEXEC_COMMAND} ${MPIEXEC_NUMPROC_FLAG} 2
      $<TARGET_FILE:adios_iotest>
        -a 1 -c ${CMAKE_CURRENT_SOURCE_DIR}/pipe2.txt
        -x ${CMAKE_CURRENT_SOURCE_DIR}/pipe2-BP.xml -d 2 1 --strong-scaling
)
set_tests_properties(Utils.IOTest.Pipe2.BP.Write PROPERTIES PROCESSORS 2)

add_test(NAME Utils.IOTest.Pipe2.BP.Write.Dump
  COMMAND ${CMAKE_COMMAND}
    -DARG1=-laD
    -DINPUT_FILE=pipe2_write.bp
    -DOUTPUT_FILE=IOTest.Pipe2.BP.Write.bpls.txt
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

add_test(NAME Utils.IOTest.Pipe2.BP.Write.Validate
  COMMAND ${DIFF_COMMAND} -u -w
    ${CMAKE_CURRENT_SOURCE_DIR}/IOTest.Pipe2.BP.Write.bpls.txt
    IOTest.Pipe2.BP.Write.bpls.txt
)

SetupTestPipeline(
  Utils.IOTest.Pipe2.BP
  "Write;Write.Dump;Write.Validate"
  TRUE
)

#------------------------------------------
#  Pipe2 BP Read
#------------------------------------------
add_test(NAME Utils.IOTest.Pipe2.BP.Read
  COMMAND ${MPIEXEC_COMMAND} ${MPIEXEC_NUMPROC_FLAG} 2
    $<TARGET_FILE:adios_iotest>
      -a 2 -c ${CMAKE_CURRENT_SOURCE_DIR}/pipe2.txt
      -x ${CMAKE_CURRENT_SOURCE_DIR}/pipe2-BP.xml -d 2 1 --strong-scaling
)
set_tests_properties(Utils.IOTest.Pipe2.BP.Read PROPERTIES PROCESSORS 2)

add_test(NAME Utils.IOTest.Pipe2.BP.Read.Dump
  COMMAND ${CMAKE_COMMAND}
    -DARG1=-laD
    -DINPUT_FILE=pipe2_read.bp
    -DOUTPUT_FILE=IOTest.Pipe2.BP.Read.bpls.txt
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

add_test(NAME Utils.IOTest.Pipe2.BP.Read.Validate
  COMMAND ${DIFF_COMMAND} -u -w
    ${CMAKE_CURRENT_SOURCE_DIR}/IOTest.Pipe2.BP.Read.bpls.txt
    IOTest.Pipe2.BP.Read.bpls.txt
)

SetupTestPipeline(
  Utils.IOTest.Pipe2.BP
  "Write;Read;Read.Dump;Read.Validate"
  FALSE
)

if(ADIOS2_HAVE_HDF5 AND HDF5_IS_PARALLEL)
  #------------------------------------------
  #  Pipe2 HDF5 Write
  #------------------------------------------
  add_test(NAME Utils.IOTest.Pipe2.HDF5.Write
    COMMAND ${MPIEXEC_COMMAND} ${MPIEXEC_NUMPROC_FLAG} 2
      $<TARGET_FILE:adios_iotest>
        -a 1 -c ${CMAKE_CURRENT_SOURCE_DIR}/pipe2-HDF5.txt
        -x ${CMAKE_CURRENT_SOURCE_DIR}/pipe2-HDF5.xml -d 2 1 --strong-scaling
  )
  set_tests_properties(Utils.IOTest.Pipe2.HDF5.Write PROPERTIES PROCESSORS 2)

  add_test(NAME Utils.IOTest.Pipe2.HDF5.Write.Dump
    COMMAND ${CMAKE_COMMAND}
      -DARG1=-laD
      -DINPUT_FILE=pipe2_write.h5
      -DOUTPUT_FILE=IOTest.Pipe2.HDF5.Write.bpls.txt
      -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
  )

  add_test(NAME Utils.IOTest.Pipe2.HDF5.Write.Validate
    COMMAND ${DIFF_COMMAND} -u -w
      ${CMAKE_CURRENT_SOURCE_DIR}/IOTest.Pipe2.HDF5.Write.bpls.txt
      IOTest.Pipe2.HDF5.Write.bpls.txt
  )

  SetupTestPipeline(
    Utils.IOTest.Pipe2.HDF5
    "Write;Write.Dump;Write.Validate"
    TRUE
  )

  #------------------------------------------
  #  Pipe2 HDF5 Read
  #------------------------------------------
  add_test(NAME Utils.IOTest.Pipe2.HDF5.Read
    COMMAND ${MPIEXEC_COMMAND} ${MPIEXEC_NUMPROC_FLAG} 2
      $<TARGET_FILE:adios_iotest>
        -a 2 -c ${CMAKE_CURRENT_SOURCE_DIR}/pipe2-HDF5.txt
        -x ${CMAKE_CURRENT_SOURCE_DIR}/pipe2-HDF5.xml -d 2 1 --strong-scaling
  )
  set_tests_properties(Utils.IOTest.Pipe2.HDF5.Read PROPERTIES PROCESSORS 2)

  add_test(NAME Utils.IOTest.Pipe2.HDF5.Read.Dump
    COMMAND ${CMAKE_COMMAND}
      -DARG1=-laD
      -DINPUT_FILE=pipe2_read.h5
      -DOUTPUT_FILE=IOTest.Pipe2.HDF5.Read.bpls.txt
      -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
  )

  add_test(NAME Utils.IOTest.Pipe2.HDF5.Read.Validate
    COMMAND ${DIFF_COMMAND} -u -w
      ${CMAKE_CURRENT_SOURCE_DIR}/IOTest.Pipe2.HDF5.Read.bpls.txt
      IOTest.Pipe2.HDF5.Read.bpls.txt
  )

  SetupTestPipeline(
    Utils.IOTest.Pipe2.HDF5
    "Write;Read;Read.Dump;Read.Validate"
    FALSE
  )
endif()

if(ADIOS2_RUN_MPI_MPMD_TESTS)
  if(MPIEXEC_MAX_NUMPROCS GREATER 1) 
    #------------------------------------------
    #  Pipe2 SSC version
    #------------------------------------------
    if(NOT MSVC)
      add_test(NAME Utils.IOTest.Pipe2.SSC
        COMMAND ${MPIEXEC_COMMAND}
          ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:adios_iotest>
            -a 1 -c ${CMAKE_CURRENT_SOURCE_DIR}/pipe2.txt
            -x ${CMAKE_CURRENT_SOURCE_DIR}/pipe2-SSC.xml -d 1 1
            --strong-scaling
          :
          ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:adios_iotest>
            -a 2 -c ${CMAKE_CURRENT_SOURCE_DIR}/pipe2.txt
            -x ${CMAKE_CURRENT_SOURCE_DIR}/pipe2-SSC.xml -d 1 1
            --strong-scaling 
      )
      set_tests_properties(Utils.IOTest.Pipe2.SSC PROPERTIES PROCESSORS 2)

      add_test(NAME Utils.IOTest.Pipe2.SSC.Dump
        COMMAND ${CMAKE_COMMAND}
          -DARG1=-laD
          -DINPUT_FILE=pipe2_read.bp
          -DOUTPUT_FILE=IOTest.Pipe2.SSC.bpls.txt
          -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
      )

      add_test(NAME Utils.IOTest.Pipe2.SSC.Validate
        COMMAND ${DIFF_COMMAND} -u -w
          ${CMAKE_CURRENT_SOURCE_DIR}/IOTest.Pipe2.SSC.bpls.txt
          IOTest.Pipe2.SSC.bpls.txt
      )

      SetupTestPipeline(
        Utils.IOTest.Pipe2.SSC
        ";Dump;Validate"
        TRUE
      )
    endif()

    #------------------------------------------
    #  Pipe2 SST version 
    #------------------------------------------
    if(ADIOS2_HAVE_SST)
      add_test(NAME Utils.IOTest.Pipe2.SST
        COMMAND ${MPIEXEC_COMMAND}
          ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:adios_iotest>
            -a 1 -c ${CMAKE_CURRENT_SOURCE_DIR}/pipe2.txt
            -x ${CMAKE_CURRENT_SOURCE_DIR}/pipe2-SST.xml -d 1 1 --strong-scaling
          :
          ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:adios_iotest>
            -a 2 -c ${CMAKE_CURRENT_SOURCE_DIR}/pipe2.txt
            -x ${CMAKE_CURRENT_SOURCE_DIR}/pipe2-SST.xml -d 1 1 --strong-scaling 
      )
  	
      add_test(NAME Utils.IOTest.Pipe2.SST.Dump
        COMMAND ${CMAKE_COMMAND}
          -DARG1=-laD
          -DINPUT_FILE=pipe2_read.bp
          -DOUTPUT_FILE=IOTest.Pipe2.SST.bpls.txt
          -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
      )
  
      add_test(NAME Utils.IOTest.Pipe2.SST.Validate
        COMMAND ${DIFF_COMMAND} -u -w
          ${CMAKE_CURRENT_SOURCE_DIR}/IOTest.Pipe2.SST.bpls.txt
          IOTest.Pipe2.SST.bpls.txt
      )
  
      SetupTestPipeline(
        Utils.IOTest.Pipe2.SST
        ";Dump;Validate"
        TRUE
      )
    endif()
  else()
    message(WARNING "Not enough (2) MPI procs for insitu tests with adios_iotest")
  endif() 
endif()

#------------------------------------------
#  BurstBuffer BP4 tests
#------------------------------------------

#------------------------------------------
#  BurstBuffer BP4 test nproc = 1
#------------------------------------------
add_test(NAME Utils.IOTest.BurstBuffer.1to1.BP4.Write
  COMMAND ${MPIEXEC_COMMAND} ${MPIEXEC_NUMPROC_FLAG} 1
      $<TARGET_FILE:adios_iotest>
        -a 1 -c ${CMAKE_CURRENT_SOURCE_DIR}/burstbuffer-32KB.txt
        -x ${CMAKE_CURRENT_SOURCE_DIR}/burstbuffer-BP4-nsub1.xml -d 1 1 --strong-scaling
)
set_tests_properties(Utils.IOTest.BurstBuffer.1to1.BP4.Write PROPERTIES PROCESSORS 1)

add_test(NAME Utils.IOTest.BurstBuffer.1to1.BP4.Write.Dump
  COMMAND ${CMAKE_COMMAND}
    -DARG1=-laD
    -DINPUT_FILE=burstbuffer-32KB.bp
    -DOUTPUT_FILE=IOTest.BurstBuffer.1to1.BP4.Write.bpls.txt
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

add_test(NAME Utils.IOTest.BurstBuffer.1to1.BP4.Write.Validate
  COMMAND ${DIFF_COMMAND} -u -w
    ${CMAKE_CURRENT_SOURCE_DIR}/IOTest.BurstBuffer.nproc1.BP4.Write.bpls.txt
    IOTest.BurstBuffer.1to1.BP4.Write.bpls.txt
)

SetupTestPipeline(
  Utils.IOTest.BurstBuffer.1to1.BP4
  "Write;Write.Dump;Write.Validate"
  TRUE
)


#-------------------------------------------------
#  BurstBuffer BP4 test nproc = 2 aggregate 2 to 1
#-------------------------------------------------
add_test(NAME Utils.IOTest.BurstBuffer.2to1.BP4.Write
  COMMAND ${MPIEXEC_COMMAND} ${MPIEXEC_NUMPROC_FLAG} 2
      $<TARGET_FILE:adios_iotest>
        -a 1 -c ${CMAKE_CURRENT_SOURCE_DIR}/burstbuffer-32KB.txt
        -x ${CMAKE_CURRENT_SOURCE_DIR}/burstbuffer-BP4-nsub1.xml -d 2 1 --weak-scaling
)
set_tests_properties(Utils.IOTest.BurstBuffer.2to1.BP4.Write PROPERTIES PROCESSORS 2)

add_test(NAME Utils.IOTest.BurstBuffer.2to1.BP4.Write.Dump
  COMMAND ${CMAKE_COMMAND}
    -DARG1=-laD
    -DINPUT_FILE=burstbuffer-32KB.bp
    -DOUTPUT_FILE=IOTest.BurstBuffer.2to1.BP4.Write.bpls.txt
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

add_test(NAME Utils.IOTest.BurstBuffer.2to1.BP4.Write.Validate
  COMMAND ${DIFF_COMMAND} -u -w
    ${CMAKE_CURRENT_SOURCE_DIR}/IOTest.BurstBuffer.nproc2.BP4.Write.bpls.txt
    IOTest.BurstBuffer.2to1.BP4.Write.bpls.txt
)

SetupTestPipeline(
  Utils.IOTest.BurstBuffer.2to1.BP4
  "Write;Write.Dump;Write.Validate"
  TRUE
)


#-------------------------------------------------
#  BurstBuffer BP4 test nproc = 2 no aggregation
#-------------------------------------------------
add_test(NAME Utils.IOTest.BurstBuffer.2to2.BP4.Write
  COMMAND ${MPIEXEC_COMMAND} ${MPIEXEC_NUMPROC_FLAG} 2
      $<TARGET_FILE:adios_iotest>
        -a 1 -c ${CMAKE_CURRENT_SOURCE_DIR}/burstbuffer-32KB.txt
        -x ${CMAKE_CURRENT_SOURCE_DIR}/burstbuffer-BP4-nsub2.yaml -d 2 1 --weak-scaling
)
set_tests_properties(Utils.IOTest.BurstBuffer.2to2.BP4.Write PROPERTIES PROCESSORS 2)

add_test(NAME Utils.IOTest.BurstBuffer.2to2.BP4.Write.Dump
  COMMAND ${CMAKE_COMMAND}
    -DARG1=-laD
    -DINPUT_FILE=burstbuffer-32KB.bp
    -DOUTPUT_FILE=IOTest.BurstBuffer.2to2.BP4.Write.bpls.txt
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

add_test(NAME Utils.IOTest.BurstBuffer.2to2.BP4.Write.Validate
  COMMAND ${DIFF_COMMAND} -u -w
    ${CMAKE_CURRENT_SOURCE_DIR}/IOTest.BurstBuffer.nproc2.BP4.Write.bpls.txt
    IOTest.BurstBuffer.2to2.BP4.Write.bpls.txt
)

SetupTestPipeline(
  Utils.IOTest.BurstBuffer.2to2.BP4
  "Write;Write.Dump;Write.Validate"
  TRUE
)

#-------------------------------------------------------
#  BurstBuffer BP4 test nproc = 1, NO Draining to target
#-------------------------------------------------------
add_test(NAME Utils.IOTest.BurstBufferNoDrain.1to1.BP4.Write
  COMMAND ${MPIEXEC_COMMAND} ${MPIEXEC_NUMPROC_FLAG} 1
      $<TARGET_FILE:adios_iotest>
        -a 1 -c ${CMAKE_CURRENT_SOURCE_DIR}/burstbuffer-32KB.txt
        -x ${CMAKE_CURRENT_SOURCE_DIR}/burstbuffer-nodrain-BP4-nsub1.xml -d 1 1 --strong-scaling
)
set_tests_properties(Utils.IOTest.BurstBufferNoDrain.1to1.BP4.Write PROPERTIES PROCESSORS 1)

add_test(NAME Utils.IOTest.BurstBufferNoDrain.1to1.BP4.Write.DumpBB
  COMMAND ${CMAKE_COMMAND}
    -DARG1=-laD
    -DINPUT_FILE=bb/burstbuffer-32KB.bp
    -DOUTPUT_FILE=IOTest.BurstBufferNoDrain.1to1.BP4.Write.bpls-bb.txt
    -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
)

add_test(NAME Utils.IOTest.BurstBufferNoDrain.1to1.BP4.Write.ValidateBB
  COMMAND ${DIFF_COMMAND} -u -w
    ${CMAKE_CURRENT_SOURCE_DIR}/IOTest.BurstBuffer.nproc1.BP4.Write.bpls.txt
    IOTest.BurstBufferNoDrain.1to1.BP4.Write.bpls-bb.txt
)

SetupTestPipeline(
  Utils.IOTest.BurstBufferNoDrain.1to1.BP4
  "Write;Write.DumpBB;Write.ValidateBB"
  TRUE
)

#-------------------------------------------------------------
#  BurstBuffer BP4 streaming tests (write/read simultaniously
#-------------------------------------------------------------


#-------------------------------------------------
#  BurstBuffer BP4 streaming test nproc = 2 
#-------------------------------------------------
if(ADIOS2_RUN_MPI_MPMD_TESTS AND NOT WIN32)
  add_test(NAME Utils.IOTest.BurstBuffer.1to1.BP4.Stream.Run
    COMMAND ${MPIEXEC_COMMAND} 
        ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:adios_iotest>
          -a 1 -c ${CMAKE_CURRENT_SOURCE_DIR}/pipe2.txt
          -x ${CMAKE_CURRENT_SOURCE_DIR}/pipe2-BP4-burstbuffer.xml -d 1 1 --strong-scaling
        :
        ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:adios_iotest>
          -a 2 -c ${CMAKE_CURRENT_SOURCE_DIR}/pipe2.txt
          -x ${CMAKE_CURRENT_SOURCE_DIR}/pipe2-BP4-burstbuffer.xml -d 1 1 --strong-scaling
  
  )
  set_tests_properties(Utils.IOTest.BurstBuffer.1to1.BP4.Stream.Run PROPERTIES PROCESSORS 2)
  
  add_test(NAME Utils.IOTest.BurstBuffer.1to1.BP4.Stream.WriteDump
    COMMAND ${CMAKE_COMMAND}
      -DARG1=-laD
      -DINPUT_FILE=pipe2_write.bp
      -DOUTPUT_FILE=IOTest.BurstBuffer.1to1.BP4.Stream.Write.bpls.txt
      -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
  )
  
  add_test(NAME Utils.IOTest.BurstBuffer.1to1.BP4.Stream.ReadDump
    COMMAND ${CMAKE_COMMAND}
      -DARG1=-laD
      -DINPUT_FILE=pipe2_read.bp
      -DOUTPUT_FILE=IOTest.BurstBuffer.1to1.BP4.Stream.Read.bpls.txt
      -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
  )
  
  add_test(NAME Utils.IOTest.BurstBuffer.1to1.BP4.Stream.WriteValidate
    COMMAND ${DIFF_COMMAND} -u -w
      ${CMAKE_CURRENT_SOURCE_DIR}/IOTest.BurstBuffer.1to1.BP4.Stream.Write.bpls.txt
      IOTest.BurstBuffer.1to1.BP4.Stream.Write.bpls.txt
  )
  
  add_test(NAME Utils.IOTest.BurstBuffer.1to1.BP4.Stream.ReadValidate
    COMMAND ${DIFF_COMMAND} -u -w
      ${CMAKE_CURRENT_SOURCE_DIR}/IOTest.BurstBuffer.1to1.BP4.Stream.Read.bpls.txt
      IOTest.BurstBuffer.1to1.BP4.Stream.Read.bpls.txt
  )
  
  SetupTestPipeline(
    Utils.IOTest.BurstBuffer.1to1.BP4.Stream
    "Run;WriteDump;WriteValidate;ReadDump;ReadValidate"
    TRUE
  )
  
  
  #-------------------------------------------------
  #  BurstBuffer BP4 streaming test nproc = 4 
  #-------------------------------------------------
  add_test(NAME Utils.IOTest.BurstBuffer.3to1.BP4.Stream.Run
    COMMAND ${MPIEXEC_COMMAND} 
        ${MPIEXEC_NUMPROC_FLAG} 3 $<TARGET_FILE:adios_iotest>
          -a 1 -c ${CMAKE_CURRENT_SOURCE_DIR}/pipe2.txt
          -x ${CMAKE_CURRENT_SOURCE_DIR}/pipe2-BP4-burstbuffer-nsub2.xml -d 3 1 --strong-scaling
        :
        ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:adios_iotest>
          -a 2 -c ${CMAKE_CURRENT_SOURCE_DIR}/pipe2.txt
          -x ${CMAKE_CURRENT_SOURCE_DIR}/pipe2-BP4-burstbuffer-nsub2.xml -d 1 1 --strong-scaling
  
  )
  set_tests_properties(Utils.IOTest.BurstBuffer.3to1.BP4.Stream.Run PROPERTIES PROCESSORS 4)
  
  add_test(NAME Utils.IOTest.BurstBuffer.3to1.BP4.Stream.WriteDump
    COMMAND ${CMAKE_COMMAND}
      -DARG1=-laD
      -DINPUT_FILE=pipe2_write.bp
      -DOUTPUT_FILE=IOTest.BurstBuffer.3to1.BP4.Stream.Write.bpls.txt
      -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
  )
  
  add_test(NAME Utils.IOTest.BurstBuffer.3to1.BP4.Stream.ReadDump
    COMMAND ${CMAKE_COMMAND}
      -DARG1=-laD
      -DINPUT_FILE=pipe2_read.bp
      -DOUTPUT_FILE=IOTest.BurstBuffer.3to1.BP4.Stream.Read.bpls.txt
      -P "${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake"
  )
  
  add_test(NAME Utils.IOTest.BurstBuffer.3to1.BP4.Stream.WriteValidate
    COMMAND ${DIFF_COMMAND} -u -w
      ${CMAKE_CURRENT_SOURCE_DIR}/IOTest.BurstBuffer.3to1.BP4.Stream.Write.bpls.txt
      IOTest.BurstBuffer.3to1.BP4.Stream.Write.bpls.txt
  )
  
  add_test(NAME Utils.IOTest.BurstBuffer.3to1.BP4.Stream.ReadValidate
    COMMAND ${DIFF_COMMAND} -u -w
      ${CMAKE_CURRENT_SOURCE_DIR}/IOTest.BurstBuffer.3to1.BP4.Stream.Read.bpls.txt
      IOTest.BurstBuffer.3to1.BP4.Stream.Read.bpls.txt
  )
  
  SetupTestPipeline(
    Utils.IOTest.BurstBuffer.3to1.BP4.Stream
    "Run;WriteDump;WriteValidate;ReadDump;ReadValidate"
    TRUE
  )
endif()

#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

find_package(Threads REQUIRED)

if(ADIOS2_HAVE_SST)
  gtest_add_tests_helper(StagingMPMD MPI_ONLY "" Engine.Staging. ".SST.FFS" EXTRA_ARGS "SST" "MarshalMethod=FFS")
  gtest_add_tests_helper(StagingMPMD MPI_ONLY "" Engine.Staging. ".SST.BP" EXTRA_ARGS "SST" "MarshalMethod=BP")
  gtest_add_tests_helper(OnDemandMPI MPI_ONLY "" Engine.Staging. ".SST" EXTRA_ARGS "SST" "OnDemandContact")
  gtest_add_tests_helper(Threads MPI_NONE "" Engine.Staging. ".SST.FFS" EXTRA_ARGS "SST"  "--engine_params" "MarshalMethod=FFS")
  gtest_add_tests_helper(Threads MPI_NONE "" Engine.Staging. ".SST.BP" EXTRA_ARGS "SST"  "--engine_params" "MarshalMethod=BP")
  gtest_add_tests_helper(Threads MPI_NONE "" Engine.Staging. ".BP4_stream" EXTRA_ARGS "BP4"  "--engine_params" "OpenTimeoutSecs=5")
  gtest_add_tests_helper(Threads MPI_NONE "" Engine.Staging. ".FileStream" EXTRA_ARGS "FileStream")
endif()

set (PROGS TestCommonWrite TestShapeChangingWrite
    TestCommonWriteModes
    TestCommonWriteAttrs
    TestCommonWriteLocal
    TestCommonWriteShared
    TestCommonSpanWrite
    TestDistributionWrite
    TestDefSyncWrite
    TestCommonRead
    TestCommonSpanRead
    TestCommonReadR64
    TestCommonReadLocal
    TestCommonReadShared
    TestCommonReadAttrs
    TestCommonServer
    TestCommonClient
    TestDistributionRead
    TestLocalWrite
    TestLocalRead
    TestStructWrite
    TestStructRead
    TestWriteJoined
    TestReadJoined
    TestOnDemandWrite
    TestOnDemandRead)

foreach(helper ${PROGS} )
  add_executable(${helper} ${helper}.cpp)

  # Workaround for multiple versions of FindSst
  if(SST_INCLUDE_DIRS)
    target_include_directories(${helper} PRIVATE ${SST_INCLUDE_DIRS})
  endif()

  target_link_libraries(${helper} adios2::thirdparty::gtest ${Sst_LIBRARY})
  if(ADIOS2_HAVE_MPI)
    target_link_libraries(${helper} adios2::cxx11_mpi adios2_core_mpi MPI::MPI_C)
  else()
    target_link_libraries(${helper} adios2::cxx11 adios2_core)
  endif()
endforeach()

if(ADIOS2_HAVE_Fortran)
  add_library(TestCommonData_f OBJECT TestData_mod.F90)
  add_executable(TestCommonWrite_f
    TestCommonWriteF.F90
    $<TARGET_OBJECTS:TestCommonData_f>
  )
  add_executable(TestCommonRead_f
    TestCommonReadF.F90
    $<TARGET_OBJECTS:TestCommonData_f>
  )
  target_compile_definitions(TestCommonRead_f PRIVATE
    $<$<BOOL:${ADIOS2_HAVE_FORTRAN_F03_ARGS}>:ADIOS2_HAVE_FORTRAN_F03_ARGS>
    $<$<BOOL:${ADIOS2_HAVE_FORTRAN_GNU_ARGS}>:ADIOS2_HAVE_FORTRAN_GNU_ARGS>
  )
  set_target_properties(TestCommonWrite_f TestCommonRead_f PROPERTIES
    LINKER_LANGUAGE Fortran
  )
  target_compile_definitions(TestCommonWrite_f PRIVATE
    $<$<BOOL:${ADIOS2_HAVE_FORTRAN_F03_ARGS}>:ADIOS2_HAVE_FORTRAN_F03_ARGS>
    $<$<BOOL:${ADIOS2_HAVE_FORTRAN_GNU_ARGS}>:ADIOS2_HAVE_FORTRAN_GNU_ARGS>
  )
  if(ADIOS2_HAVE_MPI)
    target_link_libraries(TestCommonWrite_f adios2_fortran_mpi MPI::MPI_Fortran)
    target_link_libraries(TestCommonRead_f adios2_fortran_mpi MPI::MPI_Fortran)
  else()
    target_link_libraries(TestCommonWrite_f adios2_fortran)
    target_link_libraries(TestCommonRead_f adios2_fortran)
  endif()
endif()

configure_file(
  run_test.py.gen.in
  ${CMAKE_CURRENT_BINARY_DIR}/run_test.py.gen
  @ONLY
)
file(GENERATE
  OUTPUT ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/run_test.py.$<CONFIG>
  INPUT ${CMAKE_CURRENT_BINARY_DIR}/run_test.py.gen)

include (TestSupp.cmake )

if(ADIOS2_HAVE_MPI AND MPIEXEC_EXECUTABLE)
   # we want to know if mpiexec is a shell script or a binary executable
   file(READ ${MPIEXEC_EXECUTABLE} MPIFILEDATA LIMIT 1024 HEX)
   string(LENGTH ${MPIFILEDATA} DATALENGTH)
   math(EXPR last_hex_index "(${DATALENGTH} / 2) - 1")
   foreach(hex_index RANGE ${last_hex_index})
      math (EXPR str_loc "2*${hex_index}")
      string(SUBSTRING "${MPIFILEDATA}" "${str_loc}" "2" char)
      #math (EXPR dec_char "0x${char}" )
      from_hex(${char} dec_char)
      if ((${dec_char} LESS 9) OR (${dec_char} GREATER 126) OR ((${dec_char} GREATER 13) AND (${dec_char} LESS 32)))
          #  message ("$character ${char} is outside of the printable/character ascii range")
	  set (MPIEXEC_IS_BINARY TRUE)
	  break()
      endif()
   endforeach()
endif()

set (SIMPLE_TESTS "1x1;1x1.Attrs;1x1.ModAttrs;1x1DefSync;1x1VarDestruction;1x1.SpanMinMax;1x1.Local2;1x1Struct;1x1Joined")
set (BP5File_ONLY_TESTS "1x1DataWrite;1x1.ShapeChange")

set (SIMPLE_FORTRAN_TESTS "")
if(ADIOS2_HAVE_Fortran)
  set (SIMPLE_FORTRAN_TESTS "FtoC.1x1;CtoF.1x1;FtoF.1x1")
endif()

set (SPECIAL_TESTS "TimeoutReader.1x1;LatestReader.1x1;LatestReaderHold.1x1;DiscardWriter.1x1;1x1.NoPreload;1x1.ForcePreload;1x1LockGeometry")
if (MPIEXEC_IS_BINARY)
    # run_test.py can only kill readers/writers if mpiexec is not a shell script
    list(APPEND SPECIAL_TESTS "KillReadersSerialized.3x2;KillReaders3Max.3x6;KillWriter_2x2;KillWriterTimeout_2x2")
endif()

set (SIMPLE_MPI_TESTS "")
set (SIMPLE_MPI_FORTRAN_TESTS "")
if (ADIOS2_HAVE_MPI)
  set (SIMPLE_MPI_TESTS "2x1;1x2;3x5;5x3;DelayedReader_3x5;3x5LockGeometry;2x1.Local;1x2.Local;3x5.Local;5x3.Local;1x1.LocalVarying;5x3.LocalVarying;2x1ZeroDataVar;2x1ZeroDataR64;2x2.NoData;2x2.HalfNoData;2x1.SharedNothing;2x1.SharedIO;2x1.SharedVar;2x1.SharedNothingSync;2x1.SharedIOSync;2x1.SharedVarSync;3x5EarlyExit;;5x1Joined")
  list (APPEND SPECIAL_TESTS "2x1.NoPreload;2x3.ForcePreload;PreciousTimestep.3x2;PreciousTimestepDiscard.3x2")
  if (ADIOS2_HAVE_Fortran)
    set (SIMPLE_MPI_FORTRAN_TESTS "FtoC.3x5;CtoF.3x5;FtoF.3x5")
  endif()
endif()
 
set (SIMPLE_ZFP_TESTS "")
if (ADIOS2_HAVE_ZFP)
   set (SIMPLE_ZFP_TESTS "ZFPCompression.1x1;ZFPCompression.3x5")
endif()

set (ALL_SIMPLE_TESTS "")
list (APPEND ALL_SIMPLE_TESTS ${SIMPLE_TESTS} ${SIMPLE_FORTRAN_TESTS} ${SIMPLE_MPI_TESTS} ${SIMPLE_ZFP_TESTS})

set (SST_SPECIFIC_TESTS  "")
import_bp_test(WriteMemorySelectionRead 1 1)
list (APPEND SST_SPECIFIC_TESTS  "WriteMemorySelectionRead.1x1")
list (APPEND SST_SPECIFIC_TESTS  "1x1.SstRUDP;1x1.LocalMultiblock;RoundRobinDistribution.1x1x3;AllToAllDistribution.1x1x3;OnDemandSingle.1x1")

if (ADIOS2_HAVE_MPI)
  import_bp_test(WriteMemorySelectionRead 3 3)
  list (APPEND SST_SPECIFIC_TESTS  "WriteMemorySelectionRead.3x3")
  list (APPEND SST_SPECIFIC_TESTS  "2x3.SstRUDP;2x1.LocalMultiblock;5x3.LocalMultiblock;")
endif()

#
#   Setup tests for SST engine
#
SET (BASIC_SST_TESTS "")
if(ADIOS2_HAVE_SST)
  list (APPEND BASIC_SST_TESTS ${ALL_SIMPLE_TESTS} ${SPECIAL_TESTS} ${SST_SPECIFIC_TESTS})
  list (REMOVE_ITEM BASIC_SST_TESTS 1x1DefSync 1x1Flush)
  list (REMOVE_ITEM BASIC_SST_TESTS 1x1.SpanMinMax)
endif()


#  For the moment, only test the default comm pattern (Min)
MutateTestSet( COMM_MIN_SST_TESTS "CommMin" writer "CPCommPattern=Min" "${BASIC_SST_TESTS}" )
#MutateTestSet( COMM_PEER_SST_TESTS "CommPeer" writer "CPCommPattern=Peer" "${BASIC_SST_TESTS}" )

# temporarily remove PreciousTimestep CommPeer tests
list (REMOVE_ITEM COMM_PEER_SST_TESTS "PreciousTimestep")

MutateTestSet( FFS_SST_TESTS "FFS" writer "MarshalMethod=FFS" "${COMM_MIN_SST_TESTS};${COMM_PEER_SST_TESTS}" )
MutateTestSet( BP5_SST_TESTS "BP5" writer "MarshalMethod=BP5" "${COMM_MIN_SST_TESTS};${COMM_PEER_SST_TESTS}" )
MutateTestSet( BP_SST_TESTS "BP" writer "MarshalMethod=BP" "${COMM_MIN_SST_TESTS};${COMM_PEER_SST_TESTS}" )

# no SSE engine does Joined
list (FILTER BP_SST_TESTS EXCLUDE REGEX "Joined*")

set (SST_TESTS "")
LIST (APPEND SST_TESTS ${BP5_SST_TESTS} ${BP_SST_TESTS})

# Zero Data tests are unreliable with SST and BP marshaling
list (FILTER SST_TESTS EXCLUDE REGEX "2x1ZeroData.*BP")

# BP Marshalling doesn't support Struct data
list (FILTER SST_TESTS EXCLUDE REGEX "1x1Struct.*BP$")

foreach(test ${SST_TESTS})
    add_common_test(${test} SST)
endforeach()

if(ADIOS2_HAVE_MPI AND ADIOS2_RUN_MPI_MPMD_TESTS AND NOT MSVC)
    set (SSC_BASE_TESTS "1x1;1x1.Attrs;2x1;1x2;2x1ZeroDataVar;2x1ZeroDataR64;3x5;5x3;3x5LockGeometry")
    set (SSC_TESTS  ${SIMPLE_FORTRAN_TESTS} ${SIMPLE_MPI_FORTRAN_TESTS} ${SSC_BASE_TESTS})
    foreach(test ${SSC_TESTS})
        add_common_test(${test} SSC)
    endforeach()
endif()

#
#   Setup tests for BP engines
#

if(NOT WIN32)    # not on windows
    set (BP_TESTS ${ALL_SIMPLE_TESTS})
    # Delayed reader not worth testing on file engines
    list (FILTER BP_TESTS EXCLUDE REGEX "DelayedReader")
    # SharedVars fail with BP*
    list (FILTER BP_TESTS EXCLUDE REGEX ".*SharedVar$")
    # The nobody-writes-data-in-a-timestep tests don't work for any BP file engine
    list (FILTER BP_TESTS EXCLUDE REGEX ".*NoData$")
    # BP3 and BP4 don't do PerformDataWrite
    list (FILTER BP_TESTS EXCLUDE REGEX ".*DataWrite.*")
    # BP3 and BP4 don't do Modifiable Attributes
    list (FILTER BP_TESTS EXCLUDE REGEX ".*ModAttr.*")
    # BP3 and BP4 don't do VarDestruction
    list (FILTER BP_TESTS EXCLUDE REGEX ".*VarDestruction.*")
    # BP3 and BP4 don't do structs
    list (FILTER BP_TESTS EXCLUDE REGEX ".*1x1Struct.*")
    foreach(test ${BP_TESTS})
        add_common_test(${test} BP4)
    endforeach()
    # BP3 fails to handle lots of timesteps
    list (FILTER BP_TESTS EXCLUDE REGEX "Bulk.*")
    # BP3 doesn't do Joined
    list (FILTER BP_TESTS EXCLUDE REGEX "Joined.*")
    foreach(test ${BP_TESTS})
        add_common_test(${test} BP3)
    endforeach()
endif()

#  BP5 tests
set (BP5_TESTS ${ALL_SIMPLE_TESTS} ${SPECIAL_TESTS} ${BP5File_ONLY_TESTS})
# Delayed reader not worth testing on file engines
list (FILTER BP5_TESTS EXCLUDE REGEX "DelayedReader")
# Discard not a feature of BP5
list (FILTER BP5_TESTS EXCLUDE REGEX ".*DiscardWriter.1x1")
# PreciousTimestep not a feature of BP5
list (FILTER BP5_TESTS EXCLUDE REGEX ".*PreciousTimestep")
# LatestTimestep not a feature of BP5
list (FILTER BP5_TESTS EXCLUDE REGEX ".*LatestReader")
# KillWriter fails with BP5
list (FILTER BP5_TESTS EXCLUDE REGEX ".*KillWriter")
# KillReaders  We swear this isn't necessary for BP5 streaming
list (FILTER BP5_TESTS EXCLUDE REGEX ".*KillReaders")
foreach(test ${BP5_TESTS})
    add_common_test(${test} BP5)
endforeach()


#
#   Setup streaming tests for BP4 engine
#
if(NOT MSVC)    # not on windows
    # BP4 streaming tests start with all the simple tests, but with a timeout added on open
    LIST (APPEND BP4_STREAM_TESTS ${ALL_SIMPLE_TESTS} ${SPECIAL_TESTS})
    MutateTestSet( BP4_STREAM_TESTS "BPS" reader "OpenTimeoutSecs=10,BeginStepPollingFrequencySecs=0.1" "${BP4_STREAM_TESTS}")
    # no BP4 streaming doesn't do Joined?
    list (FILTER BP4_STREAM_TESTS EXCLUDE REGEX "Joined*")
    # SharedVars fail with BP4_streaming*
    list (FILTER BP4_STREAM_TESTS EXCLUDE REGEX ".*SharedVar.BPS$")
   # Discard not a feature of BP4
    list (FILTER BP4_STREAM_TESTS EXCLUDE REGEX ".*DiscardWriter.1x1.*BPS$")
    # PreciousTimestep not a feature of BP4
    list (FILTER BP4_STREAM_TESTS EXCLUDE REGEX ".*Precious.*BPS$")
    # Timeout on Reader BeginStep failing in BP4
    list (FILTER BP4_STREAM_TESTS EXCLUDE REGEX ".*TimeoutReader.1x1.*BPS$")
    # LatestTimestep not a feature of BP4
    list (FILTER BP4_STREAM_TESTS EXCLUDE REGEX ".*LatestReader.*BPS$")
    # Preload not a feature of BP4_streaming
    list (FILTER BP4_STREAM_TESTS EXCLUDE REGEX ".*Preload.*BPS$")
    # KillWriter fails with BP4
    list (FILTER BP4_STREAM_TESTS EXCLUDE REGEX ".*KillWriter.*BPS$")
    # KillReaders  We swear this isn't necessary for BP4 streaming
    list (FILTER BP4_STREAM_TESTS EXCLUDE REGEX ".*KillReaders.*BPS$") 
    # SharedVars fail with BP4_streaming*
    list (FILTER BP4_STREAM_TESTS EXCLUDE REGEX ".*SharedVar.BPS$")
    # Local fail with BP4_streaming*
    list (FILTER BP4_STREAM_TESTS EXCLUDE REGEX ".*Local.BPS$")
    # BP3 and BP4 don't do VarDestruction
    list (FILTER BP4_STREAM_TESTS EXCLUDE REGEX ".*VarDestruction.*")
    # Fortran scalar reads don't work for BP4 streaming
    list (FILTER BP4_STREAM_TESTS EXCLUDE REGEX "(FtoF|CtoF).*")
    # Local Vars have a heisen failure for BP4 streaming
#list (FILTER BP4_STREAM_TESTS EXCLUDE REGEX ".*LocalVarying.BPS$")
    # The nobody-writes-data-in-a-timestep tests don't work for any BP file engine
    list (FILTER BP4_STREAM_TESTS EXCLUDE REGEX ".*NoData.BPS$")
    list (FILTER BP4_STREAM_TESTS EXCLUDE REGEX ".*1x1DefSync.*")
    list (FILTER BP4_STREAM_TESTS EXCLUDE REGEX ".*Flush.*")
    # BP3 and BP4 don't do Modifiable Attributes
    list (FILTER BP4_STREAM_TESTS EXCLUDE REGEX ".*ModAttr.*")
    # BP3 and BP4 don't do structs
    list (FILTER BP4_STREAM_TESTS EXCLUDE REGEX ".*1x1Struct.*")
    
    foreach(test ${BP4_STREAM_TESTS})
        add_common_test(${test} BP4_stream)
    endforeach()

    # add burst buffer tests
    MutateTestSet( BP4_BBSTREAM_TESTS "BB" writer "BurstBufferPath=bb,BurstBufferVerbose=2" "${BP4_STREAM_TESTS}")
    foreach(test ${BP4_BBSTREAM_TESTS})
        add_common_test(${test} BP4_stream)
    endforeach()
    
endif()

#
#   Setup streaming tests for FileStream virtual engine (BP4+StreamReader=true)
#
if(NOT MSVC)    # not on windows
    # FileStream streaming tests start with all the simple tests, but with a timeout added on open
    LIST (APPEND FILESTREAM_TESTS ${SIMPLE_TESTS} ${SIMPLE_MPI_TESTS})
    MutateTestSet( FILESTREAM_TESTS "FS" reader "OpenTimeoutSecs=10,BeginStepPollingFrequencySecs=0.1" "${FILESTREAM_TESTS}")
    # SharedVars fail with file_streaming*
    list (FILTER FILESTREAM_TESTS EXCLUDE REGEX ".*SharedVar.FS$")
    # SharedVars fail with file_streaming*
    list (FILTER FILESTREAM_TESTS EXCLUDE REGEX ".*SharedVar.FS$")
    # Local fail with file_streaming*
    list (FILTER FILESTREAM_TESTS EXCLUDE REGEX ".*Local.FS$")
    # The nobody-writes-data-in-a-timestep tests don't work for any BP-file based engine
    list (FILTER FILESTREAM_TESTS EXCLUDE REGEX ".*NoData.FS$")
    # Don't need to repeat tests that are identical for BP4 and FileStream
    list (FILTER FILESTREAM_TESTS EXCLUDE REGEX ".*NoReaderNoWait.FS$")
    list (FILTER FILESTREAM_TESTS EXCLUDE REGEX ".*TimeoutOnOpen.FS$")
    list (FILTER FILESTREAM_TESTS EXCLUDE REGEX ".*NoReaderNoWait.FS$")
    list (FILTER FILESTREAM_TESTS EXCLUDE REGEX ".*1x1DefSync.*")
    list (FILTER FILESTREAM_TESTS EXCLUDE REGEX ".*Flush.*")
    # BP3 and BP4 don't do VarDestruction
    list (FILTER FILESTREAM_TESTS EXCLUDE REGEX ".*VarDestruction.*")
    # BP3 and BP4 don't do structs
    list (FILTER FILESTREAM_TESTS EXCLUDE REGEX ".*1x1Struct.*")
    # no FileStream doesn't do Joined
    list (FILTER FILESTREAM_TESTS EXCLUDE REGEX "Joined*")

    foreach(test ${FILESTREAM_TESTS})
        add_common_test(${test} FileStream)
    endforeach()

    MutateTestSet( FileStream_BBSTREAM_TESTS "BB" writer "BurstBufferPath=bb,BurstBufferVerbose=2" "${FILESTREAM_TESTS}")
    foreach(test ${FileStream_BBSTREAM_TESTS})
        add_common_test(${test} FileStream)
    endforeach()

endif()

#
#   Setup tests for HDF5 engine
#
if(NOT MSVC AND ADIOS2_HAVE_HDF5 AND (NOT ADIOS2_HAVE_MPI OR HDF5_IS_PARALLEL))    # not on windows and only if we have HDF5
    set (HDF5_TESTS ${ALL_SIMPLE_TESTS})
    # Delayed reader not worth testing on file engines
    list (FILTER HDF5_TESTS EXCLUDE REGEX "DelayedReader")
    # Local variables don't work on hdf5 
    list (FILTER HDF5_TESTS EXCLUDE REGEX ".*Local")
    # Oddly, some multi-writer tests fail
    list (FILTER HDF5_TESTS EXCLUDE REGEX "[2345]x")
    # The nobody-writes-data-in-a-timestep tests don't work with HDF5
    list (FILTER HDF5_TESTS EXCLUDE REGEX ".*HalfNoData$")
    # The variables with HDF5 reads persist
    list (FILTER HDF5_TESTS EXCLUDE REGEX ".*VarDestruction.*")
    # The No span supoprt with HDF5
    list (FILTER HDF5_TESTS EXCLUDE REGEX ".*SpanMinMax.*")
    # HDF5 doesn't do Modifiable Attributes
    list (FILTER HDF5_TESTS EXCLUDE REGEX ".*ModAttr.*")
    # HDF5 engine doesn't do structs (but HDF5 does)
    list (FILTER HDF5_TESTS EXCLUDE REGEX ".*1x1Struct.*")
    # HDF5 doesn't do Joined
    list (FILTER HDF5_TESTS EXCLUDE REGEX "Joined*")

    foreach(test ${HDF5_TESTS})
        add_common_test(${test} HDF5)
    endforeach()
endif()

# Tests particularities
if (WIN32)
    set_property(TEST Staging.1x1DefSync.BP5 PROPERTY TIMEOUT 500)
    set_property(TEST Staging.1x1DataWrite.BP5 PROPERTY TIMEOUT 500)
endif()

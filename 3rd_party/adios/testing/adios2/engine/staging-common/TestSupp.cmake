#
#
# This Cmake testing specification is unfortunately complex, and will
# likely be hard to maintain, but at the moment I don't know a better
# way to do it (contributions welcome).  The idea is that we have a set of tests that we'd
# like to run on staging engines.  Some tests may be applicable to
# some engines and not to others.  Some tests we want to run multiple
# times, but with different engine parameters set to test different
# modes of operation for the engines.
#
# The approach here is as follows:
#   - A test is specified by setting a CMake variable of the form
#     <testname>_CMD.  
#   - If the test requires a different time out than is usually
#     specified (30 seconds), then you can also define a variable of
#     the form <testname>_TIMEOUT to control the timeout.
#   - If the test requires special properties you can also define a variable of 
#     the form <testname>_PROPERTIES.
# 
#     BASIC TEST ADDITION
# Given the existence of the above variables, the cmake function
# "add_common_test(testname engine)" adds a CTest which runs the given
# test with a particular ADIOS2 engine, adding the engine's name to
# the test name to differentiate it from other tests.  For example,
# given the variable 1x1_CMD, set to
#  "run_staging_test -nw 1 -nr 1 -v -p TestCommon"
# and 1x1_PROPERTIES set to "RUN_SERIAL;1" then
# add_common_test(1x1 SST) ends up doing:
# 
#  add_test(NAME "Staging.1x1.SST"
#          COMMAND "run_staging_test -e SST -f Staging.1x1.SST -nw 1 -nr 1")
#  set_tests_properties(${testname} PROPERTIES TIMEOUT 30 RUN_SERIAL 1)
#
#    RUNNING TESTS WITH DIFFERENT ENGINE PARAMETERS
# 
# One of the design goals here is to enable running the same tests
# with different engine parameters, *_CMD strings can also contain the
# string "WENGINE_PARAMS".  This string is treated specially by the
# function MutateTestSet().  This function is takes a list of tests

# (I.E. things with _CMD strings defined like above) and creates a new
# set of tests where a specified engine parameter gets added to the in
# the location of the WENGINE_PARAMS string.  MutateTestSet takes 4 parameters:
# output_test_list, param_name, param_spec, and input test list.  For example
# MutateTestSet( COMM_MIN_SST_TESTS "CommMin" "CPCommPattern=Min" "${BASIC_SST_TESTS}" )
# If BASIC_SST_TESTS contains "1x1" as defined above, MutateTestSet
# will add the test "1x1.CommMin", by defining the variable
# 1x1.CommMin_CMD, using the original value of 1x1_CMD buth with
# "CPCommPattern=Min: added to the WENGINE_PARAMS location (if
# present).  Any 1x1_TIMEOUT and 1x1_PROPERTIES values will also be
# propogated to 1x1.CommMin_TIMEOUT and 1x1.CommMin_PROPERTIES.
# "1x1.CommMin" will also be added to the output test list.
# 
# Note that MutateTestSet() can be used multiple times to add multiple
# engine params.  The WENGINE_PARAMS string is retained in the
# resulting _CMD strings until add_common_test() which removes it.
#

find_package(Python3 REQUIRED)

# Change the STAGING_COMMON_TEST_SUPP_VERBOSE value to ON for debugging output
#
set (STAGING_COMMON_TEST_SUPP_VERBOSE OFF)

set (1x1_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1")
set (1x1Struct_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1 -r $<TARGET_FILE:TestStructRead> -w $<TARGET_FILE:TestStructWrite> ")
set (1x1Joined_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1 -r $<TARGET_FILE:TestReadJoined> -w $<TARGET_FILE:TestWriteJoined> ")
set (1x1.ShapeChange_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1 -w $<TARGET_FILE:TestShapeChangingWrite>")
set (1x1GetSync_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1 --rarg=--read_mode --rarg=sync")
set (1x1DontCloseWriter_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1 --warg=--dont_close")
set (1x1DontCloseReader_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1 --rarg=--dont_close")
set (1x1DefSync_CMD "TestDefSyncWrite --data_size 200 --engine_params ChunkSize=500,MinDeferredSize=150")
set (1x1DefSync_TIMEOUT 180)
set (1x1VarDestruction_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1 --rarg=--var_destruction")
set (1x1DataWrite_CMD "TestDefSyncWrite --flush --data_size 200 --engine_params ChunkSize=500,MinDeferredSize=150")
set (1x1DataWrite_TIMEOUT 180)
set (1x1.NoPreload_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1 --rarg=PreloadMode=SstPreloadNone,RENGINE_PARAMS")
set (1x1.SstRUDP_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1 --rarg=DataTransport=WAN,WANDataTransport=enet,RENGINE_PARAMS --warg=DataTransport=WAN,WANDataTransport=enet,WENGINE_PARAMS")
set (1x1.NoData_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1 --warg=--no_data --rarg=--no_data")
set (2x2.NoData_CMD "run_test.py.$<CONFIG> -nw 2 -nr 2 --warg=--no_data --rarg=--no_data")
set (2x2.HalfNoData_CMD "run_test.py.$<CONFIG> -nw 2 -nr 2 --warg=--no_data --warg=--no_data_node --warg=1 --rarg=--no_data --rarg=--no_data_node --rarg=1" )
set (1x1.ForcePreload_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1 --rarg=PreloadMode=SstPreloadOn,RENGINE_PARAMS")
set (1x1Bulk_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1 --warg=--nx --warg=10000 --warg=--num_steps --warg=101 --rarg=--num_steps --rarg=101")
set (1x1LockGeometry_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1  --warg=--num_steps --warg=101  --warg=--nx --warg=50 --rarg=--num_steps --rarg=101 --warg=--lock_geometry --rarg=--lock_geometry --rarg=PreloadMode=SstPreloadNone,RENGINE_PARAMS")
set (2x1_CMD "run_test.py.$<CONFIG> -nw 2 -nr 1")
set (2x1ZeroDataVar_CMD "run_test.py.$<CONFIG> -nw 2 -nr 1 --warg=--zero_data_var")
set (2x1ZeroDataR64_CMD "run_test.py.$<CONFIG> -nw 2 -nr 1  -r $<TARGET_FILE:TestCommonReadR64> --warg=--zero_data_var")
set (2x1.NoPreload_CMD "run_test.py.$<CONFIG> -nw 2 -nr 1 --rarg=PreloadMode=SstPreloadNone,RENGINE_PARAMS")
set (2x3.ForcePreload_CMD "run_test.py.$<CONFIG> -nw 2 -nr 3 --rarg=PreloadMode=SstPreloadOn,RENGINE_PARAMS")
set (2x3.SstRUDP_CMD "run_test.py.$<CONFIG> -nw 2 -nr 3 --rarg=DataTransport=WAN,WANDataTransport=enet,RENGINE_PARAMS --warg=DataTransport=WAN,WANDataTransport=enet,WENGINE_PARAMS")
set (1x2_CMD "run_test.py.$<CONFIG> -nw 1 -nr 2")
set (3x5_CMD "run_test.py.$<CONFIG> -nw 3 -nr 5")
set (3x5LockGeometry_CMD "run_test.py.$<CONFIG> -nw 3 -nr 5 --warg=--num_steps --warg=50 --warg=--ms_delay --warg=10 --rarg=--num_steps --rarg=50 --warg=--lock_geometry --rarg=--lock_geometry")
set (1x1EarlyExit_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1 --warg=--num_steps --warg=50 --rarg=--num_steps --rarg=5 --rarg=--early_exit")
set (3x5EarlyExit_CMD "run_test.py.$<CONFIG> -nw 3 -nr 5 --warg=--num_steps --warg=50 --rarg=--num_steps --rarg=5 --rarg=--early_exit")
set (3x5LockGeometry_TIMEOUT 60)
set (5x1Joined_CMD "run_test.py.$<CONFIG> -nw 5 -nr 1 -r $<TARGET_FILE:TestReadJoined> -w $<TARGET_FILE:TestWriteJoined> ")
set (5x3_CMD "run_test.py.$<CONFIG> -nw 5 -nr 3")
set (5x3DontClose_CMD "run_test.py.$<CONFIG> -nw 5 -nr 3 --rarg=--dont_close --warg=--dont_close")
set (1x1.Local_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1  -w $<TARGET_FILE:TestCommonWriteLocal> -r $<TARGET_FILE:TestCommonReadLocal>")
set (1x1.Local2_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1  -w $<TARGET_FILE:TestLocalWrite> -r $<TARGET_FILE:TestLocalRead>")
set (1x1.SpanMinMax_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1  -w $<TARGET_FILE:TestCommonSpanWrite> -r $<TARGET_FILE:TestCommonSpanRead>")
set (2x1.Local_CMD "run_test.py.$<CONFIG> -nw 2 -nr 1  -w $<TARGET_FILE:TestCommonWriteLocal> -r $<TARGET_FILE:TestCommonReadLocal>")
set (1x2.Local_CMD "run_test.py.$<CONFIG> -nw 1 -nr 2  -w $<TARGET_FILE:TestCommonWriteLocal> -r $<TARGET_FILE:TestCommonReadLocal>")
set (3x5.Local_CMD "run_test.py.$<CONFIG> -nw 3 -nr 5  -w $<TARGET_FILE:TestCommonWriteLocal> -r $<TARGET_FILE:TestCommonReadLocal>")
set (5x3.Local_CMD "run_test.py.$<CONFIG> -nw 5 -nr 3  -w $<TARGET_FILE:TestCommonWriteLocal> -r $<TARGET_FILE:TestCommonReadLocal>")
set (1x1.LocalVarying_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1  -w $<TARGET_FILE:TestCommonWriteLocal> --warg=--nx --warg=20 --warg=--varying_data_size -r $<TARGET_FILE:TestCommonReadLocal> --rarg=--nx --rarg=20 --rarg=--varying_data_size ")
set (5x3.LocalVarying_CMD "run_test.py.$<CONFIG> -nw 5 -nr 3  -w $<TARGET_FILE:TestCommonWriteLocal> --warg=--nx --warg=20 --warg=--varying_data_size -r $<TARGET_FILE:TestCommonReadLocal> --rarg=--nx --rarg=20 --rarg=--varying_data_size ")
set (1x1.LocalMultiblock_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1  -w $<TARGET_FILE:TestCommonWriteLocal> -r $<TARGET_FILE:TestCommonReadLocal> --warg=--local_count --warg=2 --rarg=--local_count --rarg=2")
set (2x1.LocalMultiblock_CMD "run_test.py.$<CONFIG> -nw 2 -nr 1  -w $<TARGET_FILE:TestCommonWriteLocal> -r $<TARGET_FILE:TestCommonReadLocal> --warg=--local_count --warg=3 --rarg=--local_count --rarg=3")
set (5x3.LocalMultiblock_CMD "run_test.py.$<CONFIG> -nw 5 -nr 3  -w $<TARGET_FILE:TestCommonWriteLocal> -r $<TARGET_FILE:TestCommonReadLocal> --warg=--local_count --warg=5 --rarg=--local_count --rarg=5")
set (DelayedReader_3x5_CMD "run_test.py.$<CONFIG> -rd 5 -nw 3 -nr 5")
set (FtoC.3x5_CMD "run_test.py.$<CONFIG> -nw 3 -nr 5  -w $<TARGET_FILE:TestCommonWrite_f>")
set (FtoF.3x5_CMD "run_test.py.$<CONFIG> -nw 3 -nr 5  -w $<TARGET_FILE:TestCommonWrite_f> -r $<TARGET_FILE:TestCommonRead_f>")
set (1x1.SharedNothing_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1  -w $<TARGET_FILE:TestCommonWriteShared> -r $<TARGET_FILE:TestCommonReadShared> --warg=--write_mode --warg=deferred")
set (1x1.SharedIO_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1  -w $<TARGET_FILE:TestCommonWriteShared> -r $<TARGET_FILE:TestCommonReadShared> --warg=--shared_io --rarg=--shared_io --warg=--write_mode --warg=deferred")
set (1x1.SharedVar_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1  -w $<TARGET_FILE:TestCommonWriteShared> -r $<TARGET_FILE:TestCommonReadShared> --warg=--shared_var --rarg=--shared_var --warg=--write_mode --warg=deferred")
set (1x1.SharedNothingSync_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1  -w $<TARGET_FILE:TestCommonWriteShared> -r $<TARGET_FILE:TestCommonReadShared> --warg=--write_mode --warg=sync")
set (1x1.SharedIOSync_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1  -w $<TARGET_FILE:TestCommonWriteShared> -r $<TARGET_FILE:TestCommonReadShared> --warg=--shared_io --rarg=--shared_io --warg=--write_mode --warg=sync")
set (1x1.SharedVarSync_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1  -w $<TARGET_FILE:TestCommonWriteShared> -r $<TARGET_FILE:TestCommonReadShared> --warg=--shared_var --rarg=--shared_var --warg=--write_mode --warg=sync")

set (2x1.SharedNothing_CMD "run_test.py.$<CONFIG> -nw 2 -nr 1  -w $<TARGET_FILE:TestCommonWriteShared> -r $<TARGET_FILE:TestCommonReadShared> --warg=--write_mode --warg=deferred")
set (2x1.SharedIO_CMD "run_test.py.$<CONFIG> -nw 2 -nr 1  -w $<TARGET_FILE:TestCommonWriteShared> -r $<TARGET_FILE:TestCommonReadShared> --warg=--shared_io --rarg=--shared_io --warg=--write_mode --warg=deferred")
set (2x1.SharedVar_CMD "run_test.py.$<CONFIG> -nw 2 -nr 1  -w $<TARGET_FILE:TestCommonWriteShared> -r $<TARGET_FILE:TestCommonReadShared> --warg=--shared_var --rarg=--shared_var --warg=--write_mode --warg=deferred")
set (2x1.SharedNothingSync_CMD "run_test.py.$<CONFIG> -nw 2 -nr 1  -w $<TARGET_FILE:TestCommonWriteShared> -r $<TARGET_FILE:TestCommonReadShared> --warg=--write_mode --warg=sync")
set (2x1.SharedIOSync_CMD "run_test.py.$<CONFIG> -nw 2 -nr 1  -w $<TARGET_FILE:TestCommonWriteShared> -r $<TARGET_FILE:TestCommonReadShared> --warg=--shared_io --rarg=--shared_io --warg=--write_mode --warg=sync")
set (2x1.SharedVarSync_CMD "run_test.py.$<CONFIG> -nw 2 -nr 1  -w $<TARGET_FILE:TestCommonWriteShared> -r $<TARGET_FILE:TestCommonReadShared> --warg=--shared_var --rarg=--shared_var --warg=--write_mode --warg=sync")


# NoReaderNoWait runs a writer with the RendezvousReaderCount = 0 and then never spawns a reader.  The test should run to termination and execute cleanly
set (NoReaderNoWait_CMD "run_test.py.$<CONFIG> -nw 1 -nr 0 --warg=RendezvousReaderCount=0,QueueLimit=3,QueueFullPolicy=discard,WENGINE_PARAMS")

# TimeoutOnOpen runs a reader but never spawns a writer.  The Open should run to timeout and throw an exception
set (TimeoutOnOpen_CMD "run_test.py.$<CONFIG> -nw 0 -nr 1 --rarg=--expect_timeout --rarg=OpenTimeoutSecs=10,RENGINE_PARAMS")

# The Modes test checks to see that we can mix and match Put Sync and Deferred modes and still get good data
set (1x1.Modes_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1  -w $<TARGET_FILE:TestCommonWriteModes>")

# 1x1.Attrs tests writing and reading of attributes defined before Open
set (1x1.Attrs_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1  -w $<TARGET_FILE:TestCommonWriteAttrs> -r $<TARGET_FILE:TestCommonReadAttrs>")
set (1x1.ModAttrs_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1  -w $<TARGET_FILE:TestCommonWriteAttrs> -r $<TARGET_FILE:TestCommonReadAttrs> --rarg=--modifiable_attributes --warg=--modifiable_attributes")

# Basic Fortran tests, Fortran to C, C to Fortran and Fortran to Fortran
set (FtoC.1x1_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1  -w $<TARGET_FILE:TestCommonWrite_f>")
set (CtoF.1x1_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1  -r $<TARGET_FILE:TestCommonRead_f>")
set (FtoF.1x1_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1  -w $<TARGET_FILE:TestCommonWrite_f> -r $<TARGET_FILE:TestCommonRead_f>")

# MPI Fortran tests, Fortran to C, C to Fortran and Fortran to Fortran
set (FtoC.3x5_CMD "run_test.py.$<CONFIG> -nw 3 -nr 5  -w $<TARGET_FILE:TestCommonWrite_f>")
set (CtoF.3x5_CMD "run_test.py.$<CONFIG> -nw 3 -nr 5  -r $<TARGET_FILE:TestCommonRead_f>")
set (FtoF.3x5_CMD "run_test.py.$<CONFIG> -nw 3 -nr 5  -w $<TARGET_FILE:TestCommonWrite_f> -r $<TARGET_FILE:TestCommonRead_f>")

# Tests for ZFP compression (where supported by an engine param)
set (ZFPCompression.1x1_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1 --warg=CompressionMethod=zfp,WENGINE_PARAMS" )
set (ZFPCompression.3x5_CMD "run_test.py.$<CONFIG> -nw 3 -nr 5 --warg=CompressionMethod=zfp,WENGINE_PARAMS" )

# Test if writer will survive readers departing unexpectedly
set (KillReadersSerialized.3x2_CMD "run_test.py.$<CONFIG> --test_protocol kill_readers  -nw 3 -nr 2 --max_readers 1 --warg=RendezvousReaderCount=0,WENGINE_PARAMS --rarg=--ignore_time_gap")
set (KillReadersSerialized.3x2_TIMEOUT "300")
set (KillReaders3Max.3x6_CMD "run_test.py.$<CONFIG> --test_protocol kill_readers  -nw 3 -nr 2 --max_readers 3 --warg=RendezvousReaderCount=0,WENGINE_PARAMS --rarg=--ignore_time_gap")
set (KillReaders3Max.3x6_TIMEOUT "300")

set (KillWriter_2x2_CMD "run_test.py.$<CONFIG> --test_protocol kill_writer   -nw 2 -nr 2 --interval 10 --warg=RendezvousReaderCount=1,WENGINE_PARAMS --rarg=--expect_writer_failure --rarg=--num_steps --rarg=1000")
set (KillWriterTimeout_2x2_CMD "run_test.py.$<CONFIG> --test_protocol kill_writer -nw 2 -nr 2 --interval 10 --warg=RendezvousReaderCount=1,WENGINE_PARAMS --rarg=--expect_writer_failure --rarg=--num_steps --rarg=1000 --rarg=--non_blocking")

set (PreciousTimestep.3x2_CMD "run_test.py.$<CONFIG> --test_protocol kill_readers  -nw 3 -nr 2 --max_readers 2 --warg=FirstTimestepPrecious=True,RendezvousReaderCount=0,WENGINE_PARAMS --rarg=--ignore_time_gap --rarg=--precious_first")

set (PreciousTimestep.3x2_TIMEOUT "300")

set (PreciousTimestepDiscard.3x2_CMD "run_test.py.$<CONFIG> --test_protocol kill_readers  -nw 3 -nr 2 --max_readers 2 --warg=FirstTimestepPrecious=On,RendezvousReaderCount=0,QueueLimit=3,QueueFullPolicy=discard,WENGINE_PARAMS --rarg=--ignore_time_gap --rarg=--precious_first --rarg=--discard --warg=--ms_delay --warg=500")
set (PreciousTimestepDiscard.3x2_TIMEOUT "300")

# Writer StepDistributionModes.  Here we run the writer and three clients
set (AllToAllDistribution.1x1x3_CMD "run_test.py.$<CONFIG> --test_protocol multi_client -nw 1 -nr 1 -w $<TARGET_FILE:TestDistributionWrite> -r $<TARGET_FILE:TestDistributionRead> --warg=RendezvousReaderCount=3,WENGINE_PARAMS")
set (RoundRobinDistribution.1x1x3_CMD "run_test.py.$<CONFIG> --test_protocol multi_client -nw 1 -nr 1 -w $<TARGET_FILE:TestDistributionWrite> -r $<TARGET_FILE:TestDistributionRead> --warg=RendezvousReaderCount=3,WENGINE_PARAMS --warg=--round_robin --rarg=--round_robin")
set (OnDemandSingle.1x1_CMD "run_test.py.$<CONFIG> -w $<TARGET_FILE:TestOnDemandWrite> -r $<TARGET_FILE:TestOnDemandRead>")
set (OnDemandDistribution.1x1x3_CMD "run_test.py.$<CONFIG> --test_protocol multi_client -nw 1 -nr 1 -w $<TARGET_FILE:TestDistributionWrite> -r $<TARGET_FILE:TestDistributionRead> --warg=RendezvousReaderCount=3,WENGINE_PARAMS --warg=--on_demand --rarg=--on_demand --warg=--num_steps --warg=20")

# Readers using BeginStep with timeout.  Here we run the writer with a longer delay to make the reader timeout
set (TimeoutReader.1x1_CMD "run_test.py.$<CONFIG> --test_protocol one_client -nw 1 -nr 1 --rarg=--non_blocking --warg=--ms_delay --warg=2000")
set (TimeoutReader.1x1_TIMEOUT "60")

# Readers using LatestAvailable   Writer runs faster than reader, so we expect misses
set (LatestReader.1x1_CMD "run_test.py.$<CONFIG> --test_protocol one_client -nw 1 -nr 1 --warg=--ms_delay --warg=250 --rarg=--latest --rarg=--long_first_delay")

set (LatestReaderHold.1x1_CMD "run_test.py.$<CONFIG> --test_protocol one_client -nw 1 -nr 1 --warg=--ms_delay --warg=250 --rarg=--latest --rarg=--long_first_delay --rarg=--delay_while_holding")

# A faster writer and a queue policy that will cause timesteps to be discarded
set (DiscardWriter.1x1_CMD "run_test.py.$<CONFIG> --test_protocol one_client -nw 1 -nr 1 --warg=--engine_params --warg=QueueLimit=1,QueueFullPolicy=discard,WENGINE_PARAMS --warg=--ms_delay --warg=250 --rarg=--discard")

# Readers using Advancing attributes
set (CumulativeAttr.1x1_CMD "run_test.py.$<CONFIG> -nw 1 -nr 1 --warg=--advancing_attributes --rarg=--advancing_attributes")

function(remove_engine_params_placeholder dst_str src_str )
    string(REGEX REPLACE "([^ 		  ]*),WENGINE_PARAMS" "\\1" src_str "${src_str}")
    string(REGEX REPLACE "([^ 		  ]*),RENGINE_PARAMS" "\\1" src_str "${src_str}")
    if ("${src_str}" MATCHES "ENGINE_PARAMS")
       # empty engine params remains
       string(REGEX REPLACE "--warg=--engine_params --warg=WENGINE_PARAMS" "" src_str "${src_str}")       
       string(REGEX REPLACE "--rarg=--engine_params --rarg=RENGINE_PARAMS" "" src_str "${src_str}")       
       string(REGEX REPLACE "--warg=ENGINE_PARAMS" "" src_str "${src_str}")       
       string(REGEX REPLACE "--rarg=ENGINE_PARAMS" "" src_str "${src_str}")       
    endif()
  set(${dst_str} ${src_str} PARENT_SCOPE)
endfunction()

function(add_engine_param participant dst_string src_string engine_param)
  if (participant STREQUAL "writer") 
      set(ARG "--warg")
      set(L "W")
  elseif (participant STREQUAL "reader") 
      set(ARG "--rarg")
      set(L "R")
  else ()
      message(WARNING "add_engine_param participant value not valid \"${participant}\", assuming writer")
      set(ARG "--warg")
      set(L "W")
  endif ()
  set (tmp_string "")
  if ("${src_string}" MATCHES ",${L}ENGINE_PARAMS")
      string(REGEX REPLACE ",${L}ENGINE_PARAMS" ",${engine_param},${L}ENGINE_PARAMS" tmp_string ${src_string})
  elseif ("${src_string}" MATCHES "${L}ENGINE_PARAMS")
      string(REGEX REPLACE "${L}ENGINE_PARAMS" "${engine_param},${L}ENGINE_PARAMS" tmp_string ${src_string})
  else ()
       string(APPEND tmp_string "${src_string}" " ${ARG}=${engine_param},${L}ENGINE_PARAMS")
  endif ()

  set(${dst_string} "${tmp_string}" PARENT_SCOPE)
endfunction()

function(MutateTestSet out_list name participant param in_list )
  set (tmp_list "")
  FOREACH (basename ${in_list})
    set (modname "${basename}.${name}")
    if ("${${basename}_CMD}" STREQUAL "") 
       message(SEND_ERROR "Staging-Common MutateTestSet ${basename} has no defined ${basename}_CMD")
    endif()
    add_engine_param(${participant} tmp_cmd "${${basename}_CMD}" ${param})
    set(${modname}_CMD ${tmp_cmd} PARENT_SCOPE)
    if (NOT "${${basename}_TIMEOUT}" STREQUAL "") 
        set(${modname}_TIMEOUT ${${basename}_TIMEOUT} PARENT_SCOPE)
    endif()
    if (NOT "${${basename}_PROPERTIES}" STREQUAL "") 
        set(${modname}_PROPERTIES ${${basename}_PROPERTIES} PARENT_SCOPE)
    endif()
    LIST(APPEND tmp_list ${modname})
  endforeach()
  if (STAGING_COMMON_TEST_SUPP_VERBOSE)
    message (STATUS "Setting ${out_list} to ${tmp_list}")
  endif()
  set(${out_list} ${tmp_list} PARENT_SCOPE)
endfunction()

function(add_common_test basename engine)
    set(testname "Staging.${basename}.${engine}")
    if ("${${basename}_CMD}" STREQUAL "") 
       message(SEND_ERROR "Staging-Common test ${basename} has no defined ${basename}_CMD")
    endif()
    string(FIND ${${basename}_CMD} "run_test.py" pos)
    if (NOT ${pos} EQUAL -1 )
        set(command "${PYTHON_EXECUTABLE} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${${basename}_CMD}")
        set(place 2)
    else()
        set(command "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${${basename}_CMD}")
	set(place 1)
    endif()
    remove_engine_params_placeholder(command  "${command}")
    separate_arguments(command)
    list(INSERT command ${place} "${engine}" "${testname}")
    if(NOT ADIOS2_RUN_MPI_MPMD_TESTS)
      list(APPEND command "--disable_mpmd")
    endif()
    add_test(NAME ${testname} COMMAND ${command})
    if(testname MATCHES "([1-9][0-9]*)x([1-9][0-9]*)")
      math(EXPR nprocs "${CMAKE_MATCH_1} + ${CMAKE_MATCH_2}")
      set_tests_properties(${testname} PROPERTIES PROCESSORS ${nprocs})
    endif()
    if (STAGING_COMMON_TEST_SUPP_VERBOSE)
	message ( STATUS "Adding test \"${testname}\" COMMAND \"${command}\"")
    endif()
    set (timeout "${${basename}_TIMEOUT}")
    if ("${timeout}" STREQUAL "")
      if (DEFINED MPIEXEC_EXECUTABLE AND "${MPIEXEC_EXECUTABLE}" MATCHES "srun")
        set (timeout "1000")
      else()
        set (timeout "60")
      endif()
    endif()

    set_tests_properties(${testname} PROPERTIES
        TIMEOUT ${timeout} ${${basename}_PROPERTIES}
    )
endfunction()

function(from_hex HEX DEC)
    string(TOUPPER "${HEX}" HEX)
    set(_res 0)
    string(LENGTH "${HEX}" _strlen)

    while(_strlen GREATER 0)
        math(EXPR _res "${_res} * 16")
        string(SUBSTRING "${HEX}" 0 1 NIBBLE)
        string(SUBSTRING "${HEX}" 1 -1 HEX)
        if(NIBBLE STREQUAL "A")
            math(EXPR _res "${_res} + 10")
        elseif(NIBBLE STREQUAL "B")
            math(EXPR _res "${_res} + 11")
        elseif(NIBBLE STREQUAL "C")
            math(EXPR _res "${_res} + 12")
        elseif(NIBBLE STREQUAL "D")
            math(EXPR _res "${_res} + 13")
        elseif(NIBBLE STREQUAL "E")
            math(EXPR _res "${_res} + 14")
        elseif(NIBBLE STREQUAL "F")
            math(EXPR _res "${_res} + 15")
        else()
            math(EXPR _res "${_res} + ${NIBBLE}")
        endif()

        string(LENGTH "${HEX}" _strlen)
    endwhile()
    
    set(${DEC} ${_res} PARENT_SCOPE)
endfunction()

function(import_bp_test BASENAME WRITE_SCALE READ_SCALE)
    set (WRITER_POSTFIX "Serial")
    set (READER_POSTFIX "Serial")
    if(ADIOS2_HAVE_MPI)
        set (WRITER_POSTFIX "MPI")
    endif()
    if(ADIOS2_HAVE_MPI)
        set (READER_POSTFIX "MPI")
    endif()
  set (${BASENAME}.${WRITE_SCALE}x${READ_SCALE}_CMD "run_test.py.$<CONFIG> -nw ${WRITE_SCALE} -nr ${READ_SCALE} --warg=-do_write --rarg=-do_read -w $<TARGET_FILE:Test.Engine.BP.${BASENAME}.${WRITER_POSTFIX}> -r $<TARGET_FILE:Test.Engine.BP.${BASENAME}.${READER_POSTFIX}>" PARENT_SCOPE)

endfunction()

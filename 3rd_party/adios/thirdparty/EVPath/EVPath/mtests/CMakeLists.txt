set(TESTS cmtest cmping cmconn bulktest take_test)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/cmtests.tsts "" )
foreach(T ${TESTS})
  add_executable(mtests_${T} ${T}.c)
  target_link_libraries(mtests_${T} evpath atl::atl)
  add_test(NAME mtests_${T} COMMAND mtests_${T})
  list(APPEND CMTESTS mtests_${T})
  file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmtests.tsts "mtests_${T}\n")
endforeach()

add_test(NAME mtests_cmtest_UDP COMMAND mtests_cmtest -t udp)
add_test(NAME mtests_cmping_UDP COMMAND mtests_cmping -t udp)
add_test(NAME mtests_non_blocking_bulk COMMAND mtests_bulktest -size 1024000)
set_tests_properties(mtests_non_blocking_bulk PROPERTIES
  ENVIRONMENT "CMNonBlockWrite=1")
add_executable(trans_test trans_test.c)
target_link_libraries(trans_test evpath atl::atl)

if(MPI_FOUND)
  add_executable(mpi_trans_test trans_test.c)
  target_compile_options(mpi_trans_test PRIVATE ${MPI_C_COMPILE_FLAGS})
  target_link_libraries(mpi_trans_test evpath atl::atl ${MPI_C_LIBRARIES})
  target_include_directories(mpi_trans_test PRIVATE ${MPI_C_INCLUDE_PATH})
  target_compile_defintion(mpi_trans_test PRIVATE BUILD_WITH_MPI=1)
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/transport_test.pl
  ${CMAKE_CURRENT_BINARY_DIR}/transport_test.pl
  COPYONLY)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/test_spec
  ${CMAKE_CURRENT_BINARY_DIR}/test_spec
  COPYONLY)
file(WRITE "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/correctness_spec" 
"%macro = (\n"
"        basic_size_list  => [ \"size:24 vectors:1 msg_count:10\",\n"
"			      \"size:102400 vectors:1 msg_count:3\" ],\n"
"        tiny_size_list  => [ \"size:24 vectors:1 msg_count:10 \"],\n"
"	all_read_buffer_styles  => [ \"take_receive_buffer:0\", \"take_receive_buffer:1\" ],\n"
"	all_write_buffer_styles  => [ \"reuse_write_buffers:0\", \"reuse_write_buffers:1\" ],\n"
"      );\n")
file (APPEND "${CMAKE_CURRENT_BINARY_DIR}/correctness_spec" 
"\$test_set{\"sockets\"} = \"basic_size_list all_read_buffer_styles all_write_buffer_styles\";\n")
if(RUN_NNTI_TESTS)
  file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/correctness_spec" 
    "\$test_set{\"nnti\"} = \"basic_size_list all_read_buffer_styles all_write_buffer_styles\";\n")
endif()
if(RUN_UDT4_TESTS)
  file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/correctness_spec" 
    "\$test_set{\"udt4\"} = \"basic_size_list all_read_buffer_styles all_write_buffer_styles\";\n")
endif()
if(RUN_ENET_TESTS)
  file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/correctness_spec" 
    "\$test_set{\"enet\"} = \"basic_size_list all_read_buffer_styles all_write_buffer_styles\";\n")
endif()
if(RUN_ZPL_ENET_TESTS)
  file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/correctness_spec" 
    "\$test_set{\"zplenet\"} = \"basic_size_list all_read_buffer_styles all_write_buffer_styles\";\n")
endif()
if(RUN_IB_TESTS)
  file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/correctness_spec" 
    "\$test_set{\"ib\"} = \"basic_size_list all_read_buffer_styles all_write_buffer_styles\";\n")
endif()

add_test(NAME all_transports
  COMMAND perl transport_test.pl -q -f ./correctness_spec)
add_test(NAME non_blocking_transports
  COMMAND perl transport_test.pl -q -f ./correctness_spec)
set_tests_properties(non_blocking_transports
  PROPERTIES ENVIRONMENT "CMNonBlockWrite=1")

if(TEST_INSTALL_DIRECTORY) 
  install(
    FILES
      ${CMAKE_CURRENT_BINARY_DIR}/cmtests.tsts
      ${CMAKE_CURRENT_SOURCE_DIR}/run_remote_tests.pl
    DESTINATION "${TEST_INSTALL_DIRECTORY}")
  install(TARGETS ${CMTESTS} DESTINATION "${TEST_INSTALL_DIRECTORY}")
endif() 

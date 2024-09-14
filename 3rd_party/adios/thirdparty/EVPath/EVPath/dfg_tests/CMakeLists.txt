if (HAVE_COD_H)
  set(TESTS evtest master_not_client multi_sink rawtest chain_test filter_test 
    tree_test multi_test dchain_test fail_chain_test reconfig_chain_test 
    anon_multi_test auto_tree_test router_test force_test self_reconfig_test 
    source_split named_stones delayed_submit metrics_test)
else()
  set(TESTS evtest master_not_client multi_sink rawtest 
    tree_test dchain_test fail_chain_test source_split)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

include_directories(${CMAKE_CURRENT_BINARY_DIR})
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/dfg_tests.tsts "" )
foreach(T ${TESTS})
  add_executable(dfg_${T} ${T}.c test_support.c)
  target_link_libraries(dfg_${T} evpath ffs::ffs)
  add_test(NAME dfg_${T} COMMAND dfg_${T})
  list(APPEND DFG_TESTS "dfg_${T}")
  file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/dfg_tests.tsts "dfg_${T}\n" )
endforeach()

if(UNIX)
  target_link_libraries(dfg_tree_test m)
  if (HAVE_COD_H)
    target_link_libraries(dfg_auto_tree_test m)
  endif()
endif()

if (HAVE_COD_H)
  add_test(NAME dfg_chain_test_one_proc  COMMAND dfg_chain_test 1)
endif()
add_test(NAME dfg_dchain_test_one_proc COMMAND dfg_dchain_test 1)

if(TEST_INSTALL_DIRECTORY) 
  install(TARGETS ${DFG_TESTS} DESTINATION "${TEST_INSTALL_DIRECTORY}")
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/dfg_tests.tsts
    DESTINATION "${TEST_INSTALL_DIRECTORY}")
endif() 

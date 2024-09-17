if (HAVE_COD_H)
  set (TESTS evtest bulktest http_test filter_test router_test router_test2 no_type_router_test filter2_test transform_test auto_test extract_test submit_test thin_test rawtest rawtest2 multi_thread multiq_test split_test executing_stone_test)
else()
  set (TESTS evtest bulktest http_test)
endif()
set (NO_INST_PROGS thin_client)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

include_directories(${CMAKE_CURRENT_BINARY_DIR})
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/evtests.tsts "" )
foreach(T ${TESTS})
  add_executable(${T} ${T}.c)
  target_link_libraries(${T} evpath atl::atl)
  add_test(NAME ${T} COMMAND ${T})
  file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/evtests.tsts "${T}\n" )
endforeach()
if (HAVE_COD_H)
  target_link_libraries(rawtest ffs::ffs)
  target_link_libraries(rawtest2 ffs::ffs)
endif()

add_test(NAME evtest_SELF COMMAND evtest)
set_property(TEST evtest_SELF PROPERTY
  ENVIRONMENT "CMSelfFormats=1"
)
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/evtests.tsts
  "CMSelfFormats=1 evtest\n" )

if (HAVE_SYS_EPOLL_H) 
   add_test(NAME evtest_SELECT COMMAND evtest -control select)
   add_test(NAME evtest_SELECT_ENV COMMAND evtest)
   set_property(TEST evtest_SELECT_ENV PROPERTY
      ENVIRONMENT "CMControlModule=select")
endif()

add_executable(thin_client thin_client.c)
target_link_libraries(thin_client evpath ffs::ffs)

add_test(NAME evtest_UDP COMMAND evtest -t udp)
if(ENET_FOUNMD)
  add_test(NAME evtest_ENET COMMAND evtest -t enet)
  add_test(NAME bulktest_ENET COMMAND bulktest -t enet)
endif()

if(BUILD_SHARED_LIBS)
  add_subdirectory(testdll)
endif()

if(TEST_INSTALL_DIRECTORY) 
  install(TARGETS ${TESTS} thin_client DESTINATION "${TEST_INSTALL_DIRECTORY}")
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/evtests.tsts
    DESTINATION "${TEST_INSTALL_DIRECTORY}")
endif() 

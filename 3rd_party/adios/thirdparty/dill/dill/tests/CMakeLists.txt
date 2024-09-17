
set(TESTS regress ctest call-test t1 pkg_test)

include_directories(BEFORE
  ${dill_BINARY_DIR}
  ${dill_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/regress.c
  COMMAND ${PERL_EXECUTABLE} -w
    ${CMAKE_CURRENT_SOURCE_DIR}/test-gen
    ${TEST_PERL_FLAGS} > regress.c
  MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/test-gen
)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/call-test.c
  COMMAND ${PERL_EXECUTABLE} -w
    ${CMAKE_CURRENT_SOURCE_DIR}/call-gen ${TEST_PERL_FLAGS} > call-test.c
  MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/call-gen 
)
foreach(TEST ${TESTS})
  add_executable(dill_${TEST} ${TEST}.c)
  target_link_libraries(dill_${TEST} dill ${DIS_LIBS})
  add_test(NAME dill_${TEST} COMMAND dill_${TEST} ${TEST_PERL_FLAGS})
endforeach()

if(NOT (CMAKE_C_COMPILER_ID MATCHES "MSVC" OR
     CMAKE_C_SIMULATE_ID MATCHES "MSVC"))
     target_link_libraries(dill_regress m)
     set_target_properties(dill_regress PROPERTIES COMPILE_FLAGS "-O0")
endif()

if(NOT (CMAKE_C_COMPILER_ID MATCHES "MSVC" OR
     CMAKE_C_SIMULATE_ID MATCHES "MSVC"))
    add_executable(cplus cplus.cc)
    target_link_libraries(cplus dill ${DIS_LIBS})
    add_test(NAME dill_cplus COMMAND cplus)
endif()

add_executable(stest stest.c)
target_link_libraries(stest dill ${DIS_LIBS})
add_test(NAME dill_stest COMMAND stest)

set_tests_properties (dill_call-test
  PROPERTIES PASS_REGULAR_EXPRESSION "No errors!")
set_tests_properties (dill_regress
  PROPERTIES PASS_REGULAR_EXPRESSION "No errors!")
set_tests_properties (dill_stest
  PROPERTIES PASS_REGULAR_EXPRESSION "hello world!
success!
hello world!
success!
hello world!
success!")
set_tests_properties (dill_ctest
  PROPERTIES PASS_REGULAR_EXPRESSION "Hello: 10 20 30 40
Hello: 10 20 30 40
Hello: 10 20 30 40 50 60 70 80 90 100
Hello: 1.000000e[+]01 2.000000e[+]01 3.000000e[+]01 4.000000e[+]01 5.000000e[+]01 6.000000e[+]01 7.000000e[+]01 8.000000e[+]01 9.000000e[+]01 1.000000e[+]02")
set_target_properties(dill_call-test PROPERTIES LINKER_LANGUAGE C)
set_target_properties(dill_ctest PROPERTIES LINKER_LANGUAGE C)
set_target_properties(dill_pkg_test PROPERTIES LINKER_LANGUAGE C)

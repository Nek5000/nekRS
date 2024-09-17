
set(TESTS basic_call general t1 opt branch prefix_test)

if(NATIVE_CG)
  list(APPEND TESTS pkg_test multi_test)
endif(NATIVE_CG)

include_directories(BEFORE
  ${dill_BINARY_DIR}
  ${dill_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/general.c
  COMMAND ${PERL_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/general.ops
  MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/general.ops
)

foreach(TEST ${TESTS})
  add_executable(v${TEST} ${TEST}.c)
  target_link_libraries(v${TEST} dill ${DIS_LIBS})
  add_test(NAME vtest_${TEST} COMMAND v${TEST})
endforeach()

if(NOT (CMAKE_C_COMPILER_ID MATCHES "MSVC" OR
     CMAKE_C_SIMULATE_ID MATCHES "MSVC"))
     set_target_properties(vgeneral PROPERTIES LINKER_LANGUAGE C COMPILE_FLAGS "-O0")
ENDIF()

set_tests_properties(vtest_basic_call
  PROPERTIES PASS_REGULAR_EXPRESSION "########## A
In ff  a=1, b=2, c=3, d=4, e=5, f=6, g=7, h=8, i=9, j=0
In ff  a=1, b=2, c=3, d=4, e=5, f=6, g=7, h=8, i=9, j=0
[*][*]135=135
########## B
In gg  a=1, b=2
[*][*]3=3
########## C
In gg  a=1, b=2
[*][*]3=3
########## D
In gg  a=1, b=2
[*][*]3=3
########## E
In gg  a=1, b=2
In gg  a=3, b=4
In gg  a=5, b=6
In gg  a=7, b=8
[*][*]48=48
########## F
expect: values are 5, 3.14159, hello!
values are 5, 3.14159, hello!
########## end
")


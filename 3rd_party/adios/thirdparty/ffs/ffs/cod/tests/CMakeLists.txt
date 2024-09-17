set (TESTS t1 t2 t3 t5 t6 t8 t9 t10 t12 compound_assignment gray general atl_test strings control structs time)
set (PROGS t7)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

ADD_CUSTOM_COMMAND(
   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/general.c
   COMMAND perl ${CMAKE_CURRENT_SOURCE_DIR}/general.ops
   MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/general.ops
)

foreach (TEST ${TESTS} )
    ADD_EXECUTABLE(${TEST} ${TEST}.c data_funcs.c)
    TARGET_LINK_LIBRARIES(${TEST} ffs )
    TARGET_LINK_LIBRARIES(${TEST} ${ATL_LIBRARIES} )
    ADD_TEST(NAME ${TEST} COMMAND ${TEST})
endforeach()

if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
   SET_SOURCE_FILES_PROPERTIES(general.c PROPERTIES COMPILE_FLAGS -O0)
endif()

set_tests_properties (t2
  PROPERTIES PASS_REGULAR_EXPRESSION "Expect -> .values are is 5, 3.14159, hello!.
values are is 5, 3.14159, hello!")

ADD_EXECUTABLE(t4 t4.c data_funcs.c)
TARGET_LINK_LIBRARIES(t4 ffs )
TARGET_LINK_LIBRARIES(t4 ${ATL_LIBRARIES} )
ADD_TEST(NAME t4 COMMAND t4 -output t4.actual_out)

add_test(NAME compare_t4
  COMMAND ${CMAKE_COMMAND}
    -E compare_files t4.actual_out ${CMAKE_CURRENT_SOURCE_DIR}/t4.out)

foreach (PROG ${PROGS} )
    ADD_EXECUTABLE(${PROG} ${PROG}.c)
    TARGET_LINK_LIBRARIES(${PROG} ffs )
    TARGET_LINK_LIBRARIES(${PROG} ${ATL_LIBRARIES} )
endforeach()


set (TESTS align_test compat_test)
set (PROGS format_test self_format_test scale_test)

foreach (TEST ${TESTS} )
    ADD_EXECUTABLE(${TEST} ${TEST}.c test_funcs.c)
    TARGET_LINK_LIBRARIES(${TEST} ffs )
    ADD_TEST(NAME ${TEST} COMMAND ${TEST})
endforeach()

foreach (PROG ${PROGS} )
    ADD_EXECUTABLE(${PROG} ${PROG}.c test_funcs.c)
    TARGET_LINK_LIBRARIES(${PROG} ffs )
endforeach()

install(TARGETS ${PROGS}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT ${FFS_RUNTIME_COMPONENT}
)

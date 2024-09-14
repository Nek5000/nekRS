set (TESTS)
if (NOT NO_SOCKETS)
    set (PROGS FFSdump FFScp FFSsort)
else ()
    set (PROGS)
endif()

foreach (TEST ${TESTS} )
    ADD_EXECUTABLE(${TEST} ${TEST}.c)
    TARGET_LINK_LIBRARIES(${TEST} ffs )
    ADD_TEST(${TEST} "${TEST}")
endforeach()

foreach (PROG ${PROGS} )
    ADD_EXECUTABLE(${PROG} ${PROG}.c)
    TARGET_LINK_LIBRARIES(${PROG} ffs )
endforeach()

install(TARGETS ${PROGS}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT ${FFS_RUNTIME_COMPONENT}
)

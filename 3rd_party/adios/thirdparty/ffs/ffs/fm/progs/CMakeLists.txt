FIND_PACKAGE (Threads)

set (TESTS)
if (NOT NO_SOCKETS)
    set (PROGS format_dump format_cmd format_info)
    ADD_EXECUTABLE(format_server format_server.c server.c)
    TARGET_LINK_LIBRARIES(format_server ffs ${CMAKE_THREAD_LIBS_INIT})
    set (FORMAT_SERVER format_server)
else ()
    set (PROGS)
    set (FORMAT_SERVER)
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

if (NOT NO_SOCKETS)
   install(TARGETS ${PROGS} ${FORMAT_SERVER}
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT ${FFS_RUNTIME_COMPONENT}
     )
endif()
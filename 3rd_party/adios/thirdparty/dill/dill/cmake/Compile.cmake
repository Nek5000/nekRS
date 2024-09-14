
FUNCTION(COMPILE RESULT SOURCE)
    MESSAGE(STATUS "Compiling ${SOURCE}")
    # Ensure SOURCE is absolute:
    IF(NOT IS_ABSOLUTE ${SOURCE})
        SET(SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE})
    ENDIF()
    # Set up CMakeLists.txt for static library:
    FILE(WRITE
        ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/compile/CMakeLists.txt
        "ADD_LIBRARY(compile STATIC ${SOURCE})"
    )
    # Configure:
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} . WORKING_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/compile
        OUTPUT_VARIABLE LOG1 ERROR_VARIABLE LOG1
    )
    # Build:
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} --build
        ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/compile
        RESULT_VARIABLE RESVAR OUTPUT_VARIABLE LOG2 ERROR_VARIABLE LOG2
    )
    # Clean up:
    FILE(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/compile)
    # Set up log:
    IF(ARGC GREATER 2)
        SET(${ARGV2} "${LOG1}${LOG2}" PARENT_SCOPE)
    ENDIF()
    # Set up result:
    IF(RESVAR EQUAL 0)
        SET(${RESULT} TRUE PARENT_SCOPE)
    ELSE()
        SET(${RESULT} FALSE PARENT_SCOPE)
    ENDIF()
ENDFUNCTION()

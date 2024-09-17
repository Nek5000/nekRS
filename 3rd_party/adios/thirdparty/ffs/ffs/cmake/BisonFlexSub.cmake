FUNCTION (SETUP_BISON_FLEX_SUB)

IF ((${CMAKE_SYSTEM_NAME} STREQUAL "Darwin") OR
   (${CMAKE_SYSTEM_NAME} STREQUAL "Linux") OR
   (${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD"))
   set (BISON_FLEX_PRECOMPILE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cod/pregen_source/Linux")
elseif (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
   set (BISON_FLEX_PRECOMPILE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cod/pregen_source/Windows")
else()
   if (NOT BISON_FOUND)
      message (FATAL_ERROR "Bison was not found and no pregenerated Bison/Flex"
              "source is available for ${CMAKE_SYSTEM_NAME}. Please install Bison or Yacc")

   else()
      message (FATAL_ERROR "Flex was not found and no pregenerated Bison/Flex" 
      	      "source is available for ${CMAKE_SYSTEM_NAME}. Please install Bison or Yacc")
   endif()
ENDIF()

ADD_CUSTOM_COMMAND(OUTPUT cod.tab.c
        COMMAND ${CMAKE_COMMAND} -E copy ${BISON_FLEX_PRECOMPILE_DIR}/cod.tab.c ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND ${CMAKE_COMMAND} -E copy ${BISON_FLEX_PRECOMPILE_DIR}/cod.tab.h ${CMAKE_CURRENT_BINARY_DIR}
	COMMENT "Using pre-generated Bison Output from ${BISON_FLEX_PRECOMPILE_DIR}")
ADD_CUSTOM_COMMAND(OUTPUT lex.yy.c
	COMMAND ${CMAKE_COMMAND} -E copy ${BISON_FLEX_PRECOMPILE_DIR}/lex.yy.c ${CMAKE_CURRENT_BINARY_DIR}
	COMMENT "Using pre-generated Flex Output from ${BISON_FLEX_PRECOMPILE_DIR}")

set (BISON_CODParser_OUTPUT_SOURCE cod.tab.c PARENT_SCOPE)
ENDFUNCTION()

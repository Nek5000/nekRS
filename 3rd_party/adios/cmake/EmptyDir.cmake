if(NOT DIR_NAME)
  message(FATAL_ERROR "DIR_NAME is not defined")
endif()

message("DEBUG: DIR_NAME |${DIR_NAME}|")
if(EXISTS "${DIR_NAME}")
  file(REMOVE_RECURSE "${DIR_NAME}")
endif()
file(MAKE_DIRECTORY "${DIR_NAME}")

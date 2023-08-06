if(ENABLE_SENSEI)
  find_package(SENSEI REQUIRED)

  MESSAGE(STATUS ${SENSEI_LIBRARIES})
  if(NOT SENSEI_FOUND)
    MESSAGE(FATAL_ERROR "Could not find SENSEI.")
  endif()
 
  find_package(Python3 COMPONENTS Interpreter Development)
 
  if(NOT Python3_FOUND)
    MESSAGE(FATAL_ERROR "Could not find Python.")
  endif()

  message("Python3_FOUND:${Python3_FOUND}")
  message("Python3_VERSION:${Python3_VERSION}")
  message("Python3_Development_FOUND:${Python3_Development_FOUND}")
  message("Python3_LIBRARIES:${PYTHON3_LIBRARIES}")

  set(SENSEI_LIBS sensei ${Python3_LIBRARIES} Python3::Python)
  set(SRC src/sensei/DataAdaptor.cxx src/sensei/Bridge.cxx)

  add_library(neksensei STATIC ${SRC})
  target_link_libraries(neksensei ${SENSEI_LIBS})
  add_definitions(-DENABLE_SENSEI)

  target_link_directories(neksensei 
    PRIVATE
    ${SENSEI_DIR}/../..)

endif()


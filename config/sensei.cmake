if(ENABLE_SENSEI)
  find_package(SENSEI REQUIRED)

  MESSAGE(STATUS ${SENSEI_LIBRARIES})
  if(NOT SENSEI_FOUND)
    MESSAGE(FATAL_ERROR "Could not find SENSEI.")
  endif()
 
  find_package(Python3 COMPONENTS Interpreter Development)
 
  if(NOT Python_FOUND)
    MESSAGE(FATAL_ERROR "Could not find Python.")
  endif()

  message("Python_FOUND:${Python_FOUND}")
  message("Python_VERSION:${Python_VERSION}")
  message("Python_Development_FOUND:${Python_Development_FOUND}")
  message("Python_LIBRARIES:${PYTHON_LIBRARIES}")

  set(SENSEI_LIBS sensei ${Python_LIBRARIES} Python3::Python)
  set(SRC src/sensei/DataAdaptor.cxx src/sensei/Bridge.cxx)

  add_library(neksensei STATIC ${SRC})
  target_link_libraries(neksensei ${SENSEI_LIBS})
  add_definitions(-DENABLE_SENSEI)

  target_link_directories(neksensei 
    PRIVATE
    ${SENSEI_DIR}/../..)

endif()


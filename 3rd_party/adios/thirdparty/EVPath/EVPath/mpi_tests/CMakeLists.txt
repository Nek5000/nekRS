set(TESTS simple_cm bulk_cm )

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

foreach(T ${TESTS})
  add_executable(${T} ${T}.c)
  target_link_libraries(${T} evpath atl::atl ${MPI_C_LIBRARIES})
  target_include_directories(${T} PRIVATE ${MPI_C_INCLUDE_PATH})
  target_compile_options(${T} PRIVATE ${MPI_C_COMPILE_FLAGS})
endforeach()

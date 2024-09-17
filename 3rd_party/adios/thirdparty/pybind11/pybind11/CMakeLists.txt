add_library(pybind11 INTERFACE)
target_include_directories(pybind11 INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
install(TARGETS pybind11 EXPORT adios2Exports)

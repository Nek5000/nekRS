Python_add_library(adios2_py MODULE
  WITH_SOABI
  py11ADIOS.cpp
  py11IO.cpp
  py11Variable.cpp
  py11Attribute.cpp
  py11Engine.cpp
  py11Operator.cpp
  py11Query.cpp
  py11glue.cpp
)
target_compile_definitions(adios2_py PRIVATE "ADIOS2_PYTHON_MODULE_NAME=adios2_bindings${ADIOS2_LIBRARY_SUFFIX}")
if(ADIOS2_HAVE_MPI)
  target_sources(adios2_py PRIVATE
    py11ADIOSMPI.cpp
    py11IOMPI.cpp
  )
  set(maybe_adios2_cxx11_mpi adios2_cxx11_mpi)
  set(maybe_adios2_core_mpi adios2_core_mpi)
  set(maybe_mpi4py Python::mpi4py)
else()
  set(maybe_adios2_cxx11_mpi)
  set(maybe_adios2_core_mpi)
  set(maybe_mpi4py)
endif()
target_link_libraries(adios2_py PRIVATE
  ${maybe_adios2_cxx11_mpi} adios2_cxx11
  ${maybe_adios2_core_mpi} adios2_core
  adios2::thirdparty::pybind11
  ${maybe_mpi4py}
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py.in
  ${CMAKE_PYTHON_OUTPUT_DIRECTORY}/adios2/bindings/__init__.py
  @ONLY
)

set_target_properties(adios2_py PROPERTIES
  CXX_VISIBILITY_PRESET hidden
  OUTPUT_NAME adios2_bindings${ADIOS2_LIBRARY_SUFFIX}
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_PYTHON_OUTPUT_DIRECTORY}/adios2/bindings
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_PYTHON_OUTPUT_DIRECTORY}/adios2/bindings
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_PYTHON_OUTPUT_DIRECTORY}/adios2/bindings
  PDB_OUTPUT_DIRECTORY ${CMAKE_PYTHON_OUTPUT_DIRECTORY}/adios2/bindings
  COMPILE_PDB_OUTPUT_DIRECTORY ${CMAKE_PYTHON_OUTPUT_DIRECTORY}/adios2/bindings
)

set(install_location "${CMAKE_INSTALL_PYTHONDIR}/adios2")
if (ADIOS2_USE_PIP)
  set(install_location ${CMAKE_INSTALL_LIBDIR})
endif()

string(REGEX REPLACE "[^/]+" ".." relative_base "${install_location}/bindings")
if(CMAKE_SYSTEM_NAME MATCHES "Linux")
  set_property(TARGET adios2_py APPEND PROPERTY
    INSTALL_RPATH "$ORIGIN/${relative_base}/${CMAKE_INSTALL_LIBDIR}"
  )
endif()

install(TARGETS adios2_py
  DESTINATION ${install_location}/bindings
  COMPONENT adios2_python-python
)
install(FILES ${CMAKE_PYTHON_OUTPUT_DIRECTORY}/adios2/bindings/__init__.py
  DESTINATION ${install_location}/bindings
  COMPONENT adios2_python-python
)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test
  DESTINATION ${install_location}
  COMPONENT adios2_python-python
)

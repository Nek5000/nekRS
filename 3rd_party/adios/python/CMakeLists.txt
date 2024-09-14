#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

add_custom_target(python_api ALL COMMAND ${CMAKE_COMMAND} -E copy_directory
                  ${CMAKE_CURRENT_SOURCE_DIR}/adios2
                  ${CMAKE_PYTHON_OUTPUT_DIRECTORY}/adios2)

set(install_location ${CMAKE_INSTALL_PYTHONDIR}/adios2)
if (ADIOS2_USE_PIP)
  set(install_location ${CMAKE_INSTALL_LIBDIR})
endif()

install(DIRECTORY ${CMAKE_PYTHON_OUTPUT_DIRECTORY}/adios2/
  DESTINATION ${install_location}
  COMPONENT adios2_python-python
)

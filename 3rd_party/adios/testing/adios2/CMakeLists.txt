#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

add_subdirectory(interface)
add_subdirectory(unit)
add_subdirectory(engine)
add_subdirectory(transports)
add_subdirectory(bindings)
add_subdirectory(xml)
add_subdirectory(yaml)
add_subdirectory(performance)
add_subdirectory(helper)
add_subdirectory(hierarchy)
add_subdirectory(backward_compatibility)
if (ADIOS2_HAVE_Derived_Variable)
add_subdirectory(derived)
endif()

if (ADIOS2_HAVE_Python)
add_subdirectory(python)
endif()

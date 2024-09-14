#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

add_subdirectory(basics)
add_subdirectory(hello)
add_subdirectory(plugins)
add_subdirectory(simulations)
add_subdirectory(useCases)
if(ADIOS2_USE_Campaign)
    add_subdirectory(campaign)
endif()

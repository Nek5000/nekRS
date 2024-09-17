#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

function(GenerateSSTHeaderConfig)
  set(SST_CONFIG_DEFINES)
  foreach(OPT IN LISTS ARGN)
    string(TOUPPER ${OPT} OPT_UPPER)
    string(APPEND SST_CONFIG_DEFINES "
/* CMake Option: SST_USE_${OPT}=OFF */
#cmakedefine SST_HAVE_${OPT_UPPER}
")
    if(ADIOS2_SST_HAVE_${OPT})
      set(SST_HAVE_${OPT_UPPER} 1)
    else()
      set(SST_HAVE_${OPT_UPPER})
    endif()
  endforeach()

  configure_file(
    ${ADIOS2_SOURCE_DIR}/source/adios2/toolkit/sst/SSTConfig.h.in
    ${ADIOS2_BINARY_DIR}/source/adios2/toolkit/sst/SSTConfig.h.in
  )
  configure_file(
    ${ADIOS2_BINARY_DIR}/source/adios2/toolkit/sst/SSTConfig.h.in
    ${ADIOS2_BINARY_DIR}/source/adios2/toolkit/sst/SSTConfig.h
  )
endfunction()

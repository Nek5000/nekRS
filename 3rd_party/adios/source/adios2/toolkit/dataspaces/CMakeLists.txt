configure_file(
    ${ADIOS2_SOURCE_DIR}/source/adios2/toolkit/dataspaces/DSpacesConfig.h.in
    ${ADIOS2_BINARY_DIR}/source/adios2/toolkit/dataspaces/DSpacesConfig.h.in
  )
  configure_file(
    ${ADIOS2_BINARY_DIR}/source/adios2/toolkit/dataspaces/DSpacesConfig.h.in
    ${ADIOS2_BINARY_DIR}/source/adios2/toolkit/dataspaces/DSpacesConfig.h
  )

set(KERNEL_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/src")

# just look for okl file to find directories with kernels
file(GLOB_RECURSE KERNEL_FILES RELATIVE ${KERNEL_SOURCE} ${KERNEL_SOURCE}/*.okl)

set (KERNEL_DIRS "")
foreach(F IN LISTS KERNEL_FILES)
  get_filename_component(KERNEL_DIR ${F} DIRECTORY)
  list(APPEND KERNEL_DIRS ${KERNEL_DIR})
endforeach()

foreach(DIR IN LISTS KERNEL_DIRS)
  file(GLOB KERNEL_DIR_FILES  RELATIVE ${KERNEL_SOURCE} ${KERNEL_SOURCE}/${DIR}/*)
  foreach(F IN LISTS KERNEL_DIR_FILES)
    STRING(REGEX REPLACE "/kernels/" "/" FNEW ${F})
    set(DST_IN "${CMAKE_INSTALL_PREFIX}/kernels/${FNEW}")
    cmake_path(REMOVE_FILENAME DST_IN OUTPUT_VARIABLE DST)
    INSTALL(FILES ${KERNEL_SOURCE}/${F} DESTINATION ${DST})
  endforeach()
endforeach()
unset(KERNEL_DIRS)
unset(DST_IN)
unset(KERNEL_SOURCE)

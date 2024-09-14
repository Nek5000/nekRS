set(_NVML_EXTRA_INC_ARGS PATHS /usr/include/nvidia/gdk)
if(GPU_DEPLOYMENT_KIT_ROOT_DIR)
  set(_NVML_EXTRA_LIB_ARGS HINT
    "${GPU_DEPLOYMENT_KIT_ROOT_DIR}/src/gdk/nvml/lib")
  list(APPEND _NVML_EXTRA_INC_ARGS
    HINTS "${GPU_DEPLOYMENT_KIT_ROOT_DIR}/include/nvidia/gdk")
endif()
  
find_library(NVML_LIBRARY nvidia-ml ${_NVML_EXTRA_LIB_ARGS})
find_path(NVML_INCLUDE_DIR nvml.h ${_NVML_EXTRA_INC_ARGS})
mark_as_advanced(NVML_LIBRARY NVML_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NVML DEFAULT_MSG
  NVML_LIBRARY NVML_INCLUDE_DIR)

if(NVML_FOUND)
  set(NVML_LIBRARIES ${NVML_LIBRARY})
  set(NVML_INCLUDE_DIRS ${NVML_INCLUDE_DIR})
  if(NOT TARGET nvml::nvml)
    add_library(nvml::nvml UNKNOWN IMPORTED)
    set_target_properties(nvml::nvml PROPERTIES
      IMPORTED_LOCATION "${NVML_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${NVML_INCLUDE_DIR}")
  endif()
endif()


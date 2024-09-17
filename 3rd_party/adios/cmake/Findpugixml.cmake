include(FindPackageHandleStandardArgs)

function(__pugixml_get_version)
  find_package(PkgConfig)
  if(PkgConfig_FOUND)
    pkg_check_modules(_pugixml_pc pugixml)
    if(_pugixml_pc_FOUND)
      set(pugixml_VERSION "${_pugixml_pc_VERSION}" PARENT_SCOPE)
    endif()
  endif()
endfunction()

# We need to disable version checking, since pugixml does not provide it.
set(_pugixml_var_suffixes VERSION VERSION_MAJOR VERSION_MINOR VERSION_PATCH VERSION_TWEAK VERSION_COUNT REQUIRED)
foreach(_suffix IN LISTS _pugixml_var_suffixes)
  set(_pugixml_save_FIND_${_suffix} ${pugixml_FIND_${_suffix}})
  unset(pugixml_FIND_${_suffix})
endforeach()
find_package(pugixml CONFIG)
foreach(_suffix IN LISTS _pugixml_var_suffixes)
  set(pugixml_FIND_${_suffix} ${_pugixml_save_FIND_${_suffix}})
endforeach()

if(pugixml_FOUND AND NOT DEFINED pugixml_VERSION)
  __pugixml_get_version()
endif()

find_package_handle_standard_args(pugixml CONFIG_MODE)

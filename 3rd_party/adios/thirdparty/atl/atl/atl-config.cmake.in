include(FindPackageHandleStandardArgs)
set(${CMAKE_FIND_PACKAGE_NAME}_CONFIG ${CMAKE_CURRENT_LIST_FILE})
find_package_handle_standard_args(${CMAKE_FIND_PACKAGE_NAME} CONFIG_MODE)

if(NOT TARGET atl::atl)
  include("${CMAKE_CURRENT_LIST_DIR}/atl-targets.cmake")
endif()

set(ATL_LIBRARIES atl::atl)
set(ATL_INCLUDE_DIRS
  $<TARGET_PROPERTY:atl::atl,INTERFACE_INCLUDE_DIRECTORIES>
)

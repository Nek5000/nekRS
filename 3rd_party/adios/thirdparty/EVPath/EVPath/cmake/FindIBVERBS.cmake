######################################################
# Once done this will define
#  IBVERBS_FOUND - System has ibverbs
#  IBVERBS_LIBRARIES - The libraries needed to use ibverbs

######################################################
include(CheckFunctionExists)
check_function_exists(ibv_create_qp HAVE_IBV)
if(HAVE_IBV)
  set(IBVERBS_LIBRARY TRUE)
else()
  check_library_exists(ibverbs ibv_create_qp "" HAVE_IBVERBS_IBV)
  if(HAVE_IBVERBS_IBV)
    set(IBVERBS_LIBRARY ibverbs)
  endif()
endif()

######################################################
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(IBVERBS DEFAULT_MSG
  IBVERBS_LIBRARY)
if(IBVERBS_FOUND)
  if(IBVERBS_LIBRARY STREQUAL "TRUE")
    unset(IBVERBS_LIBRARY)
  endif()
  set(IBVERBS_LIBRARIES ${IBVERBS_LIBRARY})
endif()

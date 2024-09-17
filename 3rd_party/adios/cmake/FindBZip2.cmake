#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

# This module is already included in new versions of CMake
if(CMAKE_VERSION VERSION_LESS 3.7)
  include(${CMAKE_CURRENT_LIST_DIR}/upstream/FindBZip2.cmake)
else()
  include(${CMAKE_ROOT}/Modules/FindBZip2.cmake)
endif()

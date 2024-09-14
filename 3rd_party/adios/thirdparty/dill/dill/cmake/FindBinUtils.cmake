#
#  BinUtils_FOUND - system has the BinUtils library
#  BinUtils_INCLUDE_DIRS - the BinUtils include directory
#  BinUtils_LIBRARIES - The libraries needed to use BinUtils
#
#  Targets:
#    binutils::opcodes
#    binutils::bfd
#

find_path(BinUtils_INCLUDE_DIR dis-asm.h)
find_library(BinUtils_opcodes_LIBRARY opcodes)
find_library(BinUtils_bfd_LIBRARY bfd)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BinUtils DEFAULT_MSG
  BinUtils_opcodes_LIBRARY BinUtils_bfd_LIBRARY BinUtils_INCLUDE_DIR
)
if(BinUtils_FOUND)
  set(BinUtils_INCLUDE_DIRS ${BinUtils_INCLUDE_DIR})
  set(BinUtils_LIBRARIES ${BinUtils_opcodes_LIBRARY} ${BinUtils_bfd_LIBRARY})

  if(NOT TARGET binutils::bfd)
    add_library(binutils::bfd UNKNOWN IMPORTED)
    set_target_properties(binutils::bfd PROPERTIES
      IMPORTED_LOCATION              ${BinUtils_bfd_LIBRARY}
      INTERFACE_INCLUDE_DIRECTORIES  ${BinUtils_INCLUDE_DIR}
    )
  endif()

  if(NOT TARGET binutils::opcodes)
    add_library(binutils::opcodes UNKNOWN IMPORTED)
    set_target_properties(binutils::opcodes PROPERTIES
      IMPORTED_LOCATION              ${BinUtils_opcodes_LIBRARY}
      INTERFACE_INCLUDE_DIRECTORIES  ${BinUtils_INCLUDE_DIR}
      INTERFACE_LINK_LIBRARIES       binutils::bfd
    )
  endif()
endif()


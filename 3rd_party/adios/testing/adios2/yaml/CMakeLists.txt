#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

gtest_add_tests_helper(YAMLConfig MPI_ALLOW "" "" "")

foreach(tgt IN LISTS Test.YAMLConfig-TARGETS)
  target_compile_definitions(${tgt} PRIVATE
    "YAML_CONFIG_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
  )
endforeach()

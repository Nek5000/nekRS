cmake_minimum_required(VERSION 3.12)

# Write a file containing the ADIOS2 version to use in pip packaging.  This is
# how we work around the issue that scikit-build-core does not yet have a way
# to get the package version from CMake.
#
# There are several issues on scikit-build-core similar to these:
#
#   - https://github.com/scikit-build/scikit-build-core/issues/172
#   - https://github.com/scikit-build/scikit-build-core/issues/176
#
# Once there is a general solution for that, we can dispense with
# invoking this to write a custom version file to be consumed both
# by ADIOSFunctions.cmake and pyproject.toml.

set(adios2_pip_package_version "UNKNOWN")

if (DEFINED ENV{ADIOS2_CUSTOM_VERSION_OVERRIDE} AND NOT "$ENV{ADIOS2_CUSTOM_VERSION_OVERRIDE}" STREQUAL "")
  # Use the override version from env var
  set(adios2_pip_package_version $ENV{ADIOS2_CUSTOM_VERSION_OVERRIDE})
else()
  # Compute the version using git describe
  if(NOT GIT_COMMAND)
    find_program(GIT_COMMAND git)
  endif()
  if(GIT_COMMAND)
    execute_process(
      COMMAND git describe
      RESULT_VARIABLE res
      OUTPUT_VARIABLE out
      ERROR_QUIET
    )
    if(res EQUAL 0)
      if (out MATCHES "^v([^-]*)$")
        set(ver_tag ${CMAKE_MATCH_1})
        set(adios2_pip_package_version "${ver_tag}")
      elseif (out MATCHES "^v([^-]*)-([^-]*)-g[a-f0-9]*")
        set(ver_tag ${CMAKE_MATCH_1})
        set(ver_ncommits ${CMAKE_MATCH_2})
        math(EXPR ver_tweak "100000+${ver_ncommits}")
        set(adios2_pip_package_version "${ver_tag}.${ver_tweak}")
      elseif (out MATCHES "^v([^-]*)-rc([0-9]*)(-([^-]*)-g[a-f0-9]*)?")
        set(ver_tag ${CMAKE_MATCH_1})
        set(ver_rc ${CMAKE_MATCH_2})
        set(ver_ncommits 0)
        if (CMAKE_MATCH_COUNT EQUAL 4)
          set(ver_ncommits ${CMAKE_MATCH_4})
        endif()
        math(EXPR ver_tweak "(1000*${ver_rc})+${ver_ncommits}")
        set(adios2_pip_package_version "${ver_tag}.${ver_tweak}")
      endif()
    endif()
    if (NOT res EQUAL 0 OR out MATCHES "fatal: not a git repository")
      message(FATAL_ERROR "Must be in a git repository to run this script")
    endif()
  endif()
endif()

message(${adios2_pip_package_version})
file(WRITE "VERSION.TXT" ${adios2_pip_package_version})

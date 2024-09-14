#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#
#
# ADIOS Common Dashboard Script
#
# This script contains basic dashboard driver code common to all
# clients.
#
#   # Client maintainer: me@mydomain.net
#   set(CTEST_SITE "machine.site")
#   set(CTEST_BUILD_NAME "Platform-Compiler")
#   set(CTEST_CONFIGURATION_TYPE Debug)
#   set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
#   include(${CTEST_SCRIPT_DIRECTORY}/adios_common.cmake)
#
# Then run a scheduled task (cron job) with a command line such as
#
#   ctest -S ~/Dashboards/Scripts/my_dashboard.cmake -V
#
# By default the source and build trees will be placed in the path
# "../My Tests/" relative to your script location.
#
# The following variables may be set before including this script
# to configure it:
#
#   dashboard_model           = Nightly | Experimental
#   dashboard_root_name       = Change name of "MyTests" directory
#   dashboard_source_name     = Name of source directory (adios)
#   dashboard_binary_name     = Name of binary directory (adios-build)
#   dashboard_cache           = Initial CMakeCache.txt file content

#   dashboard_do_checkout  = True to enable source checkout via git
#   dashboard_do_update    = True to enable source update
#   dashboard_do_configure = True to enable the Configure step
#   dashboard_do_build     = True to enable the Build step
#   dashboard_do_test      = True to enable the Test step
#   dashboard_do_coverage  = True to enable coverage (ex: gcov)
#   dashboard_do_memcheck  = True to enable memcheck (ex: valgrind)

#   CTEST_GIT_COMMAND     = path to git command-line client
#   CTEST_BUILD_FLAGS     = build tool arguments (ex: -j2)
#   CTEST_DASHBOARD_ROOT  = Where to put source and build trees
#   CTEST_TEST_CTEST      = Whether to run long CTestTest* tests
#   CTEST_TEST_TIMEOUT    = Per-test timeout length
#   CTEST_TEST_ARGS       = ctest_test args (ex: PARALLEL_LEVEL 4)
#   CMAKE_MAKE_PROGRAM    = Path to "make" tool to use
#
# Options to configure Git:
#   dashboard_git_url      = Custom git clone url
#   dashboard_git_branch   = Custom remote branch to track
#   dashboard_git_crlf     = Value of core.autocrlf for repository
#

# For Makefile generators the script may be executed from an
# environment already configured to use the desired compilers.
# Alternatively the environment may be set at the top of the script:
#
#   set(ENV{CC}  /path/to/cc)   # C compiler
#   set(ENV{CXX} /path/to/cxx)  # C++ compiler
#   set(ENV{FC}  /path/to/fc)   # Fortran compiler (optional)
#   set(ENV{LD_LIBRARY_PATH} /path/to/vendor/lib) # (if necessary)

set(CTEST_PROJECT_NAME "ADIOS2")
set(CTEST_DROP_SITE "open.cdash.org")
if(NOT dashboard_git_url)
  set(dashboard_git_url "https://github.com/ornladios/ADIOS2.git")
endif()

if(NOT DEFINED CTEST_TEST_TIMEOUT)
  set(CTEST_TEST_TIMEOUT 120)
endif()
if(NOT DEFINED CTEST_SUBMIT_NOTES) 
  set(CTEST_SUBMIT_NOTES TRUE)
endif()
if(NOT dashboard_root_name)
  set(dashboard_root_name "Builds/My Tests")
endif()
if(NOT dashboard_source_name)
  set(dashboard_source_name "ADIOS2")
endif()
if(NOT dashboard_model)
  set(dashboard_model Experimental)
endif()

if(NOT DEFINED ADIOS_TEST_REPEAT)
  set(ADIOS_TEST_REPEAT 5)
endif()

if((CMAKE_VERSION VERSION_GREATER 3.16.20191201 ) AND
   (ADIOS_TEST_REPEAT GREATER 0) AND
   NOT "REPEAT" IN_LIST CTEST_TEST_ARGS)
  list(APPEND CTEST_TEST_ARGS REPEAT "UNTIL_PASS:${ADIOS_TEST_REPEAT}")
endif()

list(APPEND CTEST_UPDATE_NOTES_FILES "${CMAKE_CURRENT_LIST_FILE}")
include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)

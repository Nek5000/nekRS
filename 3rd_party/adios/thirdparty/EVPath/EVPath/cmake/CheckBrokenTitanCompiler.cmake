#
#  CHECK_BROKEN_TITAN_COMPILER()
#
#
#

FUNCTION (CHECK_BROKEN_TITAN_COMPILER )
  if (NOT ${CMAKE_CXX_COMPILER_ID} STREQUAL "PGI")
     RETURN()
  endif()
  if (NOT EXISTS "${CMAKE_PLATFORM_INFO_DIR}/CMakeCCompiler.cmake")
     RETURN()
  endif()
  FILE(READ "${CMAKE_PLATFORM_INFO_DIR}/CMakeCCompiler.cmake" contents)

  # Convert file contents into a CMake list (where each element in the list
  # is one line of the file)
  #

  if ("${contents}" MATCHES ".*stdc[+][+].*") 
    message (FATAL_ERROR 
"Hello.  I'm afraid I have bad news.  You seem to be configuring this "
"package on a machine with a broken 'cc' wrapper script.  Most likely "
"you are on ORNL's Titan or a similar Cray machine which is improperly "
"adding '-lstdc++' to the link line whenever the module 'cray-libsci' "
"is loaded.  stdc++ is a GCC library used by some non-GNU compilers, "
"but not PGI.  At this point, cmake has done a 'cc -v', snagged what "
"it assumes to be a valid implicit dependency and squirreled away that "
"knowledge deep in its bowels.  I have no way to remove it, and if we "
"were to continue nothing would link properly.  So, we're going to "
"stop here.  You should do 'module unload cray-libsci', delete all "
"the temporary CMake files (including CMakeCache.txt and the "
"directory CMakeFiles), and start over. You should also email the "
"system admins to get them to fix the linux-cc wrapper script. Sorry..."
"(If you're *not* on a broken Cray system, we've written an "
"overly-broad test for this and you should Email eisen@cc.gatech.edu "
"instead.)")
 endif()
ENDFUNCTION()


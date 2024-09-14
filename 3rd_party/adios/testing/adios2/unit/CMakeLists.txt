#------------------------------------------------------------------------------#
#Distributed under the OSI - approved Apache License, Version 2.0. See
#accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

gtest_add_tests_helper(ChunkV MPI_NONE "" Unit. "")
if(UNIX)
  gtest_add_tests_helper(PosixTransport MPI_NONE "" Unit. "")
endif()


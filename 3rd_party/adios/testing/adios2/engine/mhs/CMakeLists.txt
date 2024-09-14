#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

gtest_add_tests_helper(SingleRank MPI_NONE Mhs Engine.MHS. "")
gtest_add_tests_helper(MultiRank MPI_ONLY Mhs Engine.MHS. "")
gtest_add_tests_helper(MultiReader MPI_ONLY Mhs Engine.MHS. "")

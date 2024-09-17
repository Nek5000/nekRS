#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

gtest_add_tests_helper(WriteTypes MPI_ALLOW BP Bindings.C. "")
gtest_add_tests_helper(WriteReadMultiblock MPI_ALLOW BP Bindings.C. "")
gtest_add_tests_helper(NullWriteRead MPI_ALLOW "" Bindings.C. "")
gtest_add_tests_helper(WriteAggregateReadLocal MPI_ONLY BP Bindings.C. "")
gtest_add_tests_helper(AvailableVariablesAttribites MPI_ONLY BP Bindings.C. "")
gtest_add_tests_helper(MemorySpace MPI_NONE BP Bindings.C. "")

#!/usr/bin/env python
#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# TestBPWriteReadString.py: test writing/reading Python string type
# in ADIOS2 File Write
#  Created on: Oct 19, 2020
#      Author: Dmitry Ganyushin ganyushindi@ornl.gov
import unittest
from mpi4py import MPI
from adios2 import Stream

N_STEPS = 3


class test_adios_write_read_string_full_api(unittest.TestCase):
    def test_write_read_string_high_api(self):
        comm = MPI.COMM_WORLD
        theString = "hello adios"
        bpFilename = "string_test_highAPI.bp"
        varname = "mystringvar"

        with Stream(bpFilename, "w", comm=comm) as s:
            for step in s.steps(N_STEPS):
                s.write(varname, theString + str(step.current_step()))

        with Stream(bpFilename, "r", comm=comm) as s:
            for _ in s.steps():
                step = s.current_step()
                result = s.read(varname)
                self.assertEqual(result, theString + str(step))

    def test_read_strings_all_steps(self):
        comm = MPI.COMM_WORLD
        fileName = "string_test_all.bp"
        with Stream(fileName, "w", comm=comm) as s:
            i = 0
            for _ in s.steps(N_STEPS):
                s.write("string_variable", "written {}".format(i))
                i += 1

        # with Stream(fileName, "rra", comm = comm) as s:
        #    n = s.num_steps()
        #    name = "string_variable"
        #    result = s.read_string(name, 0, n)
        #    expected_str = ["written {}".format(i) for i in range(n)]
        #    self.assertEqual(result, expected_str)


if __name__ == "__main__":
    unittest.main()

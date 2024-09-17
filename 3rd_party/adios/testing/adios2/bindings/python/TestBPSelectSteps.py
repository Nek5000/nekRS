#!/usr/bin/env python
#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# TestBPSelectSteps_nompi.py: test step selection by reading in Python
# in ADIOS2 File Write
#  Created on: Jan 29, 2021
#      Author: Dmitry Ganyushin ganyushindi@ornl.gov
import unittest
import numpy as np
from mpi4py import MPI
import adios2.bindings as adios2

TESTDATA_FILENAME = "steps_int32.bp"

# MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
Nx = 8
shape = [size * Nx]
start = [rank * Nx]
count = [Nx]


class TestAdiosSelectSteps(unittest.TestCase):
    def setUp(self):
        total_steps = 10
        with adios2.open(TESTDATA_FILENAME, "w", comm) as fh:
            for i in range(total_steps):
                fh.write("step", np.full((Nx), i, dtype=np.int32), shape, start, count)
                fh.end_step()

    def test_select_steps_reading_fullAPI(self):
        selected_steps = [3, 5, 7]
        param_string = ",".join([str(i) for i in selected_steps])
        adios = adios2.ADIOS()
        ioReadBP = adios.DeclareIO("hellopy")
        ioReadBP.SetParameter(TESTDATA_FILENAME, param_string)
        fh = ioReadBP.Open(TESTDATA_FILENAME, adios2.Mode.ReadRandomAccess, comm)
        var = ioReadBP.InquireVariable("step")
        var.SetSelection([[0], [size * Nx]])
        var.SetStepSelection([0, len(selected_steps)])
        data = np.zeros((len(selected_steps), size * Nx), dtype=np.int32)
        fh.Get(var, data)
        fh.PerformGets()
        self.assertTrue(
            all(
                [
                    list(data[i]) == [selected_steps[i] for x in range(len(data[i]))]
                    for i in range(len(selected_steps))
                ]
            )
        )


if __name__ == "__main__":
    unittest.main()

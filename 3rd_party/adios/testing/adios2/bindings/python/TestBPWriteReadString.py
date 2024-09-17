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
import adios2.bindings as adios2

N_STEPS = 3


class TestAdiosWriteReadStringfullAPI(unittest.TestCase):
    def test_write_read_string_fullAPI(self):
        comm = MPI.COMM_WORLD
        theString = "hello adios"
        bpFilename = "string_test_fullAPI.bp"
        varname = "mystringvar"
        adios = adios2.ADIOS(comm)
        ioWrite = adios.DeclareIO("ioWriter")
        adEngine = ioWrite.Open(bpFilename, adios2.Mode.Write)
        varMyString = ioWrite.DefineVariable(varname)
        for step in range(N_STEPS):
            adEngine.BeginStep()
            adEngine.Put(varMyString, theString + str(step))
            adEngine.EndStep()
        adEngine.Close()

        ioRead = adios.DeclareIO("ioReader")
        adEngine = ioRead.Open(bpFilename, adios2.Mode.Read)
        for step in range(N_STEPS):
            adEngine.BeginStep()
            varReadMyString = ioRead.InquireVariable(varname)
            result = adEngine.Get(varReadMyString)
            adEngine.EndStep()
            self.assertEqual(result, theString + str(step))
        adEngine.Close()


if __name__ == "__main__":
    unittest.main()

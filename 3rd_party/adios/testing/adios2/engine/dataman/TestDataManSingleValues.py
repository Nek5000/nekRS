#!/usr/bin/env python
#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# TestDataMan1D.py: test for 1D data transfer by reading in Python
#  Created on: March 3, 2023
#      Author: Dmitry Ganyushin ganyushindi@ornl.gov
from multiprocessing import Process
import unittest
import numpy as np
import adios2


class TestDataMan1D(unittest.TestCase):

    def setUp(self):
        self.conf = {
            "IPAddress": "127.0.0.1",
            "Port": "12306",
            "Timeout": "5",
            "TransportMode": "reliable",
            "RendezvousReaderCount": "1",
        }
        self.Nx = 1
        self.fill_value = 1.0
        self.shape = [self.Nx]

    def test_run(self):

        s = Process(target=self.thread_send)
        r = Process(target=self.thread_receive)

        s.start()
        r.start()

        r.join()
        s.join()

    def thread_send(self):
        data = np.full(shape=self.shape, fill_value=self.fill_value)
        shape = data.shape
        count = shape
        start = (0,) * len(shape)

        adios_io = adios2.ADIOS()
        wan = adios_io.DeclareIO("Server")
        wan.SetEngine("Dataman")

        wan.SetParameters(self.conf)
        writer = wan.Open("testdata", adios2.Mode.Write)
        sendbuffer = wan.DefineVariable("np_data", data, shape,
                                        start, count, adios2.ConstantDims)
        self.assertIsNotNone(sendbuffer)
        if sendbuffer:
            writer.BeginStep()
            writer.Put(sendbuffer, data, adios2.Mode.Deferred)
            writer.EndStep()
        else:
            raise ValueError("DefineVariable failed")

        writer.Close()

    def thread_receive(self):
        data = np.zeros(shape=self.shape)
        adios_io = adios2.ADIOS()
        wan = adios_io.DeclareIO("Client")
        wan.SetEngine("Dataman")
        wan.SetParameters(self.conf)
        reader = wan.Open("testdata", adios2.Mode.Read)
        while True:
            stepStatus = reader.BeginStep()
            if stepStatus == adios2.StepStatus.OK:
                recvar = wan.InquireVariable("np_data")
                self.assertIsNotNone(recvar)
                reader.Get(recvar, data, adios2.Mode.Sync)

            elif stepStatus == adios2.StepStatus.EndOfStream:
                break
            else:
                raise StopIteration()
            reader.EndStep()
        reader.Close()
        self.assertTrue(all([data[i] == self.fill_value for i
                             in range(len(data))]))


if __name__ == '__main__':
    unittest.main()

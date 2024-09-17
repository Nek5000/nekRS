from adios2.adios import Adios

import adios2.bindings as bindings

import unittest
import numpy as np


class Test_attribute(unittest.TestCase):
    def test_create_write(self):
        adios = Adios()
        with adios.declare_io("BPWriter") as io:
            ts = io.define_attribute("timestamp", "20231122")
            self.assertEqual(ts.name(), "timestamp")
            self.assertEqual(ts.data_string(), "20231122")

    def test_create_reader(self):
        adios = Adios()
        with adios.declare_io("BPWriter") as io:
            ts = io.define_attribute("timestamp", "20231122")
            self.assertEqual(ts.name(), "timestamp")
            self.assertEqual(ts.data_string(), "20231122")

    def test_create_write_ndarray(self):
        adios = Adios()
        with adios.declare_io("BPWriter") as io:
            arr = np.array([2023, 11, 22])
            ts = io.define_attribute("timestamp", arr)
            self.assertEqual(ts.name(), "timestamp")
            self.assertTrue(np.array_equal(ts.data(), [2023, 11, 22]))


if __name__ == "__main__":
    unittest.main()

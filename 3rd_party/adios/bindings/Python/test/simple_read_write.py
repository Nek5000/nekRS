from adios2 import *
import adios2.bindings as bindings

import sys
import unittest


DATA_FILENAME = "hello-world-py.bp"

class TestSimpleReadWrite(unittest.TestCase):

    def _write(self, ad, greeting):
        """write a string to a bp file"""
        io = ad.declare_io("hello-world-writer")
        var_greeting = io.define_variable("Greeting")
        w = io.open(DATA_FILENAME, bindings.Mode.Write)
        w.begin_step()
        w.put(var_greeting, greeting)
        w.end_step()
        w.close()
        return 0

    def _read(self, ad):
        """read a string from to a bp file"""
        io = ad.declare_io("hello-world-reader")
        r = io.open(DATA_FILENAME, bindings.Mode.Read)
        r.begin_step()
        var_greeting = io.inquire_variable("Greeting")
        message = r.get(var_greeting)
        r.end_step()
        r.close()
        return message

    def test_simple_read_write(self):
        """driver function"""
        print("ADIOS2 version {0}".format(adios2.__version__))
        ad = adios2.Adios()
        greeting = "Hello World from ADIOS2"
        self._write(ad, greeting)
        message = self._read(ad)
        print("{}".format(message))
        self.assertEqual(greeting, message)
        return 0


if __name__ == '__main__':
    unittest.main()

from adios2.adios import Adios

import adios2.bindings as bindings

import unittest
import numpy as np


class TestOperator(unittest.TestCase):
    def test_operator_basic(self):
        adios = Adios()
        op1 = adios.define_operator("noop", "null")
        with adios.declare_io("BPWriter") as io:
            temps = io.define_variable("temps", np.empty([4], dtype=np.int64))
            temps.add_operation(op1)
            with io.open("pythontestvariable.bp", bindings.Mode.Write) as engine:
                temps_measures = np.array([35, 40, 30, 45], dtype=np.int64)
                engine.put(temps, temps_measures)

        op2 = adios.define_operator("noop2", "null")
        with adios.declare_io("BPReader") as reader:
            with reader.open("pythontestvariable.bp", bindings.Mode.Read) as engine:
                engine.begin_step()
                temps = reader.inquire_variable("temps")
                temps.add_operation(op2)
                engine.end_step()

    def test_operator_params(self):
        adios = Adios()
        op = adios.define_operator("noop", "null")
        op.set_parameter("speed", "best")
        op.get_parameters()
        self.assertTrue("speed" in op.get_parameters())
        self.assertEqual(op.get_parameters()["speed"], "best")


if __name__ == "__main__":
    unittest.main()

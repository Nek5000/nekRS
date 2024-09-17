from adios2.adios import Adios

import adios2.bindings as bindings

import unittest
import numpy as np


class TestVariable(unittest.TestCase):
    def test_create_write(self):
        adios = Adios()
        with adios.declare_io("BPWriter") as io:
            temps = io.define_variable("temps", np.empty([4], dtype=np.int64))
            with io.open("pythontestvariable.bp", bindings.Mode.Write) as engine:
                temps_measures = np.array([35, 40, 30, 45], dtype=np.int64)
                engine.put(temps, temps_measures)
                self.assertEqual(temps.name(), "temps")
                self.assertEqual(temps.block_id(), 0)
                self.assertEqual(temps.count(), [])
                self.assertEqual(temps.shape(), [])
                self.assertEqual(temps.sizeof(), 8)

    def test_create_reader(self):
        adios = Adios()
        with adios.declare_io("BPWriter") as io:
            temps = io.define_variable("temps", np.empty([4], dtype=np.int64))
            with io.open("pythontestvariable.bp", bindings.Mode.Write) as engine:
                temps_measures = np.array([35, 40, 30, 45], dtype=np.int64)
                engine.put(temps, temps_measures)

        with adios.declare_io("BPReader") as reader:
            with reader.open("pythontestvariable.bp", bindings.Mode.Read) as engine:
                engine.begin_step()
                temps = reader.inquire_variable("temps")
                engine.end_step()

                self.assertEqual(temps.name(), "temps")
                self.assertEqual(temps.block_id(), 0)
                self.assertEqual(temps.count(), [])
                self.assertEqual(temps.sizeof(), 8)
                self.assertEqual(temps.steps(), 1)
                self.assertEqual(temps.steps_start(), 0)

    def test_operators(self):
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


if __name__ == "__main__":
    unittest.main()

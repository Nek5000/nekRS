from adios2.adios import Adios

import adios2.bindings as bindings

import unittest
import numpy as np


class TestEngine(unittest.TestCase):
    def test_close(self):
        adios = Adios()
        io = adios.declare_io("BPWriter")
        with io.open("pythontestengine.bp", bindings.Mode.Write):
            pass

    def test_put(self):
        adios = Adios()
        with adios.declare_io("BPWriter") as io:
            pressure = io.define_variable("pressure")
            temps = io.define_variable("temps", np.empty([1], dtype=np.int64))
            with io.open("pythontestengine.bp", bindings.Mode.Write) as engine:
                engine.put(pressure, "35PSI")
                temps_measures = np.array([35, 40, 30, 45], dtype=np.int64)
                engine.put(temps, temps_measures)

    def test_get(self):
        adios = Adios()
        with adios.declare_io("BPWriter") as io:
            pressure = io.define_variable("pressure")
            temps = io.define_variable(
                name="temps",
                content=np.empty([1], dtype=np.int64),
                start=[0],
                shape=[4],
                count=[4],
            )
            with io.open("pythontestengine.bp", bindings.Mode.Write) as engine:
                engine.put(pressure, "35PSI")
                temps_measures = np.array([35, 40, 30, 45], dtype=np.int64)
                engine.put(temps, temps_measures)

        with adios.declare_io("BPReader") as reader:
            with reader.open("pythontestengine.bp", bindings.Mode.Read) as engine:
                engine.begin_step()
                pressure = reader.inquire_variable("pressure")
                temps = reader.inquire_variable("temps")
                pressure_reading = engine.get(pressure)
                temps_reading = np.empty([4], dtype=np.int64)
                engine.get(temps, temps_reading)
                info = engine.blocks_info("temps", 0)
                engine.end_step()
                self.assertEqual(pressure_reading, "35PSI")
                self.assertTrue(np.array_equal(temps_reading, np.array([35, 40, 30, 45])))

                self.assertIsNot(info, None)
                self.assertEqual(info[0]["Start"], "0")
                self.assertEqual(info[0]["Count"], "4")
                self.assertEqual(info[0]["WriterID"], "0")

    def test_steps(self):
        adios = Adios()
        with adios.declare_io("BPWriter") as io:
            pressure = io.define_variable("pressure")
            with io.open("pythontestengine.bp", bindings.Mode.Write) as engine:
                for step in range(0, 10):
                    engine.begin_step()
                    self.assertTrue(engine.between_step_pairs())
                    engine.put(pressure, f"{step}PSI")
                    engine.end_step()
                    self.assertFalse(engine.between_step_pairs())

        with adios.declare_io("BPReader") as reader:
            with reader.open("pythontestengine.bp", bindings.Mode.Read) as engine:
                for i in range(0, engine.steps()):
                    engine.begin_step()
                    pressure = reader.inquire_variable("pressure")
                    pressure_reading = engine.get(pressure)
                    self.assertEqual(pressure_reading, f"{i}PSI")
                    engine.end_step()

    def test_blockinfo(self):
        adios = Adios()
        with adios.declare_io("BPWriter") as io:
            temps = io.define_variable("temps", np.empty([4], dtype=np.int64))
            with io.open("pythontestengine.bp", bindings.Mode.Write) as engine:
                temps_measures = np.array([35, 40, 30, 45], dtype=np.int64)
                engine.put(temps, temps_measures)

        with adios.declare_io("BPReader") as reader:
            with reader.open("pythontestengine.bp", bindings.Mode.Read) as engine:
                engine.begin_step()
                all_blocks = engine.all_blocks_info("temps")
                temps = reader.inquire_variable("temps")
                info = engine.blocks_info("temps", 0)
                engine.end_step()

                self.assertIsNot(info, None)
                self.assertEqual(info[0]["Start"], "0")
                self.assertEqual(info[0]["Count"], "0")
                self.assertEqual(info[0]["WriterID"], "0")
                self.assertEqual(info, all_blocks[0])


if __name__ == "__main__":
    unittest.main()

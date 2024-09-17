from adios2.adios import Adios
import adios2.bindings as bindings

import unittest


class TestIO(unittest.TestCase):
    def test_io_empty(self):
        adios = Adios()
        adios.declare_io("BPWriter")

    def test_io_define_attribute(self):
        adios = Adios()
        io = adios.declare_io("BPWriter")
        ts = io.define_attribute("timestamp", "20231122")
        self.assertIsNot(ts, None)

    def test_io_inquire_attribute(self):
        adios = Adios()
        io = adios.declare_io("BPWriter")
        ts = io.define_attribute("timestamp", "20231122")
        coords = io.define_attribute("coords", "43N74W")
        x = io.inquire_attribute("coords")
        self.assertNotEqual(ts, coords)
        self.assertNotEqual(ts, x)
        self.assertEqual(coords, x)
        self.assertIs(io.inquire_attribute("notanattribute"), None)

    def test_available_attributes(self):
        adios = Adios()
        io = adios.declare_io("BPWriter")
        io.define_attribute("timestamp", "20231122")
        io.define_attribute("stringarray", ["one", "two", "three"])
        io.define_attribute("afloat", 3.14)
        io.define_attribute("floatarray", [3.14, 6.28])
        attrs = io.available_attributes()
        print("Available attributes:")
        for aname, ainfo in attrs.items():
            print(f"    {aname}:\t{ainfo}")
        self.assertEqual(attrs["timestamp"]["Value"], '"20231122"')
        self.assertEqual(attrs["stringarray"]["Value"], '{ "one", "two", "three" }')
        self.assertFalse("coords" in attrs)

    def test_remove_attribute(self):
        adios = Adios()
        io = adios.declare_io("BPWriter")
        io.define_attribute("timestamp", "20231122")
        io.remove_attribute("timestamp")
        self.assertIs(io.inquire_attribute("timestamp"), None)

    def test_remove_all_attribute(self):
        adios = Adios()
        io = adios.declare_io("BPWriter")
        io.define_attribute("timestamp", "20231122")
        io.define_attribute("coords", "43N74W")
        io.remove_all_attributes()
        self.assertIs(io.inquire_attribute("timestamp"), None)
        self.assertIs(io.inquire_attribute("coords"), None)

    def test_io_define_variable(self):
        adios = Adios()
        io = adios.declare_io("BPWriter")
        temp = io.define_variable("temp")
        self.assertNotEqual(temp, None)

    def test_io_inquire_variable(self):
        adios = Adios()
        io = adios.declare_io("BPWriter")
        temp = io.define_variable("temp")
        presure = io.define_variable("pressure")
        x = io.inquire_variable("pressure")
        self.assertNotEqual(temp, presure)
        self.assertNotEqual(temp, x)
        self.assertEqual(presure, x)

    def test_available_variable(self):
        adios = Adios()
        io = adios.declare_io("BPWriter")
        io.define_variable("temp")
        io.inquire_variable("temp")
        self.assertIs(io.inquire_attribute("pressure"), None)

    def test_remove_variable(self):
        adios = Adios()
        io = adios.declare_io("BPWriter")
        io.define_variable("temp")
        io.remove_variable("temp")
        self.assertIs(io.inquire_attribute("temp"), None)

    def test_remove_all_variable(self):
        adios = Adios()
        io = adios.declare_io("BPWriter")
        io.define_variable("temp")
        io.define_variable("pressure")
        io.remove_all_variables()
        self.assertIs(io.inquire_attribute("pressure"), None)
        self.assertIs(io.inquire_attribute("temp"), None)

    def test_open_engine(self):
        adios = Adios()
        io = adios.declare_io("BPWriter")
        io.set_engine("BPFile")
        io.set_parameter("threads", "2")
        io.set_parameters({"AsyncOpen": "On", "MaxOpenFilesAtOnce": "512"})
        io.add_transport("File", {"Library": "POSIX"})
        engine = io.open("pythontest.bp", bindings.Mode.Write)
        self.assertNotEqual(engine, None)
        engine.close()


if __name__ == "__main__":
    unittest.main()

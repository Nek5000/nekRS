from adios2.adios import Adios
import unittest


class Test_adios(unittest.TestCase):
    def test_define_operator(self):
        adios = Adios()
        op = adios.define_operator("op", "null")
        self.assertNotEqual(op, None)

    def test_inquiry_operator(self):
        adios = Adios()
        op1 = adios.define_operator("op1", "null")
        op2 = adios.define_operator("op2", "null")
        op_x = adios.inquire_operator("op2")
        self.assertNotEqual(op1, op2)
        self.assertEqual(op2, op_x)
        self.assertNotEqual(op1, op_x)

    def test_declare_io(self):
        adios = Adios()
        io = adios.declare_io("BPWriter")
        self.assertNotEqual(io, None)

    def test_at_io(self):
        adios = Adios()
        io = adios.declare_io("BPWriter")
        reader = adios.declare_io("BPReader")
        x = adios.at_io("BPReader")
        self.assertNotEqual(io, reader)
        self.assertEqual(reader, x)
        self.assertNotEqual(io, x)

    def test_remove_io(self):
        adios = Adios()
        adios.declare_io("BPWriter")
        adios.remove_io("BPWriter")

    def test_remove_all_io(self):
        adios = Adios()
        adios.declare_io("BPWriter")
        adios.declare_io("BPReader")
        adios.remove_all_ios()


if __name__ == "__main__":
    unittest.main()

from adios2 import FileReader, Stream, LocalValueDim
from random import randint

import unittest


class TestFileReader(unittest.TestCase):
    def test_basic(self):
        with Stream("pythonfiletest.bp", "w") as s:
            for _ in s.steps(10):
                if s.current_step() == 0:
                    # String
                    s.write("Outlook", "Good")
                    # Global Value single step
                    s.write("LocationID", 42)
                # Global Array
                s.write(
                    "Temp",
                    content=[randint(15, 35), randint(15, 35), randint(15, 35)],
                    shape=[3],
                    start=[0],
                    count=[3],
                )
                # Local Value
                s.write("Wind", [5], shape=[LocalValueDim])
                # Local Array
                s.write("Coords", [38, -46], [], [], [2])
                # Global Value every step
                s.write("Hour", 8 + s.current_step())

        with FileReader("pythonfiletest.bp") as f:
            self.assertEqual(len(f.all_blocks_info("temp")), 10)
            outlook_var = f.inquire_variable("Outlook")
            outlook = f.read(outlook_var)
            self.assertEqual(outlook, "Good")

            for name in f.available_variables():
                var = f.inquire_variable(name)
                if not var.type() == "string":
                    if var.steps() > 1:
                        self.assertEqual(var.steps(), 10)
                    var.set_step_selection([0, var.steps()])
                    output = f.read(var)
                    print(f"var:{var.name()} output:{output}")

        with FileReader("pythonfiletest.bp") as f:
            output = f.read("LocationID")
            self.assertEqual(output.ndim, 0)
            self.assertEqual(output.size, 1)
            self.assertEqual(output, 42)

            output = f.read("Hour", step_selection=[0, 10])
            self.assertEqual(output.ndim, 1)
            self.assertEqual(output.size, 10)
            self.assertTrue((output == [8, 9, 10, 11, 12, 13, 14, 15, 16, 17]).all())

            output = f.read("Hour", step_selection=[5, 1])
            self.assertEqual(output.ndim, 1)
            self.assertEqual(output.size, 1)
            self.assertEqual(output, [13])

            output = f.read("Temp", step_selection=[0, 10])
            self.assertEqual(len(output), 30)
            print(f"var:Temp output:{output}")

            output = f.read("Temp", step_selection=[0, 5])
            self.assertEqual(len(output), 15)
            print(f"var:Temp output:{output}")

            output = f.read("Temp", start=[0], count=[2], step_selection=[0, 10])
            self.assertEqual(len(output), 20)
            print(f"var:Temp output:{output}")


if __name__ == "__main__":
    unittest.main()

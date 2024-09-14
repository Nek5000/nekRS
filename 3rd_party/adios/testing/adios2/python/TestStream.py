from adios2 import Stream, LocalValueDim
from random import randint

import unittest


class TestStream(unittest.TestCase):
    def test_basic(self):
        with Stream("pythonstreamtest.bp", "w") as s:
            for _ in s.steps(10):
                # Single value string
                s.write("Outlook", "Good")
                # Global array
                s.write(
                    "temp",
                    content=[randint(15, 35), randint(15, 35), randint(15, 35)],
                    shape=[3, 1],
                    start=[0, 0],
                    count=[3, 1],
                )
                # Local Value
                s.write("Wind", [5], shape=[LocalValueDim])
                # Local Array
                s.write("Coords", [38, -46], [], [], [2])

        with Stream("pythonstreamtest.bp", "r") as s:
            for _ in s.steps():
                for var_name in s.available_variables():
                    print(f"var:{var_name}\t{s.read(var_name)}")
                    self.assertEqual(s.read("Wind", block_id=0), 5)
                    self.assertEqual(s.read("Coords", block_id=0)[0], 38)
                    self.assertEqual(s.read("Coords", block_id=0)[1], -46)

    def test_start_count(self):
        with Stream("pythonstreamtest.bp", "w") as s:
            for _ in s.steps(10):
                # Global array
                s.write(
                    "temp",
                    content=[randint(15, 35), randint(15, 35), randint(15, 35)],
                    shape=[3],
                    start=[0],
                    count=[3],
                )

        with Stream("pythonstreamtest.bp", "r") as s:
            print(s.all_blocks_info("temp"))
            self.assertEqual(len(s.all_blocks_info("temp")), 10)
            for _ in s.steps():
                for var_name in s.available_variables():
                    print(f"var:{var_name}\t{s.read(var_name)}")
                    output = s.read("temp", start=[0], count=[2])
                    self.assertEqual(len(output), 2)


if __name__ == "__main__":
    unittest.main()

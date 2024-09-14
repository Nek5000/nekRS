#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# nptypes.py small test data for np types
#  Created on: Feb 2, 2017
#      Author: William F Godoy godoywf@ornl.gov

import numpy as np


class SmallTestData:
    def __init__(self):
        self.Nx = 10
        self.Str = "Hello ADIOS2 Python"
        self.I8 = np.array([0, 1, -2, 3, -4, 5, -6, 7, -8, 9], dtype=np.int8)
        self.I16 = np.array([512, 513, -510, 515, -508, 517, -506, 519, -504, 521], dtype=np.int16)
        self.I32 = np.array(
            [131072, 131073, -131070, 131075, -131068, 131077, -131066, 131079, -131064, 131081],
            dtype=np.int32,
        )
        self.I64 = np.array(
            [
                8589934592,
                8589934593,
                -8589934590,
                8589934595,
                -8589934588,
                8589934597,
                -8589934586,
                8589934599,
                -8589934584,
                8589934601,
            ],
            dtype=np.int64,
        )

        self.U8 = np.array([128, 129, 130, 131, 132, 133, 134, 135, 136, 137], dtype=np.uint8)
        self.U16 = np.array(
            [32768, 32769, 32770, 32771, 32772, 32773, 32774, 32775, 32776, 32777], dtype=np.uint16
        )
        self.U32 = np.array(
            [
                2147483648,
                2147483649,
                2147483650,
                2147483651,
                2147483652,
                2147483653,
                2147483654,
                2147483655,
                2147483656,
                2147483657,
            ],
            dtype=np.uint32,
        )
        self.U64 = np.array(
            [
                9223372036854775808,
                9223372036854775809,
                9223372036854775810,
                9223372036854775811,
                9223372036854775812,
                9223372036854775813,
                9223372036854775814,
                9223372036854775815,
                9223372036854775816,
                9223372036854775817,
            ],
            dtype=np.uint64,
        )

        self.R32 = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=np.float32)
        self.R64 = np.array([0, -1, -2, -3, -4, -5, -6, -7, -8, -9], dtype=np.float64)

    def update(self, rank, step, size):
        self.I8 += 1
        self.I16 += 1
        self.I32 += 1
        self.I64 += 1

        self.U8 += 1
        self.U16 += 1
        self.U32 += 1
        self.U64 += 1

        self.R32 += 1
        self.R64 += 1

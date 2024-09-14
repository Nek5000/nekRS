#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# nptypes.py small test data for np types
#  Created on: Feb 2, 2017
#      Author: William F Godoy godoywf@ornl.gov

import numpy as np
import cmath


class SmallTestData:
    def __init__(self, rank=0):
        self.Nx = 10
        self.Str = "Hello ADIOS2 Python"
        self.i8 = np.array([0, 1, -2, 3, -4, 5, -6, 7, -8, 9], dtype=np.int8)
        self.i16 = np.array([512, 513, -510, 515, -508, 517, -506, 519, -504, 521], dtype=np.int16)
        self.i32 = np.array(
            [131072, 131073, -131070, 131075, -131068, 131077, -131066, 131079, -131064, 131081],
            dtype=np.int32,
        )
        self.i64 = np.array(
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

        self.u8 = np.array([128, 129, 130, 131, 132, 133, 134, 135, 136, 137], dtype=np.uint8)
        self.u16 = np.array(
            [32768, 32769, 32770, 32771, 32772, 32773, 32774, 32775, 32776, 32777], dtype=np.uint16
        )
        self.u32 = np.array(
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
        self.u64 = np.array(
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

        self.r32 = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=np.float32)
        self.r64 = np.array([0, -1, -2, -3, -4, -5, -6, -7, -8, -9], dtype=np.float64)

        self.c64 = np.array(
            [
                complex(0, 1),
                complex(2, 3),
                complex(4, 5),
                complex(6, 7),
                complex(8, 9),
                complex(0, -1),
                complex(-2, -3),
                complex(-4, -5),
                complex(-6, -7),
                complex(-8, -9),
            ],
            dtype=np.complex128,
        )

        self.int_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        self.float_list = [0.0, 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9]
        self.complex_list = [
            complex(0.0, rank),
            complex(1.0, rank),
            complex(2.0, rank),
            complex(3.0, rank),
            complex(4.0, rank),
            complex(5.0, rank),
            complex(6.0, rank),
            complex(7.0, rank),
            complex(8.0, rank),
            complex(9.0, rank),
        ]

        self.i8 += rank
        self.i16 += rank
        self.i32 += rank
        self.i64 += rank
        self.u8 += rank
        self.u16 += rank
        self.u32 += rank
        self.u64 += rank
        self.r32 += rank
        self.r64 += rank
        self.c64 += complex(rank, rank)
        self.int_list = [i + rank for i in self.int_list]
        self.float_list = [f + rank for f in self.float_list]
        self.complex_list = [c + complex(rank, 0) for c in self.complex_list]

    def update(self, rank, step, size):
        self.i8 += 1
        self.i16 += 1
        self.i32 += 1
        self.i64 += 1
        self.u8 += 1
        self.u16 += 1
        self.u32 += 1
        self.u64 += 1
        self.r32 += 1.0
        self.r64 += 1.0
        self.c64 += complex(1.0, 1.0)
        self.int_list = [i + 1 for i in self.int_list]
        self.float_list = [i + 1.0 for i in self.float_list]
        self.complex_list = [i + complex(1, 0) for i in self.complex_list]

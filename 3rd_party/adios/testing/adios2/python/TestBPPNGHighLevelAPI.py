#!/usr/bin/env python

#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# TestBPPNGHighLevelAPI.py
#
#  Created on: June 7th, 2019
#      Author: William F Godoy

import numpy as np
import random
from mpi4py import MPI
from adios2 import Stream


def compress_png(compression_level):
    fname = "BPWRPNG_" + str(compression_level) + "_py.bp"
    Nx = 10
    Ny = 50
    channels = 3
    NSteps = 1

    # initialize values
    u32s = np.zeros([Nx, Ny], np.uint32)
    u8s = np.zeros([Nx, Ny, channels], np.uint8)

    value_ji = 0.0
    for i in range(0, Nx):
        for j in range(0, Ny):
            u32s[i][j] = value_ji
            u8s[i][j][0] = random.randrange(256)
            u8s[i][j][1] = random.randrange(256)
            u8s[i][j][2] = random.randrange(256)

            value_ji += 1.0

    # set global dimensions
    # MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    shape3D = [Nx * size, Ny, 3]
    start3D = [Nx * rank, 0, 0]
    count3D = [Nx, Ny, 3]

    shape2D = [Nx * size, Ny]
    start2D = [Nx * rank, 0]
    count2D = [Nx, Ny]

    # writer
    with Stream(fname, "w", comm=comm) as s:
        for step in s.steps(NSteps):
            s.write(
                "u8",
                u8s,
                shape3D,
                start3D,
                count3D,
                [
                    (
                        "PNG",
                        {
                            "bit_depth": "8",
                            "color_type": "PNG_COLOR_TYPE_RGB",
                            "compression_level": str(compression_level),
                        },
                    )
                ],
            )
            s.write(
                "u32",
                u32s,
                shape2D,
                start2D,
                count2D,
                [
                    (
                        "PNG",
                        {
                            "bit_depth": "8",
                            "color_type": "PNG_COLOR_TYPE_RGBA",
                            "compression_level": str(compression_level),
                        },
                    )
                ],
            )

    # reader
    with Stream(fname, "r", comm=comm) as s:
        for step in s.steps():
            in_u8s = step.read("u8", start3D, count3D)
            in_u32s = step.read("u32", start2D, count2D)

            for i in range(0, Nx):
                for j in range(0, Ny):
                    assert u32s[i][j] == in_u32s[i][j]
                    assert u8s[i][j][0] == in_u8s[i][j][0]
                    assert u8s[i][j][1] == in_u8s[i][j][1]
                    assert u8s[i][j][2] == in_u8s[i][j][2]


def main():
    compress_png(compression_level=1)
    compress_png(compression_level=4)
    compress_png(compression_level=9)


if __name__ == "__main__":
    main()

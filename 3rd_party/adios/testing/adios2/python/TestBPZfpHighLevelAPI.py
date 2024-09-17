#!/usr/bin/env python

#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# TestBPZfpHighLevelAPI.py
#
#  Created on: April 2nd, 2019
#      Author: William F Godoy

import numpy as np
from mpi4py import MPI
from adios2 import Stream


def CompressZfp2D(rate):
    fname = "BPWRZfp2D_" + str(rate) + "_py.bp"
    Nx = 100
    Ny = 50
    NSteps = 2

    # initialize values
    r32s = np.zeros([Ny, Nx], np.float32)
    r64s = np.zeros([Ny, Nx], np.float64)

    value_ji = 0.0
    for j in range(0, Ny):
        for i in range(0, Nx):
            r32s[j][i] = value_ji
            r64s[j][i] = value_ji
            value_ji += 1.0

    # set global dimensions
    # MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    shape = [Ny * size, Nx]
    start = [Ny * rank, 0]
    count = [Ny, Nx]

    # writer
    with Stream(fname, "w", comm=comm) as s:
        for _ in s.steps(NSteps):
            s.write("r32", r32s, shape, start, count, [("zfp", {"accuracy": str(rate)})])
            s.write("r64", r64s, shape, start, count, [("zfp", {"accuracy": str(rate)})])

    # reader
    with Stream(fname, "r", comm=comm) as s:
        for _ in s.steps():
            in_r32s = s.read("r32", start, count)
            in_r64s = s.read("r64", start, count)

            for j in range(0, Ny):
                for i in range(0, Nx):
                    assert abs(r32s[j][i] - in_r32s[j][i]) < 1e-4
                    assert abs(r64s[j][i] - in_r64s[j][i]) < 1e-4


def main():
    CompressZfp2D(rate=8)
    CompressZfp2D(rate=9)
    CompressZfp2D(rate=10)


if __name__ == "__main__":
    main()

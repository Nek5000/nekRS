#!/usr/bin/env python

#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# TestChangingShape.py
#
#  Created on: April 2nd, 2019
#      Author: Jeremy Logan

import numpy as np
from mpi4py import MPI
from adios2 import Stream

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Test data
nx = [10, 15]
data = [np.zeros(nx[0]), np.ones(nx[1])]
shape = [[size * nx[0]], [size * nx[1]]]
start = [[rank * nx[0]], [rank * nx[1]]]
count = [[nx[0]], [nx[1]]]

# Write different sized arrays as separate steps
with Stream("out.bp", "w", comm=comm) as s:
    s.begin_step()
    s.write("z", data[0], shape[0], start[0], count[0])
    s.end_step()
    s.begin_step()
    s.write("z", data[1], shape[1], start[1], count[1])
    s.end_step()

# Read back arrays
with Stream("out.bp", "r", comm=comm) as s:
    for step in s.steps():
        shape_z = int(step.available_variables()["z"]["Shape"])
        print(shape_z)
        assert shape_z == int(shape[step.current_step()][0])

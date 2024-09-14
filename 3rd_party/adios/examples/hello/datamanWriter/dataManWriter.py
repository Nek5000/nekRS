#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# dataManWriter.py
#
#  Created on: Sept 5, 2019
#      Author: Jason Wang <wangr1@ornl.gov>
#

from mpi4py import MPI
import numpy as np
from time import sleep
from adios2 import Adios, Stream

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

Nx = 10
Ny = 10
steps = 10000

count = [Nx, Ny]
start = [rank * Nx, 0]
shape = [size * Nx, Ny]

floatArray = np.zeros(count, dtype=np.float32)

for i in range(0, Nx):
    for j in range(0, Ny):
        floatArray[i, j] = (start[0] + i) * shape[1] + (j + start[1])

adios = Adios(comm)
io = adios.declare_io("whatever")
io.set_engine("DataMan")
io.set_parameters({"IPAddress": "127.0.0.1", "Port": "12306", "Timeout": "5"})

var = io.define_variable("FloatArray", floatArray, shape, start, count, True)

with Stream(io, "HelloDataMan", "w") as stream:
    for _ in stream.steps(steps):
        # imitating costly computation
        floatArray = floatArray + 1
        sleep(1.0)

        print("Step", stream.current_step(), floatArray)
        stream.write(var, floatArray)
        # Warning: the data of the current step is not published until
        # the next loop entry or the exit of the loop

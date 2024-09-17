#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# TestBPWriteRead2D.py
#
#
#  Created on: Mar 3rd, 2019
#      Author: Kai Germaschewski <kai.germaschewski@unh.edu>
#              William F Godoy godoywf@ornl.gov
#

from mpi4py import MPI
import numpy as np
import adios2.bindings as adios2

# MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# User data
Nx = 10
Ny = 10

count = [Nx, Ny]
start = [rank * Nx, 0]
shape = [size * Nx, Ny]

temperatures = np.zeros(count, dtype=np.int32)

for i in range(0, Nx):
    for j in range(0, Ny):
        temperatures[i, j] = (start[0] + i) * shape[1] + (j + start[1])

# print(temperatures)
# ADIOS2 read
adios = adios2.ADIOS(comm)
ioWrite = adios.DeclareIO("ioWriter")

varTemperature = ioWrite.DefineVariable(
    "temperature2D", temperatures, shape, start, count, adios2.ConstantDims
)

obpStream = ioWrite.Open("HeatMap2D_py.bp", adios2.Mode.Write)
obpStream.Put(varTemperature, temperatures)
obpStream.Close()


if rank == 0:
    # ADIOS2 read
    ioRead = adios.DeclareIO("ioReader")
    ibpStream = ioRead.Open("HeatMap2D_py.bp", adios2.Mode.ReadRandomAccess, MPI.COMM_SELF)
    var_inTemperature = ioRead.InquireVariable("temperature2D")

    if var_inTemperature is False:
        raise ValueError("var_inTemperature is False")

    assert var_inTemperature is not None
    readOffset = [2, 2]
    readSize = [4, 4]

    var_inTemperature.SetSelection([readOffset, readSize])
    inTemperatures = np.zeros(readSize, dtype=np.int32)
    ibpStream.Get(var_inTemperature, inTemperatures, adios2.Mode.Sync)
    ibpStream.Close()

    # print('Incoming temperature map\n', inTemperatures)
    expected = np.array(
        [[22, 23, 24, 25], [32, 33, 34, 35], [42, 43, 44, 45], [52, 53, 54, 55]], np.int32
    )
    assert np.array_equal(inTemperatures, expected)

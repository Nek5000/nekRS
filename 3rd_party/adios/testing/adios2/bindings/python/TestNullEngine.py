#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# TestNullEngine.py
#
#
# Created on : Apr 11th, 2019
# Author : William F Godoy godoywf @ornl.gov
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

# ADIOS write
adios = adios2.ADIOS(comm)
ioWrite = adios.DeclareIO("ioWriter")

varTemperature = ioWrite.DefineVariable(
    "temperature2D", temperatures, shape, start, count, adios2.ConstantDims
)
ioWrite.SetEngine("NULL")

nullWriter = ioWrite.Open("NULL_py.bp", adios2.Mode.Write)

assert nullWriter.Type() == "NullWriter"

status = nullWriter.BeginStep()
assert status == adios2.StepStatus.OK

nullWriter.Put(varTemperature, temperatures)
nullWriter.EndStep()
nullWriter.Close()

# ADIOS2 read
ioRead = adios.DeclareIO("ioReader")
ioRead.SetEngine("null")
nullReader = ioRead.Open("NULL_py.bp", adios2.Mode.Read, MPI.COMM_SELF)

assert nullReader.Type() == "NullReader"

inTemperatures = np.zeros(1, dtype=np.int32)

status = nullReader.BeginStep()
assert status == adios2.StepStatus.EndOfStream

var_inTemperature = ioRead.InquireVariable("temperature2D")

if var_inTemperature is True:
    raise ValueError("var_inTemperature is not False")

# nullReader.Get(var_inTemperature, inTemperatures)

nullReader.PerformGets()
nullReader.EndStep()
nullReader.Close()

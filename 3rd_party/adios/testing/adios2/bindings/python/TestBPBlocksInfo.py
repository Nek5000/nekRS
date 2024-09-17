#
# Distributed under the OSI - approved Apache License, Version 2.0. See
# accompanying file Copyright.txt for details.
#
# TestBPWriteRead2D.py
#
#
# Created on : Mar 3rd, 2019
# Author : Jeremy Logan <lot@ornl.gov>
#          Kai Germaschewski <kai.germaschewski@unh.edu>
#          William F Godoy godoywf@ornl.gov
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
    ioRead = adios.DeclareIO("ioReader")
    ibpStream = ioRead.Open("HeatMap2D_py.bp", adios2.Mode.ReadRandomAccess, MPI.COMM_SELF)
    var_inTemperature = ioRead.InquireVariable("temperature2D")

    info = ibpStream.BlocksInfo("temperature2D", 0)
    assert info is not None
    assert info[0]["Start"] == "0,0"
    assert info[0]["Count"] == "10,10"
    assert info[0]["WriterID"] == "0"

#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# bpWriter.py : only works with MPI version
#  Created on: Feb 2, 2017
#      Authors: William F Godoy godoywf@ornl.gov
#               Norbert Podhorszki pnorbert@ornl.gov
#
# We use adios2.Adios and adios2.IO classes to set up the IO parameters.
# For default IO, it is sufficient to use the adios2.Stream class alone.

from mpi4py import MPI
import numpy
import adios2

# MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# User data
myArray = numpy.array([0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
myArray = myArray + (rank * 10)
Nx = myArray.size

# ADIOS MPI Communicator
adios = adios2.Adios(config_file=None, comm=comm)

# ADIOS IO
bpIO = adios.declare_io("BPFile_N2N")

bpIOParams = {}
bpIOParams["Threads"] = "2"
bpIOParams["ProfileUnits"] = "Microseconds"
bpIOParams["InitialBufferSize"] = "17Kb"
bpIO.set_parameters(bpIOParams)

bpIO.add_transport("File", {"Library": "fstream"})
bpIO.set_engine("BPFile")
a = bpIO.adios()

# ADIOS output stream
with adios2.Stream(bpIO, "bpWriter-py.bp", "w", comm) as fh:
    fh.write("bpArray", myArray, [size * Nx], [rank * Nx], [Nx])
    fh.write("Nx", Nx)
    fh.write_attribute("size", Nx, "bpArray")
    fh.write_attribute("dimensions", ["Nx"], "bpArray")

# Read content:
# bpls -la bpWriter-py.bp
# bpls -la bpWriter-py.bp -d bpArray -n 10
# bpls -la bpWriter-py.bp -d bpArray -n 10 -D

from mpi4py import MPI
import numpy as np
from time import sleep
from adios2 import Stream, Adios, bindings

# MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# ADIOS MPI Communicator
adios = Adios(comm)

# ADIOS IO
io = adios.declare_io("myIO")
io.set_engine("Sst")

# Python MPI is usually not compatible with MPI used to build ADIOS
# Must use TCP based data transport, or RDMA (ucx, libfabric) if available
io.set_parameter("DataTransport", "WAN")

# note: we need to use np.float32 to be compatible with data from C++ writer
# using "float" works in Python only but leads to type mismatch with C++
myArray = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], dtype=np.float32)
myArray = 10.0 * rank + myArray
nx = len(myArray)
increment = nx * size * 1.0

with Stream(io, "helloSst", "w", comm) as stream:
    for _ in stream.steps(4):
        currentStep = stream.current_step()

        # imitating computation
        sleep(1.0)

        stream.write("bpFloats", myArray, [size * nx], [rank * nx], [nx])
        print("Rank=", rank, "loop index =", currentStep, "data =", myArray, flush=True)
        myArray += increment
        # Warning: the data of the current step is not published until
        # the next loop entry or the exit of the loop

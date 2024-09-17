from mpi4py import MPI
import numpy as np
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

with Stream(io, "helloSst", "r", comm) as stream:
    for _ in stream.steps():
        var = stream.inquire_variable("bpFloats")
        shape = var.shape()
        count = int(shape[0] / size)
        start = count * rank
        if rank == size - 1:
            count += shape[0] % size
        floatArray = stream.read("bpFloats", [start], [count])
        currentStep = stream.current_step()
        loopStep = stream.loop_index()
        print(
            "Rank=",
            rank,
            "loop index =",
            loopStep,
            "stream step =",
            currentStep,
            "data =",
            floatArray,
            flush=True,
        )

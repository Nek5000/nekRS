from mpi4py import MPI
import numpy as np
from adios2 import Stream

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

nx = 10
shape = [size * nx]
start = [rank * nx]
count = [nx]

temperature = np.zeros(nx, dtype=np.double)
pressure = np.ones(nx, dtype=np.double)
delta_time = 0.01
physical_time = 0.0
nsteps = 5

with Stream("cfd.bp", "w", comm) as s:
    # NSteps from application
    for _ in s.steps(nsteps):
        if rank == 0 and s.current_step() == 0:
            # write a Python integer
            s.write("nproc", size)

        # write a Python floating point value
        s.write("physical_time", physical_time)
        # temperature and pressure are numpy arrays
        s.write("temperature", temperature, shape, start, count)
        s.write_attribute("temperature/unit", "K")
        s.write("pressure", pressure, shape, start, count)
        s.write_attribute("pressure/unit", "Pa")
        physical_time += delta_time

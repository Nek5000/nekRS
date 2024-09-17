#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# dataManReader.py
#
#  Created on: Sept 5, 2019
#      Author: Jason Wang <wangr1@ornl.gov>
#

from mpi4py import MPI
import numpy as np
from adios2 import Adios, Stream
from adios2.bindings import StepStatus

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

adios = Adios(comm)
io = adios.declare_io("whatever")
io.set_engine("DataMan")
io.set_parameters({"IPAddress": "127.0.0.1", "Port": "12306", "Timeout": "5"})

with Stream(io, "HelloDataMan", "r") as stream:
    for _ in stream.steps():
        stepStatus = stream.step_status()
        print(f"Step status = {stepStatus}")
        var = io.inquire_variable("FloatArray")
        shape = var.shape()
        floatArray = stream.read(var)
        currentStep = stream.current_step()
        loopStep = stream.loop_index()
        print("Loop index =", loopStep, "stream step =", currentStep, "data =")
        print(floatArray)

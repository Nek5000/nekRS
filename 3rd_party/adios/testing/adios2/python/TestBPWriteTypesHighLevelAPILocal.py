#!/usr/bin/env python

# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# TestBPWriteTypesHighLevelAPILocal.py: test Python numpy types in ADIOS2 File
#                      Write/Read High-Level API for Local Arrays
#  Created on: March 12, 2018
#      Author: William F Godoy godoywf@ornl.gov

from adios2NPTypes import SmallTestData
from mpi4py import MPI
from adios2 import Stream


def check_array(np1, np2, hint):
    if not (np1 == np2).all():
        print("InData: " + str(np1))
        print("Data: " + str(np2))
        raise ValueError("Array read failed " + str(hint))


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Test data
data = SmallTestData(rank)
nx = data.Nx

shape = []
start = []
count = [nx]


# Writer
with Stream("types_np_local.bp", "w", comm=comm) as s:
    for step in s.steps(5):
        data.update(rank, step.current_step(), size)
        s.write("varI8", data.i8, shape, start, count)
        s.write("varI16", data.i16, shape, start, count)
        s.write("varI32", data.i32, shape, start, count)
        s.write("varI64", data.i64, shape, start, count)
        s.write("varU8", data.u8, shape, start, count)
        s.write("varU16", data.u16, shape, start, count)
        s.write("varU32", data.u32, shape, start, count)
        s.write("varU64", data.u64, shape, start, count)
        s.write("varR32", data.r32, shape, start, count)
        s.write("varR64", data.r64, shape, start, count)
        s.write("an_int_list", data.int_list, shape, start, count)
        s.write("a_float_list", data.float_list, shape, start, count)
        s.write("a_complex_list", data.complex_list, shape, start, count)
# Reader
data = SmallTestData(rank)

with Stream("types_np_local.bp", "r", comm=comm) as s:
    for fr_step in s.steps():
        step = fr_step.current_step()
        data.update(rank, step, size)
        indataI8 = fr_step.read("varI8", block_id=rank)
        indataI16 = fr_step.read("varI16", block_id=rank)
        indataI32 = fr_step.read("varI32", block_id=rank)
        indataI64 = fr_step.read("varI64", block_id=rank)
        indataU8 = fr_step.read("varU8", block_id=rank)
        indataU16 = fr_step.read("varU16", block_id=rank)
        indataU32 = fr_step.read("varU32", block_id=rank)
        indataU64 = fr_step.read("varU64", block_id=rank)
        indataR32 = fr_step.read("varR32", block_id=rank)
        indataR64 = fr_step.read("varR64", block_id=rank)

        in_int_list = fr_step.read("an_int_list", block_id=rank)

        print(f"step {step} rank {rank} I16={indataI16} data.I16 = {data.i16}", flush=True)

        check_array(indataI8, data.i8, "i8")
        check_array(indataI16, data.i16, "i16")
        check_array(indataI32, data.i32, "i32")
        check_array(indataI64, data.i64, "i64")
        check_array(indataU8, data.u8, "u8")
        check_array(indataU16, data.u16, "u16")
        check_array(indataU32, data.u32, "u32")
        check_array(indataU64, data.u64, "u64")
        check_array(indataR32, data.r32, "r32")
        check_array(indataR64, data.r64, "r64")

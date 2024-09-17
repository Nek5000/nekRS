#!/usr/bin/env python

#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# TestBPWriteTypes.py: test Python numpy types in ADIOS2 File Write
#  Created on: Feb 2, 2017
#      Author: William F Godoy godoywf@ornl.gov


from adios2NPTypes import SmallTestData
from mpi4py import MPI
import numpy as np
import adios2.bindings as adios2


def check_object(adios2_object, name):
    if adios2_object is None:
        raise ValueError(str(name) + " not found")


def check_name(name, name_list):
    if name not in name_list:
        raise ValueError(str(name) + " not found in list")


# MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
Nx = 8

# Start ADIOS
adios = adios2.ADIOS(comm)
ioWriter = adios.DeclareIO("writer")

shape = [size * Nx]
start = [rank * Nx]
count = [Nx]

data = SmallTestData()

# ADIOS Variable name, shape, start, offset, constant dims
# All local variables
varI8 = ioWriter.DefineVariable("varI8", data.I8, shape, start, count, adios2.ConstantDims)
varI16 = ioWriter.DefineVariable("varI16", data.I16, shape, start, count, adios2.ConstantDims)
varI32 = ioWriter.DefineVariable("varI32", data.I32, shape, start, count, adios2.ConstantDims)
varI64 = ioWriter.DefineVariable("varI64", data.I64, shape, start, count, adios2.ConstantDims)

varU8 = ioWriter.DefineVariable("varU8", data.U8, shape, start, count, adios2.ConstantDims)
varU16 = ioWriter.DefineVariable("varU16", data.U16, shape, start, count, adios2.ConstantDims)
varU32 = ioWriter.DefineVariable("varU32", data.U32, shape, start, count, adios2.ConstantDims)
varU64 = ioWriter.DefineVariable("varU64", data.U64, shape, start, count, adios2.ConstantDims)

varR32 = ioWriter.DefineVariable("varR32", data.R32, shape, start, count, adios2.ConstantDims)

varR64 = ioWriter.DefineVariable("varR64", data.R64, shape, start, count, adios2.ConstantDims)

attString = ioWriter.DefineAttribute("attrString", ["hello attribute"])
attI8 = ioWriter.DefineAttribute("attrI8", data.I8)

# ADIOS Engine
writer = ioWriter.Open("npTypes.bp", adios2.Mode.Write)

for i in range(0, 3):
    npi8 = np.full((Nx), i, dtype=np.int8)
    npi16 = np.full((Nx), i, dtype=np.int16)
    npi32 = np.full((Nx), i, dtype=np.int32)
    npi64 = np.full((Nx), i, dtype=np.int64)
    npu8 = np.full((Nx), i, dtype=np.uint8)
    npu16 = np.full((Nx), i, dtype=np.uint16)
    npu32 = np.full((Nx), i, dtype=np.uint32)
    npu64 = np.full((Nx), i, dtype=np.uint64)
    npr32 = np.full((Nx), i, dtype=np.float32)
    npr64 = np.full((Nx), i, dtype=np.float64)

    writer.BeginStep()
    writer.Put(varI8, npi8)
    writer.Put(varI16, npi16)
    writer.Put(varI32, npi32)
    writer.Put(varI64, npi64)

    writer.Put(varU8, npu8)
    writer.Put(varU16, npu16)
    writer.Put(varU32, npu32)
    writer.Put(varU64, npu64)

    writer.Put(varR32, npr32)
    writer.Put(varR64, npr64)
    writer.EndStep()

writer.Close()

# Start reader
ioReader = adios.DeclareIO("reader")

reader = ioReader.Open("npTypes.bp", adios2.Mode.ReadRandomAccess)

attrString = ioReader.InquireAttribute("attrString")
attrI8 = ioReader.InquireAttribute("attrI8")

varI8 = ioReader.InquireVariable("varI8")
varI16 = ioReader.InquireVariable("varI16")
varI32 = ioReader.InquireVariable("varI32")
varI64 = ioReader.InquireVariable("varI64")
varU8 = ioReader.InquireVariable("varU8")
varU16 = ioReader.InquireVariable("varU16")
varU32 = ioReader.InquireVariable("varU32")
varU64 = ioReader.InquireVariable("varU64")
varR32 = ioReader.InquireVariable("varR32")
varR64 = ioReader.InquireVariable("varR64")

check_object(attrString, "attrString")
check_object(attrString, "attrI8")

check_object(varI8, "varI8")
check_object(varI16, "varI16")
check_object(varI32, "varI32")
check_object(varI64, "varI64")
check_object(varU8, "varU8")
check_object(varU16, "varU16")
check_object(varU32, "varU32")
check_object(varU64, "varU64")
check_object(varR32, "varR32")
check_object(varR64, "varR64")


attr_names = ["attrString", "attrI8"]
var_names = [
    "varStr",
    "varI8",
    "varI16",
    "varI32",
    "varI64",
    "varU8",
    "varU16",
    "varU32",
    "varU64",
    "varR32",
    "varR64",
]

attributesInfo = ioReader.AvailableAttributes()
for name, info in attributesInfo.items():
    check_name(name, attr_names)
    if rank == 0:
        print("attribute_name: " + name)
        for key, value in info.items():
            print("\t" + key + ": " + value)
        print("\n")

variablesInfo = ioReader.AvailableVariables()
for name, info in variablesInfo.items():
    check_name(name, var_names)
    if rank == 0:
        print("variable_name: " + name)
        for key, value in info.items():
            print("\t" + key + ": " + value)
        print("\n")


if varI8 is not None:
    varI8.SetSelection([[0], [size * Nx]])
    varI8.SetStepSelection([0, 3])
    inI8 = np.zeros((3, size * Nx), dtype=np.int8)
    reader.Get(varI8, inI8)

if varI16 is not None:
    varI16.SetSelection([[0], [size * Nx]])
    varI16.SetStepSelection([0, 3])
    inI16 = np.zeros((3, size * Nx), dtype=np.int16)
    reader.Get(varI16, inI16)

if varI32 is not None:
    varI32.SetSelection([[0], [size * Nx]])
    varI32.SetStepSelection([0, 3])
    inI32 = np.zeros((3, size * Nx), dtype=np.int32)
    reader.Get(varI32, inI32)

if varI64 is not None:
    varI64.SetSelection([[0], [size * Nx]])
    varI64.SetStepSelection([0, 3])
    inI64 = np.zeros((3, size * Nx), dtype=np.int64)
    reader.Get(varI64, inI64)

if varU8 is not None:
    varU8.SetSelection([[0], [size * Nx]])
    varU8.SetStepSelection([0, 3])
    inU8 = np.zeros((3, size * Nx), dtype=np.uint8)
    reader.Get(varU8, inU8)

if varU16 is not None:
    varU16.SetSelection([[0], [size * Nx]])
    varU16.SetStepSelection([0, 3])
    inU16 = np.zeros((3, size * Nx), dtype=np.uint16)
    reader.Get(varU16, inU16)

if varU32 is not None:
    varU32.SetSelection([[0], [size * Nx]])
    varU32.SetStepSelection([0, 3])
    inU32 = np.zeros((3, size * Nx), dtype=np.uint32)
    reader.Get(varU32, inU32)

if varU64 is not None:
    varU64.SetSelection([[0], [size * Nx]])
    varU64.SetStepSelection([0, 3])
    inU64 = np.zeros((3, size * Nx), dtype=np.uint64)
    reader.Get(varU64, inU64)

if varR32 is not None:
    varR32.SetSelection([[0], [size * Nx]])
    varR32.SetStepSelection([0, 3])
    inR32 = np.zeros((3, size * Nx), dtype=np.float32)
    reader.Get(varR32, inR32)

if varR64 is not None:
    varR64.SetSelection([[0], [size * Nx]])
    varR64.SetStepSelection([0, 3])
    inR64 = np.zeros((3, size * Nx), dtype=np.float64)
    reader.Get(varR64, inR64)

reader.PerformGets()

for i in range(0, 3):
    for j in range(0, Nx):
        if inI8[i][j] != i:
            raise ValueError("failed reading I8")

        if inI16[i][j] != i:
            raise ValueError("failed reading I16")

        if inI32[i][j] != i:
            raise ValueError("failed reading I32")

        if inI64[i][j] != i:
            raise ValueError("failed reading I64")

        if inU8[i][j] != i:
            raise ValueError("failed reading U8")

        if inU16[i][j] != i:
            raise ValueError("failed reading U16")

        if inU32[i][j] != i:
            raise ValueError("failed reading U32")

        if inU64[i][j] != i:
            raise ValueError("failed reading U64")

        if inR32[i][j] != i:
            raise ValueError("failed reading R32")

        if inR64[i][j] != i:
            raise ValueError("failed reading R64")


# here tests reader data
reader.Close()

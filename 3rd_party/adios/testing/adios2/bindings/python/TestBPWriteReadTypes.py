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
import adios2.bindings as adios2


def check_object(adios2_object, name):
    if adios2_object is False:
        raise ValueError(str(name) + " not found")


def check_name(name, name_list):
    if name not in name_list:
        raise ValueError(str(name) + " not found in list")


def check_array(np1, np2, hint):
    if (np1 == np2).all() is False:
        print("InData: " + str(np1))
        print("Data: " + str(np2))
        raise ValueError("Array read failed " + str(hint))


# MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
Nx = 8

# list of tested attributes and variables
attr_names = [
    "attrString",
    "attrStringArray",
    "attrI8",
    "attrI16",
    "attrI32",
    "attrI64",
    "attrU8",
    "attrU16",
    "attrU32",
    "attrU64",
    "varR32/attrR32",
    "varR64::attrR64",
]
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

# Start ADIOS
adios = adios2.ADIOS(comm)
ioWriter = adios.DeclareIO("writer")

shape = [size * Nx]
start = [rank * Nx]
count = [Nx]

data = SmallTestData()

# ADIOS Variable name, shape, start, offset, constant dims
# All local variables
varStr = ioWriter.DefineVariable("varStr")

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

attString = ioWriter.DefineAttribute("attrString", "one")
attStringArray = ioWriter.DefineAttribute("attrStringArray", ["one", "two", "three"])
attI8 = ioWriter.DefineAttribute("attrI8", data.I8)
attI16 = ioWriter.DefineAttribute("attrI16", data.I16)
attI32 = ioWriter.DefineAttribute("attrI32", data.I32)
attI64 = ioWriter.DefineAttribute("attrI64", data.I64)
attU8 = ioWriter.DefineAttribute("attrU8", data.U8)
attU16 = ioWriter.DefineAttribute("attrU16", data.U16)
attU32 = ioWriter.DefineAttribute("attrU32", data.U32)
attU64 = ioWriter.DefineAttribute("attrU64", data.U64)
# add an attribute to a variable
attR32 = ioWriter.DefineAttribute("attrR32", data.R32, "varR32")
attR64 = ioWriter.DefineAttribute("attrR64", data.R64, "varR64", "::")

ioWriter.SetEngine("BPFile")
ioParams = {"Threads": "1", "InitialBufferSize": "17Kb"}
ioWriter.SetParameters(ioParams)

engineType = ioWriter.EngineType()
if engineType != "BPFile":
    raise ValueError(str(engineType) + " incorrect engine type, should be BPFile")

ioWriter.SetParameter("profileunits", "microseconds")
ioWriter.AddTransport("file")

ioParams = ioWriter.Parameters()
print("Final IO parameters")
for key, value in ioParams.items():
    print("\t" + key + ": " + value)

nsteps = 3

# ADIOS Engine
writer = ioWriter.Open("npTypes.bp", adios2.Mode.Write)

writer.LockWriterDefinitions()

for i in range(0, nsteps):
    data.update(rank, i, size)

    writer.BeginStep()
    writer.Put(varStr, data.Str)
    writer.Put(varI8, data.I8)
    writer.Put(varI16, data.I16)
    writer.Put(varI32, data.I32)
    writer.Put(varI64, data.I64)

    writer.Put(varU8, data.U8)
    writer.Put(varU16, data.U16)
    writer.Put(varU32, data.U32)
    writer.Put(varU64, data.U64)

    writer.Put(varR32, data.R32)
    writer.Put(varR64, data.R64)
    writer.EndStep()

writer.Close()

# Start reader
ioReader = adios.DeclareIO("reader")

reader = ioReader.Open("npTypes.bp", adios2.Mode.ReadRandomAccess)

attrString = ioReader.InquireAttribute("attrString")
attrStringArray = ioReader.InquireAttribute("attrStringArray")
attrI8 = ioReader.InquireAttribute("attrI8")
attrI16 = ioReader.InquireAttribute("attrI16")
attrI32 = ioReader.InquireAttribute("attrI32")
attrI64 = ioReader.InquireAttribute("attrI64")
attrU8 = ioReader.InquireAttribute("attrU8")
attrU16 = ioReader.InquireAttribute("attrU16")
attrU32 = ioReader.InquireAttribute("attrU32")
attrU64 = ioReader.InquireAttribute("attrU64")
attrR32 = ioReader.InquireAttribute("attrR32", "varR32")
attrR64 = ioReader.InquireAttribute("attrR64", "varR64", "::")

check_object(attrString, "attrString")
check_object(attrStringArray, "attrStringArray")
check_object(attrI8, "attrI8")
check_object(attrI16, "attrI16")
check_object(attrI32, "attrI32")
check_object(attrI64, "attrI64")
check_object(attrU8, "attrU8")
check_object(attrU16, "attrU16")
check_object(attrU32, "attrU32")
check_object(attrU64, "attrU64")
check_object(attrR32, "varR32/attrR32")
check_object(attrR64, "varR64::attrR64")

# alternative inquire format
attrR32 = ioReader.InquireAttribute("varR32/attrR32")
attrR64 = ioReader.InquireAttribute("varR64::attrR64")
check_object(attrR32, "varR32/attrR32")
check_object(attrR64, "varR64::attrR64")

attrStringData = attrString.DataString()
print(f"attrString = {attrStringData}", flush=True)
if attrStringData[0] != "one":
    raise ValueError("attrString failed")

attrStringData = attrStringArray.DataString()
print(f"attrStringArray = {attrStringData}", flush=True)
if attrStringData[0] != "one":
    raise ValueError("attrStringData[0] failed")
if attrStringData[1] != "two":
    raise ValueError("attrStringData[1] failed")
if attrStringData[2] != "three":
    raise ValueError("attrStringData[2] failed")

attrI8Data = attrI8.Data()
attrI16Data = attrI16.Data()
attrI32Data = attrI32.Data()
attrI64Data = attrI64.Data()
attrU8Data = attrU8.Data()
attrU16Data = attrU16.Data()
attrU32Data = attrU32.Data()
attrU64Data = attrU64.Data()
attrR32Data = attrR32.Data()
attrR64Data = attrR64.Data()

check_array(attrI8Data, data.I8, "I8")
check_array(attrI16Data, data.I16, "I16")
check_array(attrI32Data, data.I32, "I32")
check_array(attrI64Data, data.I64, "I64")
check_array(attrU8Data, data.U8, "U8")
check_array(attrU16Data, data.U16, "U16")
check_array(attrU32Data, data.U32, "U32")
check_array(attrU64Data, data.U64, "U64")
check_array(attrR32Data, data.R32, "R32")
check_array(attrR64Data, data.R64, "R64")

print("=========== Attributes ===========")
attributesInfo = ioReader.AvailableAttributes()
for name, info in attributesInfo.items():
    check_name(name, attr_names)
    if rank == 0:
        print("attribute_name: " + name)
        for key, value in info.items():
            print("\t" + key + ": " + value)
        print("\n")

print("=========== Available attributes of varR32 ===========")
attributesInfoR32 = ioReader.AvailableAttributes("varR32")
for name, info in attributesInfoR32.items():
    check_name(name, ["attrR32"])
    if rank == 0:
        print("attribute_name: " + name)
        for key, value in info.items():
            print("\t" + key + ": " + value)
        print("\n")

print("=========== Available attributes of varR64 ===========")
attributesInfoR32 = ioReader.AvailableAttributes("varR64", "::")
for name, info in attributesInfoR32.items():
    check_name(name, ["attrR64"])
    if rank == 0:
        print("attribute_name: " + name)
        for key, value in info.items():
            print("\t" + key + ": " + value)
        print("\n")

print("=========== Variables ===========")
varStr = ioReader.InquireVariable("varStr")
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

check_object(varStr, "varStr")
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

variablesInfo = ioReader.AvailableVariables()
for name, info in variablesInfo.items():
    check_name(name, var_names)
    if rank == 0:
        print("variable_name: " + name)
        for key, value in info.items():
            print("\t" + key + ": " + value)
        print("\n")


result = adios.RemoveIO("writer")
if result is False:
    raise ValueError("Could not remove IO writer")

assert reader.Steps() == nsteps


reader.Close()

ioReader.RemoveAllVariables()
varStr = ioReader.InquireVariable("varStr")
if varStr is True:
    raise ValueError("Could remove reader variables")

adios.RemoveAllIOs()
try:
    ioWriter = adios.DeclareIO("reader")
except ValueError:
    raise ValueError("Could not re-Declare IO reader")

#!/usr/bin/env python

#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# TestBPWriteTypes.py: test Python numpy types in ADIOS2 File Write
#  Created on: Feb 2, 2017
#      Author: William F Godoy godoywf@ornl.gov


from adios2NPTypes import SmallTestData
import adios2.bindings as adios2


# Test data
data = SmallTestData()

adios = adios2.ADIOS()

bpIO = adios.DeclareIO("NPTypes")

# ADIOS Variable name, shape, start, offset, constant dims
# All local variables
varI8 = bpIO.DefineVariable("varI8", data.I8, [], [], [data.I8.size], adios2.ConstantDims)
varI16 = bpIO.DefineVariable("varI16", data.I16, [], [], [data.I16.size], adios2.ConstantDims)
varI32 = bpIO.DefineVariable("varI32", data.I32, [], [], [data.I32.size], adios2.ConstantDims)
varI64 = bpIO.DefineVariable("varI64", data.I64, [], [], [data.I64.size], adios2.ConstantDims)

varU8 = bpIO.DefineVariable("varUI8", data.U8, [], [], [data.U8.size], adios2.ConstantDims)
varU16 = bpIO.DefineVariable("varUI16", data.U16, [], [], [data.U16.size], adios2.ConstantDims)
varU32 = bpIO.DefineVariable("varUI32", data.U32, [], [], [data.U32.size], adios2.ConstantDims)
varU64 = bpIO.DefineVariable("varUI64", data.U64, [], [], [data.U64.size], adios2.ConstantDims)

varR32 = bpIO.DefineVariable("varR32", data.R32, [], [], [data.R32.size], adios2.ConstantDims)

varR64 = bpIO.DefineVariable("varR64", data.R64, [], [], [data.R64.size], adios2.ConstantDims)


# ADIOS Engine
bpFileWriter = bpIO.Open("npTypes.bp", adios2.Mode.Write)

bpFileWriter.Put(varI8, data.I8)
bpFileWriter.Put(varI16, data.I16)
bpFileWriter.Put(varI32, data.I32)
bpFileWriter.Put(varI64, data.I64)

bpFileWriter.Put(varU8, data.U8)
bpFileWriter.Put(varU16, data.U16)
bpFileWriter.Put(varU32, data.U32)
bpFileWriter.Put(varU64, data.U64)

bpFileWriter.Put(varR32, data.R32)
bpFileWriter.Put(varR64, data.R64)

bpFileWriter.Close()

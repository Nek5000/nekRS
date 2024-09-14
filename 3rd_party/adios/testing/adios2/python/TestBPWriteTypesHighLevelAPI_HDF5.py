#!/usr/bin/env python

#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# TestBPWriteTypes.py: test Python numpy types in ADIOS2 File
#                      Write/Read High-Level API
#  Created on: March 12, 2018
#      Author: William F Godoy godoywf@ornl.gov

from adios2NPTypes import SmallTestData
import numpy as np
from adios2 import Adios, Stream

rank = 0
size = 1

# Test data
data = SmallTestData(rank)
nx = data.Nx

shape = [size * nx]
start = [rank * nx]
count = [nx]

adios = Adios()
io = adios.declare_io("writeh5")
io.set_engine("HDF5")

# Writer
with Stream(io, "types_np.h5", "w") as fw:
    for i in range(0, 5):
        data.update(rank, i, size)

        if rank == 0 and i == 0:
            fw.write("tag", "Testing ADIOS2 high-level API")
            fw.write("nx", data.Nx)
            fw.write("gvarI8", np.array(data.i8[0]))
            fw.write("gvarI16", np.array(data.i16[0]))
            fw.write("gvarI32", np.array(data.i32[0]))
            fw.write("gvarI64", np.array(data.i64[0]))
            fw.write("gvarU8", np.array(data.u8[0]))
            fw.write("gvarU16", np.array(data.u16[0]))
            fw.write("gvarU32", np.array(data.u32[0]))
            fw.write("gvarU64", np.array(data.u64[0]))
            fw.write("gvarR32", np.array(data.r32[0]))
            fw.write("gvarR64", np.array(data.r64[0]))

            # single value attributes
            fw.write_attribute("attrStr", "Testing single string attribute")
            fw.write_attribute("attrNx", data.Nx)
            fw.write_attribute("attrI8", np.array(data.i8[0]))
            fw.write_attribute("attrI16", np.array(data.i16[0]))
            fw.write_attribute("attrI32", np.array(data.i32[0]))
            fw.write_attribute("attrI64", np.array(data.i64[0]))
            fw.write_attribute("attrU8", np.array(data.u8[0]))
            fw.write_attribute("attrU16", np.array(data.u16[0]))
            fw.write_attribute("attrU32", np.array(data.u32[0]))
            fw.write_attribute("attrU64", np.array(data.u64[0]))
            fw.write_attribute("attrR32", np.array(data.r32[0]))
            fw.write_attribute("attrR64", np.array(data.r64[0]))

            fw.write_attribute("attrStrArray", ["string1", "string2", "string3"])
            fw.write_attribute("attrI8Array", data.i8)
            fw.write_attribute("attrI16Array", data.i16)
            fw.write_attribute("attrI32Array", data.i32)
            fw.write_attribute("attrI64Array", data.i64)
            fw.write_attribute("attrU8Array", data.u8)
            fw.write_attribute("attrU16Array", data.u16)
            fw.write_attribute("attrU32Array", data.u32)
            fw.write_attribute("attrU64Array", data.u64)
            fw.write_attribute("attrR32Array", data.r32)
            fw.write_attribute("attrR64Array", data.r64)

        fw.write("steps", "Step:" + str(i))
        fw.write("varI8", data.i8, shape, start, count)
        fw.write("varI16", data.i16, shape, start, count)
        fw.write("varI32", data.i32, shape, start, count)
        fw.write("varI64", data.i64, shape, start, count)
        fw.write("varU8", data.u8, shape, start, count)
        fw.write("varU16", data.u16, shape, start, count)
        fw.write("varU32", data.u32, shape, start, count)
        fw.write("varU64", data.u64, shape, start, count)
        fw.write("varR32", data.r32, shape, start, count)
        fw.write("varR64", data.r64, shape, start, count)

        if rank == 0 and i == 0:
            fw.write_attribute("varattrStrArray", ["varattr1", "varattr2", "varattr3"], "steps")
            fw.write_attribute("varattrI8Array", data.i8, "varI8")
            fw.write_attribute("varattrI16Array", data.i16, "varI16")
            fw.write_attribute("varattrI32Array", data.i32, "varI32")
            fw.write_attribute("varattrI64Array", data.i64, "varI64")
            fw.write_attribute("varattrU8Array", data.u8, "varU8")
            fw.write_attribute("varattrU16Array", data.u16, "varU16")
            fw.write_attribute("varattrU32Array", data.u32, "varU32")
            fw.write_attribute("varattrU64Array", data.u64, "varU64")
            fw.write_attribute("varattrR32Array", data.r32, "varR32")
            fw.write_attribute("varattrR64Array", data.r64, "varR64")
            fw.write_attribute("varattrR64Value", data.r64, "varR64")

        fw.end_step()

# Reader
data = SmallTestData(rank)

io = adios.declare_io("readh5")
io.set_engine("HDF5")

with Stream("types_np.h5", "r") as fr:
    for fr_step in fr:
        step = fr_step.current_step()
        data.update(rank, step, size)

        step_vars = fr_step.available_variables()

        #         for name, info in step_vars.items():
        #             print("variable_name: " + name)
        #             for key, value in info.items():
        #                 print("\t" + key + ": " + value)
        #             print("\n")

        if rank == 0 and step == 0:
            inTag = fr_step.read("tag")
            inNx = fr_step.read("nx")
            inI8 = fr_step.read("gvarI8")
            inI16 = fr_step.read("gvarI16")
            inI32 = fr_step.read("gvarI32")
            inI64 = fr_step.read("gvarI64")
            inU8 = fr_step.read("gvarU8")
            inU16 = fr_step.read("gvarU16")
            inU32 = fr_step.read("gvarU32")
            inU64 = fr_step.read("gvarU64")
            inR32 = fr_step.read("gvarR32")
            inR64 = fr_step.read("gvarR64")

            if inTag != "Testing ADIOS2 high-level API":
                print("InTag: " + str(inTag))
                raise ValueError("tag variable read failed")

            if inNx != nx:
                raise ValueError("tag variable read failed")

            if inI8 != data.i8[0]:
                raise ValueError("gvarI8 read failed")

            if inI16 != data.i16[0]:
                raise ValueError("gvarI16 read failed")

            if inI32 != data.i32[0]:
                raise ValueError("gvarI32 read failed")

            if inI64 != data.i64[0]:
                raise ValueError("gvarI64 read failed")

            if inU8 != data.u8[0]:
                raise ValueError("gvarU8 read failed")

            if inU16 != data.u16[0]:
                raise ValueError("gvarU16 read failed")

            if inU32 != data.u32[0]:
                raise ValueError("gvarU32 read failed")

            if inU64 != data.u64[0]:
                raise ValueError("gvarU64 read failed")

            if inR32 != data.r32[0]:
                raise ValueError("gvarR32 read failed")

            if inR64 != data.r64[0]:
                raise ValueError("gvarR64 read failed")

            # attributes
            inTag = fr_step.read_attribute("attrStr")
            inNx = fr_step.read_attribute("attrNx")
            inI8 = fr_step.read_attribute("attrI8")
            inI16 = fr_step.read_attribute("attrI16")
            inI32 = fr_step.read_attribute("attrI32")
            inI64 = fr_step.read_attribute("attrI64")
            inU8 = fr_step.read_attribute("attrU8")
            inU16 = fr_step.read_attribute("attrU16")
            inU32 = fr_step.read_attribute("attrU32")
            inU64 = fr_step.read_attribute("attrU64")
            inR32 = fr_step.read_attribute("attrR32")
            inR64 = fr_step.read_attribute("attrR64")

            if inTag != "Testing single string attribute":
                raise ValueError("attr string read failed")

            if inNx != nx:
                raise ValueError("attrI8 read failed")

            if inI8 != data.i8[0]:
                raise ValueError("attrI8 read failed")

            if inI16 != data.i16[0]:
                raise ValueError("attrI16 read failed")

            if inI32 != data.i32[0]:
                raise ValueError("attrI32 read failed")

            if inI64 != data.i64[0]:
                raise ValueError("attrI64 read failed")

            if inU8 != data.u8[0]:
                raise ValueError("attrU8 read failed")

            if inU16 != data.u16[0]:
                raise ValueError("attrU16 read failed")

            if inU32 != data.u32[0]:
                raise ValueError("attrU32 read failed")

            if inU64 != data.u64[0]:
                raise ValueError("attrU64 read failed")

            if inR32 != data.r32[0]:
                raise ValueError("attrR32 read failed")

            if inR64 != data.r64[0]:
                raise ValueError("attrR64 read failed")

            # Array attribute
            inTag = fr_step.read_attribute_string("attrStrArray")
            inI8 = fr_step.read_attribute("attrI8Array")
            inI16 = fr_step.read_attribute("attrI16Array")
            inI32 = fr_step.read_attribute("attrI32Array")
            inI64 = fr_step.read_attribute("attrI64Array")
            inU8 = fr_step.read_attribute("attrU8Array")
            inU16 = fr_step.read_attribute("attrU16Array")
            inU32 = fr_step.read_attribute("attrU32Array")
            inU64 = fr_step.read_attribute("attrU64Array")
            inR32 = fr_step.read_attribute("attrR32Array")
            inR64 = fr_step.read_attribute("attrR64Array")

            if inTag != ["string1", "string2", "string3"]:
                raise ValueError("attrStrArray read failed")

            if (inI8 == data.i8).all() is False:
                raise ValueError("attrI8 array read failed")

            if (inI16 == data.i16).all() is False:
                raise ValueError("attrI16 array read failed")

            if (inI32 == data.i32).all() is False:
                raise ValueError("attrI32 array read failed")

            if (inI64 == data.i64).all() is False:
                raise ValueError("attrI64 array read failed")

            if (inU8 == data.u8).all() is False:
                raise ValueError("attrU8 array read failed")

            if (inU16 == data.u16).all() is False:
                raise ValueError("attrU16 array read failed")

            if (inU32 == data.u32).all() is False:
                raise ValueError("attrU32 array read failed")

            if (inU64 == data.u64).all() is False:
                raise ValueError("attrU64 array read failed")

            if (inR32 == data.r32).all() is False:
                raise ValueError("attrR32 array read failed")

            if (inR64 == data.r64).all() is False:
                raise ValueError("attrR64 array read failed")

            inTags = fr_step.read_attribute_string("varattrStrArray", "steps")
            inI8 = fr_step.read_attribute("varattrI8Array", "varI8")
            in16 = fr_step.read_attribute("varattrI16Array", "varI16")
            inI32 = fr_step.read_attribute("varattrI32Array", "varI32")
            inI64 = fr_step.read_attribute("varattrI64Array", "varI64")
            inU8 = fr_step.read_attribute("varattrU8Array", "varU8")
            inU16 = fr_step.read_attribute("varattrU16Array", "varU16")
            inU32 = fr_step.read_attribute("varattrU32Array", "varU32")
            inU64 = fr_step.read_attribute("varattrU64Array", "varU64")
            inR32 = fr_step.read_attribute("varattrR32Array", "varR32")
            inR64 = fr_step.read_attribute("varattrR64Array", "varR64")

            if inTags != ["varattr1", "varattr2", "varattr3"]:
                print(inTags)
                raise ValueError("var attrStrArray read failed")

            if (inI8 == data.i8).all() is False:
                raise ValueError("var attrI8 array read failed")

            if (inI16 == data.i16).all() is False:
                raise ValueError("var attrI16 array read failed")

            if (inI32 == data.i32).all() is False:
                raise ValueError("var attrI32 array read failed")

            if (inI64 == data.i64).all() is False:
                raise ValueError("var attrI64 array read failed")

            if (inU8 == data.u8).all() is False:
                raise ValueError("var attrU8 array read failed")

            if (inU16 == data.u16).all() is False:
                raise ValueError("var attrU16 array read failed")

            if (inU32 == data.u32).all() is False:
                raise ValueError("var attrU32 array read failed")

            if (inU64 == data.u64).all() is False:
                raise ValueError("var attrU64 array read failed")

            if (inR32 == data.r32).all() is False:
                raise ValueError("var attrR32 array read failed")

            if (inR64 == data.r64).all() is False:
                raise ValueError("var attrR64 array read failed")

        stepStr = "Step:" + str(step)

        instepStr = fr_step.read("steps")
        if instepStr != stepStr:
            raise ValueError("steps variable read failed: " + instepStr + " " + stepStr)

        indataI8 = fr_step.read("varI8", start, count)
        indataI16 = fr_step.read("varI16", start, count)
        indataI32 = fr_step.read("varI32", start, count)
        indataI64 = fr_step.read("varI64", start, count)
        indataU8 = fr_step.read("varU8", start, count)
        indataU16 = fr_step.read("varU16", start, count)
        indataU32 = fr_step.read("varU32", start, count)
        indataU64 = fr_step.read("varU64", start, count)
        indataR32 = fr_step.read("varR32", start, count)
        indataR64 = fr_step.read("varR64", start, count)
        fr_step.end_step()

        if (indataI8 == data.i8).all() is False:
            raise ValueError("i8 array read failed")

        if (indataI16 == data.i16).all() is False:
            raise ValueError("i16 array read failed")

        if (indataI32 == data.i32).all() is False:
            raise ValueError("i32 array read failed")

        if (indataI64 == data.i64).all() is False:
            raise ValueError("i64 array read failed")

        if (indataU8 == data.u8).all() is False:
            raise ValueError("u8 array read failed")

        if (indataU16 == data.u16).all() is False:
            raise ValueError("u16 array read failed")

        if (indataU32 == data.u32).all() is False:
            raise ValueError("u32 array read failed")

        if (indataU64 == data.u64).all() is False:
            raise ValueError("u64 array read failed")

        if (indataR32 == data.r32).all() is False:
            raise ValueError("r32 array read failed")

        if (indataR64 == data.r64).all() is False:
            raise ValueError("r64 array read failed")

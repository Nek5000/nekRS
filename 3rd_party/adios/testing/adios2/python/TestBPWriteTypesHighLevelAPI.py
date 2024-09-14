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
from mpi4py import MPI
import numpy as np
from adios2 import Stream, LocalValueDim

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Test data
data = SmallTestData(rank)
nx = data.Nx

shape = [size * nx]
start = [rank * nx]
count = [nx]

# Writer
with Stream("types_np.bp", "w", comm=comm) as s:
    for step in s.steps(5):
        data.update(rank, step.current_step(), size)
        s.write("rank", np.array(rank), shape=[LocalValueDim])
        if rank == 0 and step.current_step() == 0:
            s.write("tag", "Testing ADIOS2 high-level API")
            s.write("nx", data.Nx)
            s.write("gvarI8", np.array(data.i8[0]))
            s.write("gvarI16", np.array(data.i16[0]))
            s.write("gvarI32", np.array(data.i32[0]))
            s.write("gvarI64", np.array(data.i64[0]))
            s.write("gvarU8", np.array(data.u8[0]))
            s.write("gvarU16", np.array(data.u16[0]))
            s.write("gvarU32", np.array(data.u32[0]))
            s.write("gvarU64", np.array(data.u64[0]))
            s.write("gvarR32", np.array(data.r32[0]))
            s.write("gvarR64", np.array(data.r64[0]))
            s.write("gvarC64", np.array(data.c64[0]))
            i = data.int_list[0]
            f = data.float_list[0]
            c = data.complex_list[0]
            print(f"type of i = {type(i)}, type of f = {type(f)}, type of c = {type(c)}")
            s.write("an_int_value", i)
            s.write("a_float_value", f)
            s.write("a_complex_value", c)
            s.write("an_int_list", data.int_list)
            s.write("a_float_list", data.float_list)
            s.write("a_complex_list", data.complex_list)

            # single value attributes with numpy variables
            s.write_attribute("attrStr", "Testing single string attribute")
            print(
                f"---- type of np.array(data.i8[0]) is {type(np.array(data.i8[0]))}"
                f" shape = {np.array(data.i8[0]).shape}"
            )
            s.write_attribute("attrI8", np.array(data.i8[0]))
            s.write_attribute("attrI16", np.array(data.i16[0]))
            s.write_attribute("attrI32", np.array(data.i32[0]))
            s.write_attribute("attrI64", np.array(data.i64[0]))
            s.write_attribute("attrU8", np.array(data.u8[0]))
            s.write_attribute("attrU16", np.array(data.u16[0]))
            s.write_attribute("attrU32", np.array(data.u32[0]))
            s.write_attribute("attrU64", np.array(data.u64[0]))
            s.write_attribute("attrR32", np.array(data.r32[0]))
            s.write_attribute("attrR64", np.array(data.r64[0]))
            s.write_attribute("attrC64", np.array(data.c64[0]))

            # single value attributes with Python variables
            s.write_attribute("attrNx", data.Nx)
            s.write_attribute("attr_int_value", i)
            s.write_attribute("attr_float_value", f)
            s.write_attribute("attr_complex_value", c)

            # array attributes with Python lists
            s.write_attribute("attr_int_list", data.int_list)
            s.write_attribute("attr_float_list", data.float_list)
            s.write_attribute("attr_complex_list", data.complex_list)
            s.write_attribute("attrStrArray1", ["string1"])
            s.write_attribute("attrStrArray3", ["string1", "string2", "string3"])

            # array attributes with numpy arrays
            s.write_attribute("attrI8Array", data.i8)
            s.write_attribute("attrI16Array", data.i16)
            s.write_attribute("attrI32Array", data.i32)
            s.write_attribute("attrI64Array", data.i64)
            s.write_attribute("attrU8Array", data.u8)
            s.write_attribute("attrU16Array", data.u16)
            s.write_attribute("attrU32Array", data.u32)
            s.write_attribute("attrU64Array", data.u64)
            s.write_attribute("attrR32Array", data.r32)
            s.write_attribute("attrR64Array", data.r64)
            s.write_attribute("attrC64Array", data.c64)

        s.write("steps", "Step:" + str(step.current_step()))
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
        s.write("varC64", data.c64, shape, start, count)

        if rank == 0 and step.current_step() == 0:
            # attribute assigned to variable
            s.write_attribute("size", data.Nx, "varI8")
            s.write_attribute("size", data.Nx, "varI16", "::")
            s.write_attribute("varattrStrArray", ["varattr1", "varattr2", "varattr3"], "steps")
            s.write_attribute("varattrI8Array", data.i8, "varI8")
            s.write_attribute("varattrI16Array", data.i16, "varI16")
            s.write_attribute("varattrI32Array", data.i32, "varI32")
            s.write_attribute("varattrI64Array", data.i64, "varI64")
            s.write_attribute("varattrU8Array", data.u8, "varU8")
            s.write_attribute("varattrU16Array", data.u16, "varU16")
            s.write_attribute("varattrU32Array", data.u32, "varU32")
            s.write_attribute("varattrU64Array", data.u64, "varU64")
            s.write_attribute("varattrR32Array", data.r32, "varR32")
            s.write_attribute("varattrR64Array", data.r64, "varR64")
            s.write_attribute("varattrR64Value", data.r64, "varR64")

comm.Barrier()

# Reader
data = SmallTestData(rank)

with Stream("types_np.bp", "r", comm=comm) as fr:
    # file only
    assert fr.num_steps() == 5

    for fr_step in fr.steps():
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
            print(f"tag = {inTag}")
            nx = fr_step.read("nx")
            print(f"nx = {nx}")
            assert nx == data.Nx
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
            print(f"attrStr = {inTag}")
            inNx = fr_step.read_attribute("attrNx")
            inI8 = fr_step.read_attribute("attrI8")
            inI16 = fr_step.read_attribute("attrI16")
            print(f"attrI16 = {inI16}")
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

            if inNx != data.Nx:
                raise ValueError("attrNx read failed")

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

            in_an_int_value = fr_step.read("an_int_value")
            print(f"an_int_value = {in_an_int_value} of type {type(in_an_int_value)}")
            assert in_an_int_value == data.int_list[0]

            in_a_float_value = fr_step.read("a_float_value")
            print(f"a_float_value = {in_a_float_value} of type {type(in_a_float_value)}")
            assert in_a_float_value == data.float_list[0]

            inC64 = fr_step.read("gvarC64")
            print(f"inC64 = {inC64} of type {type(inC64)}")
            assert inC64 == data.c64[0]

            in_a_complex_value = fr_step.read("a_complex_value")
            print(f"a_complex_value = {in_a_complex_value} of type {type(in_a_complex_value)}")
            assert in_a_complex_value == data.complex_list[0]

            an_int_list = fr_step.read("an_int_list")
            if not (an_int_list == data.int_list).all():
                raise ValueError("an_int_list array read failed")
            print(f"an_int_list = {an_int_list} of type {type(an_int_list)}")

            a_float_list = fr_step.read("a_float_list")
            if not (a_float_list == data.float_list).all():
                raise ValueError("a_float_list array read failed")
            print(f"a_float_list = {a_float_list} of type {type(a_float_list)}")

            a_complex_list = fr_step.read("a_complex_list")
            if not (a_complex_list == data.complex_list).all():
                raise ValueError("a_complex_list array read failed")
            print(f"a_complex_list = {a_complex_list} of type {type(a_complex_list)}")

            # Array attribute
            inStr1 = fr_step.read_attribute("attrStrArray1")
            print(f"attrStrArray1 = {inStr1}")
            inStr3 = fr_step.read_attribute("attrStrArray3")
            print(f"attrStrArray3 = {inStr3}")
            inI8 = fr_step.read_attribute("attrI8Array")
            inI16 = fr_step.read_attribute("attrI16Array")
            print(f"attrI16Array = {inI16}")
            inI32 = fr_step.read_attribute("attrI32Array")
            inI64 = fr_step.read_attribute("attrI64Array")
            inU8 = fr_step.read_attribute("attrU8Array")
            inU16 = fr_step.read_attribute("attrU16Array")
            inU32 = fr_step.read_attribute("attrU32Array")
            inU64 = fr_step.read_attribute("attrU64Array")
            inR32 = fr_step.read_attribute("attrR32Array")
            inR64 = fr_step.read_attribute("attrR64Array")

            if inStr1 != ["string1"]:
                raise ValueError("attrStrArray1 read failed")

            if inStr3 != ["string1", "string2", "string3"]:
                raise ValueError("attrStrArray3 read failed")

            if not (inI8 == data.i8).all():
                raise ValueError("attrI8 array read failed")

            if not (inI16 == data.i16).all():
                raise ValueError("attrI16 array read failed")

            if not (inI32 == data.i32).all():
                raise ValueError("attrI32 array read failed")

            if not (inI64 == data.i64).all():
                raise ValueError("attrI64 array read failed")

            if not (inU8 == data.u8).all():
                raise ValueError("attrU8 array read failed")

            if not (inU16 == data.u16).all():
                raise ValueError("attrU16 array read failed")

            if not (inU32 == data.u32).all():
                raise ValueError("attrU32 array read failed")

            if not (inU64 == data.u64).all():
                raise ValueError("attrU64 array read failed")

            if not (inR32 == data.r32).all():
                raise ValueError("attrR32 array read failed")

            if not (inR64 == data.r64).all():
                raise ValueError("attrR64 array read failed")

            # Array attributes (written as List)
            in_attr_float_list = fr_step.read_attribute("attr_float_list")
            in_attr_int_list = fr_step.read_attribute("attr_int_list")

            print(f"attr_float_list = {in_attr_float_list}")
            print(f"attr_int_list = {in_attr_int_list}")

            if not (in_attr_float_list == np.array(data.float_list)).all():
                raise ValueError("attr_float_list array read failed")

            if not (in_attr_int_list == np.array(data.int_list)).all():
                raise ValueError("attr_int_list array read failed")

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

            if not (inI8 == data.i8).all():
                raise ValueError("var attrI8 array read failed")

            if not (inI16 == data.i16).all():
                raise ValueError("var attrI16 array read failed")

            if not (inI32 == data.i32).all():
                raise ValueError("var attrI32 array read failed")

            if not (inI64 == data.i64).all():
                raise ValueError("var attrI64 array read failed")

            if not (inU8 == data.u8).all():
                raise ValueError("var attrU8 array read failed")

            if not (inU16 == data.u16).all():
                raise ValueError("var attrU16 array read failed")

            if not (inU32 == data.u32).all():
                raise ValueError("var attrU32 array read failed")

            if not (inU64 == data.u64).all():
                raise ValueError("var attrU64 array read failed")

            if not (inR32 == data.r32).all():
                raise ValueError("var attrR32 array read failed")

            if not (inR64 == data.r64).all():
                raise ValueError("var attrR64 array read failed")

            # Attributes assigned to a variable
            sizeI8 = fr_step.read_attribute("size", "varI8")
            sizeI16 = fr_step.read_attribute("size", "varI16", "::")

            if sizeI8 != data.Nx:
                raise ValueError("attribute varI8/size read failed")

            if sizeI16 != data.Nx:
                raise ValueError("attribute varI16::size read failed")

            sizeI8 = fr_step.read_attribute("varI8/size")
            sizeI16 = fr_step.read_attribute("varI16::size")

            if sizeI8 != data.Nx:
                raise ValueError("attribute varI8/size read failed")

            if sizeI16 != data.Nx:
                raise ValueError("attribute varI16::size read failed")

            step_attrs = fr_step.available_attributes()
            # for name, info in step_attrs.items():
            #     print(f"attribute {name} : {info}")

            step_attrs = fr_step.available_attributes("varI8")
            # for name, info in step_attrs.items():
            #    print(f"attribute {name} : {info}")
            if not [*step_attrs] == ["size", "varattrI8Array"]:
                raise ValueError("getting attributes of varI8 failed")

            step_attrs = fr_step.available_attributes("varI16", "::")
            for name, info in step_attrs.items():
                print(f"attribute {name} : {info}")
            if not [*step_attrs] == ["size"]:
                raise ValueError("getting attributes of varI16 with separator :: failed")

        stepStr = "Step:" + str(step)

        instepStr = fr_step.read("steps")
        if instepStr != stepStr:
            raise ValueError("steps variable read failed: " + instepStr + " " + stepStr)

        indataRanks = fr_step.read("rank", [0], [size])
        dataRanks = np.arange(0, size)
        if not (indataRanks == dataRanks).all():
            raise ValueError("Ranks read failed")

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

        if not (indataI8 == data.i8).all():
            raise ValueError("i8 array read failed")

        if not (indataI16 == data.i16).all():
            raise ValueError("i16 array read failed")

        if not (indataI32 == data.i32).all():
            raise ValueError("i32 array read failed")

        if not (indataI64 == data.i64).all():
            raise ValueError("i64 array read failed")

        if not (indataU8 == data.u8).all():
            raise ValueError("u8 array read failed")

        if not (indataU16 == data.u16).all():
            raise ValueError("u16 array read failed")

        if not (indataU32 == data.u32).all():
            raise ValueError("u32 array read failed")

        if not (indataU64 == data.u64).all():
            raise ValueError("u64 array read failed")

        if not (indataR32 == data.r32).all():
            raise ValueError("r32 array read failed")

        if not (indataR64 == data.r64).all():
            raise ValueError("r64 array read failed")

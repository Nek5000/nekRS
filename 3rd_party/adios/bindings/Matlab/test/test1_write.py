#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# Test writer to produce a test dataset for Matlab reading
# This is serial since the Matlab reader is built with serial ADIOS
#
#
#  Created on: Jan 28, 2019
#      Author: Norbert Podhorszki, pnorbert@ornl.gov
#

import numpy
import adios2

# User data
NRows = 5
NCols = 6

shape = [NRows, NCols]
start = [0, 0]
count = [NRows, NCols]

temperatures = numpy.zeros(NRows * NCols, dtype=numpy.int16)
npcols = numpy.array([NCols])
nprows = numpy.array([NRows])

value = (NRows * NCols) + 1
for i in range(0, NRows):
    for j in range(0, NCols):
        temperatures[i * NCols + j] = value
        value = value + 1

# ADIOS2 low-level API for Write
adios = adios2.ADIOS()
io = adios.DeclareIO("writer")

varStr = io.DefineVariable("note")
varCols = io.DefineVariable("ncols", npcols)
varRows = io.DefineVariable("nrows", nprows)
varT = io.DefineVariable(
    "temperature2D", temperatures, shape, start, count, adios2.ConstantDims)

io.DefineAttribute("aaa", numpy.array([3.1415]))
io.DefineAttribute("anote", "just a string")
io.DefineAttribute("adimnames", ["rows", "columns"])
npdims = numpy.array([NRows, NCols], dtype=numpy.int32)
io.DefineAttribute("adims", npdims)


fw = io.Open("test1.bp", adios2.Mode.Write)
fw.Put(varStr, "This is an ADIOS2 output")
fw.Put(varCols, npcols)
fw.Put(varRows, nprows)
fw.Put(varT, temperatures)
fw.Close()

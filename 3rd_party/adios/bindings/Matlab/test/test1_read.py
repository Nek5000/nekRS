#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# Test reader in Python that should have the same output
# as the test1.m Matlab reader
#
#
#  Created on: Jan 28, 2019
#      Author: Norbert Podhorszki, pnorbert@ornl.gov
#

import numpy
import adios2

adios = adios2.ADIOS()
io = adios.DeclareIO("reader")


# ADIOS2 high-level API for Reading
fr = io.Open("test1.bp", adios2.Mode.Read)


anote = io.InquireAttribute("anote")
adimnames = io.InquireAttribute("adimnames")
adims = io.InquireAttribute("adims")

dimNames = adimnames.DataString()
dims = adims.Data()

print("Info based on attributes in file:")
print("  Dimensions  {0}  {1}".format(dimNames[0], dimNames[1]))
print("               {0}     {1}".format(dims[0], dims[1]))

nrows = numpy.zeros(1, dtype=numpy.int64)
varRows = io.InquireVariable("nrows")
fr.Get(varRows, nrows)

ncols = numpy.zeros(1, dtype=numpy.int64)
varCols = io.InquireVariable("ncols")
fr.Get(varCols, ncols)

note = "_______________________________________"
varNote = io.InquireVariable("note")
fr.Get(varNote, note, adios2.Mode.Sync)

fr.PerformGets()

print("# of rows = {0}".format(nrows[0]))
print("# of cols = {0}".format(ncols[0]))
print("Note = {0}".format(note))
#
# inTemperatures = fr.read("temperature2D")
# print("temperature2d array size = " + str(inTemperatures.size))
#
# for row in inTemperatures:
#     print(''.join(['{:7}'.format(item) for item in row]))
#
fr.Close()

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

value = (NRows * NCols) + 1
for i in range(0, NRows):
    for j in range(0, NCols):
        temperatures[i * NCols + j] = value
        value = value + 1

# ADIOS2 high-level API for Write
fw = adios2.open("test1.bp", "w")
fw.write("note", 'This is an ADIOS2 output')
fw.write("temperature2D", temperatures, shape, start, count)
fw.write("nrows", numpy.array([NRows]))
fw.write("ncols", numpy.array([NCols]))
fw.close()

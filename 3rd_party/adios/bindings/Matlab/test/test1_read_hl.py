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

import adios2

# ADIOS2 high-level API for Reading
fr = adios2.open("test1.bp", "r")

inNRows = fr.read("nrows")
print("# of rows = {0}".format(inNRows[0]))

inNCols = fr.read("ncols")
print("# of cols = {0}".format(inNCols[0]))

inNote = fr.readstring("note")
print("Note = {0}".format(inNote))

inTemperatures = fr.read("temperature2D")
print("temperature2d array size = " + str(inTemperatures.size))

for row in inTemperatures:
    print(''.join(['{:7}'.format(item) for item in row]))

fr.close()

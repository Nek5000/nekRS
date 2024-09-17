#!/bin/bash

# Check if output files exists
if [ ! -d balance_atom.bp ]
then
   echo "ERROR: balance_atom.bp was not produced"
   exit 1
fi

if [ ! -d balance_custom.bp ]
then
   echo "ERROR: balance_custom.bp was not produced"
   exit 2
fi

if [ ! -f balance.dump ]
then
   echo "ERROR: balance.dump was not produced"
   exit 3
fi

# Compare data from default I/O balance.dump and ADIOS output

BOUNDARY_DEFAULT=$(grep BOX balance.dump  | head -1 | cut -f 4,5,6 -d " ")
BOUNDARY_ADIOS=$(bpls -l balance_custom.bp -a boundarystr | cut -d '"' -f 2)

if [ ! "$BOUNDARY_DEFAULT" = "$BOUNDARY_ADIOS" ]
then
    echo "balance.dump:[$BOUNDARY_DEFAULT] does not match balance_custom.bp/boundarystr:[$BOUNDARY_ADIOS]"
    exit 4
fi

echo "balance.dump:[$BOUNDARY_DEFAULT] matches balance_custom.bp/boundarystr:[$BOUNDARY_ADIOS]"
exit 0

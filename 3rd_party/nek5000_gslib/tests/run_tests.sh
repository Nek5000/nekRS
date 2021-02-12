#!/bin/sh

# Set the MPI command and number of MPI ranks
MPI=mpirun
np="1 2 3 4"

# Build the tests if they have not been built
# before.
make -j -C .. CC=mpicc CFLAGS="-g" NBC=1 tests

for i in *.o fortran/*.o; do
  j=${i%.*}
  for n in $np; do
    $MPI -np $n ./$j >> test_log
    if [ "$?" -eq 0 ]; then
      echo "Running test: $j, np: $n ... Passed."
    else
      echo "Running test: $j, np: $n ... Failed."
    fi
  done
done

#!/bin/bash
: ${NEKRS_HOME:="/home/sean/.local/nekrs_next_ascent/bin/../"}
: ${ASCENT_DIR="/home/sean/ascent_repo021424/install/ascent-develop/"}
: ${COMPILE_METHOD=3}

echo "COMPILE_METHOD=$COMPILE_METHOD"


# method 1: ascent_config.mk
if [ "$COMPILE_METHOD" -eq "1" ]; then
s="helper.mak"
echo "include $ASCENT_DIR/share/ascent/ascent_config.mk" > $s
echo ".PHONY: printvars" >> $s
echo "printvars:" >> $s
echo -e "\t@echo \"ASCENT_INCLUDE_FLAGS=\\\"\$(ASCENT_INCLUDE_FLAGS)\\\"\"" >> $s
echo -e "\t@echo \"ASCENT_LINK_RPATH=\\\"\$(ASCENT_LINK_RPATH)\\\"\"" >> $s
echo -e "\t@echo \"ASCENT_MPI_LIB_FLAGS=\\\"\$(ASCENT_MPI_LIB_FLAGS)\\\"\"" >> $s
make -f helper.mak printvars > make.tmp 2>/dev/null
source make.tmp

NEKRS_UDF_LDFLAGS="${ASCENT_LINK_RPATH} ${ASCENT_MPI_LIB_FLAGS}"
NEKRS_UDF_LDFLAGS=`echo $NEKRS_UDF_LDFLAGS|awk '{$1=$1;print}'` # rm tailing white space
NEKRS_UDF_INCLUDES="$ASCENT_INCLUDE_FLAGS"
NEKRS_UDF_INCLUDES=`echo $NEKRS_UDF_INCLUDES|sed 's/-I//g'|awk '{$1=$1;print}'|tr -s " "|sed 's/ /\;/g'`
export NEKRS_UDF_INCLUDES="$NEKRS_UDF_INCLUDES"
export NEKRS_UDF_LDFLAGS="$NEKRS_UDF_LDFLAGS"


# method 2: cmake (udf/CMakeList.txt)
elif [ "$COMPILE_METHOD" -eq "2" ]; then
export NEKRS_ENABLE_ASCENT=1
export NEKRS_ASCENT_INSTALL_DIR=$ASCENT_DIR


# method 3
elif [ "$COMPILE_METHOD" -eq "3" ]; then
export NEKRS_UDF_LDFLAGS="-L${ASCENT_DIR}/lib -lascent_mpi"
export NEKRS_UDF_INCLUDES="${ASCENT_DIR}/../conduit-v0.8.8/include/conduit\;${ASCENT_DIR}/../conduit-v0.8.8/include\;${ASCENT_DIR}/include/ascent"
fi


# FIXME: this is a temporary hack on my laptop
export LD_LIBRARY_PATH+=${ASCENT_DIR}/../ascent-develop/lib:
export LD_LIBRARY_PATH+=${ASCENT_DIR}/../vtk-m-v2.1.0/lib: 


## TODO: call nrspre directly
export NEKRS_HOME=$NEKRS_HOME

if [ $# -eq 0 ] || [ $# -ne 2 ] || [ "$1" == "-h" ] || [ "$1" == "-help" ]; then
  echo "usage: ${0##*/} <casename> <#target procs>"
  exit 0
fi
if [ "$1" == "-clean" ]; then
  rm -r .cache/udf
fi

mpirun -np 1 $NEKRS_HOME/bin/nekrs --setup $1 --build-only $2



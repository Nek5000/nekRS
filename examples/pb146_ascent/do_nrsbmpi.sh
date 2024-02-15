#!/bin/bash
: ${NEKRS_HOME:="/home/sean/.local/nekrs_next_ascent/bin/../"}
: ${NEKRS_HOME:="/home/sean/.local/nekrs_next_ascent/bin/../"}
: ${ASCENT_DIR="/home/sean/ascent/install/ascent-develop/"}
#: ${ASCENT_DIR="/home/sean/ascent_repo021424/install/ascent-develop/"}

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
export NEKRS_HOME=$NEKRS_HOME

# FIXME: this is a temporary hack
export LD_LIBRARY_PATH+=/home/sean/ascent/install/vtk-m-v2.0.0/lib
#export LD_LIBRARY_PATH+=/home/sean/ascent_repo021424/install/vtk-m-v2.1.0/lib

## TODO: call nrspre directly
rm -f logfile
mv $1.log.$2 $1.log1.$2 2>/dev/null

ulimit -s unlimited 2>/dev/null

if [ $# -eq 0 ] || [ $# -lt 2 ] || [ "$1" == "-h" ] || [ "$1" == "-help" ]; then
  echo "usage: ${0##*/} casename #tasks [args]"
  exit 1
fi


nohup mpirun --mca osc ucx -np $2 $NEKRS_HOME/bin/nekrs --setup $1 ${@:3} >$1.log.$2 </dev/null 2>&1 &

ln -sf $1.log.$2 logfile
echo "started job in background, redirecting output to ./logfile ..."

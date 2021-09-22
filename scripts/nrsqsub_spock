#!/bin/bash

: ${PROJ_ID:=""}
: ${QUEUE:=""}
: ${NEKRS_HOME:="$HOME/.local/nekrs"}

export NVME_HOME="/mnt/bb/$USER/"

if [ $# -ne 3 ]; then
  echo "usage: [PROJ_ID] [QUEUE] $0 <casename> <number of compute nodes> <hh:mm:ss>"
  exit 0
fi

if [ -z "$PROJ_ID" ]; then
  echo "ERROR: PROJ_ID is empty"
  exit 1
fi

if [ -z "$QUEUE" ]; then
  echo "ERROR: QUEUE is empty"
  exit 1
fi

bin=${NEKRS_HOME}/bin/nekrs
case=$1
nodes=$2
gpu_per_node=4
cores_per_socket=16
let nn=$nodes*$gpu_per_node
let ntasks=nn
time=$3
backend=HIP


if [ ! -f $bin ]; then
  echo "Cannot find" $bin
  exit 1
fi

if [ ! -f $case.par ]; then
  echo "Cannot find" $case.par
  exit 1
fi

if [ ! -f $case.udf ]; then
  echo "Cannot find" $case.udf
  exit 1
fi

if [ ! -f $case.oudf ]; then
  echo "Cannot find" $case.oudf
  exit 1
fi

if [ ! -f $case.re2 ]; then
  echo "Cannot find" $case.re2
  exit 1
fi


# romio setup
export ROMIO_HINTS="$(pwd)/.romio_hint"
if [ ! -f "$ROMIO_HINTS" ]; then
  echo "romio_no_indep_rw true"   >$ROMIO_HINTS
  echo "romio_cb_write enable"   >>$ROMIO_HINTS
  echo "romio_ds_write enable"   >>$ROMIO_HINTS
  echo "romio_cb_read enable"    >>$ROMIO_HINTS
  echo "romio_ds_read enable"    >>$ROMIO_HINTS
  echo "cb_buffer_size 16777216" >>$ROMIO_HINTS
  echo "cb_config_list *:1"      >>$ROMIO_HINTS
fi


# sbatch
SFILE=s.bin
echo "#!/bin/bash" > $SFILE
echo "#SBATCH -A $PROJ_ID" >>$SFILE
echo "#SBATCH -J nekRS_$case" >>$SFILE
echo "#SBATCH -o %x-%j.out" >>$SFILE
echo "#SBATCH -t $time" >>$SFILE
echo "#SBATCH -N $nodes" >>$SFILE
echo "#SBATCH -p $QUEUE" >>$SFILE
echo "#SBATCH -C nvme" >>$SFILE
echo "#SBATCH --exclusive" >>$SFILE
echo "#SBATCH --ntasks-per-node=$gpu_per_node" >>$SFILE
echo "#SBATCH --gpus-per-task=1" >>$SFILE
echo "#SBATCH --gpu-bind=closest" >>$SFILE
echo "#SBATCH --cpus-per-task=$cores_per_socket" >>$SFILE

echo "module load craype-accel-amd-gfx908" >>$SFILE
echo "module load rocm" >>$SFILE
echo "module load cmake" >>$SFILE
echo "module load PrgEnv-gnu" >>$SFILE
echo "module load wget" >>$SFILE
echo "module unload cray-libsci" >>$SFILE
echo "module list" >>$SFILE
echo "rocm-smi" >>$SFILE
echo "rocm-smi --showpids" >>$SFILE

echo "export MPICH_GPU_SUPPORT_ENABLED=1" >>$SFILE

echo "## OLCFDEV-118: Parallel HDF5 failures when running on GPFS" >>$SFILE
echo "export ROMIO_FSTYPE_FORCE=\"ufs:\"" >>$SFILE

echo "ulimit -s unlimited " >>$SFILE
echo "export NEKRS_HOME=$NEKRS_HOME" >>$SFILE
echo "export NEKRS_GPU_MPI=1 " >>$SFILE

echo "export NVME_HOME=$NVME_HOME" >>$SFILE
echo "export ROMIO_HINTS=$ROMIO_HINTS" >>$SFILE

#echo "# precompile" >>$SFILE
#echo "date" >>$SFILE
#echo "srun -n 1 $bin --setup $1 --backend $backend --build-only $ntasks" >>$SFILE
#echo "srun -n $nodes cp -a .cache $NVME_HOME" >>$SFILE
#echo "export NEKRS_CACHE_DIR=$NVME_HOME/.cache" >>$SFILE

echo "# actual run" >>$SFILE
echo "date" >>$SFILE
echo "srun -n $ntasks $bin --backend $backend --device-id 0 --setup $case" >>$SFILE

sbatch $SFILE


# clean-up
#rm -rf $SFILE $ROMIO_HINTS

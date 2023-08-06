#!/bin/bash

: ${PROJ_ID:=""}
: ${QUEUE:="batch"}
: ${NEKRS_HOME:="$HOME/.local/nekrs"}
: ${NEKRS_CACHE_BCAST:=1}
: ${NEKRS_SKIP_BUILD_ONLY:=0}

export NVME_HOME="/mnt/bb/$USER/"

if [ $# -ne 3 ]; then
  echo "usage: [PROJ_ID=value] [QUEUE=value] $0 <casename> <number of compute nodes> <hh:mm:ss>"
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
gpu_per_node=8
cores_per_numa=7
let nn=$nodes*$gpu_per_node
let ntasks=nn
time=$3
backend=HIP

time_fmt=`echo $time|tr ":" " "|awk '{print NF}'`
if [ "$time_fmt" -ne "3" ]; then
  echo "Warning: time is not in the format <hh:mm:ss>"
  echo $time
  exit 1
fi

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

if [ ! -f $case.re2 ]; then
  echo "Cannot find" $case.re2
  exit 1
fi

striping_unit=16777216
max_striping_factor=400
let striping_factor=$nodes/2
if [ $striping_factor -gt $max_striping_factor ]; then
  striping_factor=$max_striping_factor
fi
if [ $striping_factor -lt 1 ]; then
  striping_factor=1
fi


MPICH_MPIIO_HINTS="*:cray_cb_write_lock_mode=2:cray_cb_nodes_multiplier=4:striping_unit=${striping_unit}:striping_factor=${striping_factor}:romio_cb_write=enable:romio_ds_write=disable:romio_no_indep_rw=true"

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
echo "#SBATCH --cpus-per-task=$cores_per_numa" >>$SFILE

echo "module load PrgEnv-gnu" >> $SFILE
echo "module load craype-accel-amd-gfx90a" >> $SFILE
echo "module load cray-mpich" >> $SFILE
echo "module load rocm" >> $SFILE
echo "module unload cray-libsci" >> $SFILE
echo "module list" >> $SFILE

echo "rocm-smi" >>$SFILE
echo "rocm-smi --showpids" >>$SFILE

echo "squeue -u \$USER" >>$SFILE

echo "export MPICH_GPU_SUPPORT_ENABLED=1" >>$SFILE

## These must be set before compiling so the executable picks up GTL
echo "export PE_MPICH_GTL_DIR_amd_gfx90a=\"-L${CRAY_MPICH_ROOTDIR}/gtl/lib\"" >> $SFILE
echo "export PE_MPICH_GTL_LIBS_amd_gfx90a=\"-lmpi_gtl_hsa\"" >> $SFILE

echo "ulimit -s unlimited " >>$SFILE
echo "export NEKRS_HOME=$NEKRS_HOME" >>$SFILE
echo "export NEKRS_GPU_MPI=1 " >>$SFILE

echo "export NVME_HOME=$NVME_HOME" >>$SFILE

echo "export MPICH_MPIIO_HINTS=$MPICH_MPIIO_HINTS" >>$SFILE
echo "export MPICH_MPIIO_STATS=1" >>$SFILE

echo "export MPICH_OFI_NIC_POLICY=NUMA" >>$SFILE
echo "export NEKRS_CACHE_BCAST=$NEKRS_CACHE_BCAST" >> $SFILE

echo "if [ \$NEKRS_CACHE_BCAST -eq 1 ]; then" >> $SFILE
echo "  export NEKRS_LOCAL_TMP_DIR=\$NVME_HOME" >> $SFILE
echo "fi" >> $SFILE
echo "" >> $SFILE
echo "date" >>$SFILE
echo "" >> $SFILE

bin_nvme=$NVME_HOME"nekrs-bin"
bin_nvme_libs=$bin_nvme"_libs"
echo "sbcast -fp --send-libs $bin $bin_nvme" >> $SFILE
echo "if [ ! \"\$?\" == \"0\" ]; then"  >> $SFILE
echo "    echo \"SBCAST failed!\"" >> $SFILE
echo "    exit 1" >> $SFILE
echo "fi" >> $SFILE

echo "export LD_LIBRARY_PATH=$bin_nvme_libs:${LD_LIBRARY_PATH}" >> $SFILE
echo "export LD_PRELOAD=$bin_nvme_libs/libnekrs.so:$bin_nvme_libs/libocca.so:$bin_nvme_libs/libnekrs-hypre-device.so:$bin_nvme_libs/libnekrs-hypre.so" >> $SFILE

echo "ls -ltra $NVME_HOME" >> $SFILE
echo "ls -ltra $bin_nvme_libs" >> $SFILE
echo "ldd $bin_nvme" >> $SFILE

if [ $NEKRS_SKIP_BUILD_ONLY -eq 0 ]; then
echo "# precompilation" >>$SFILE
echo "srun -N 1 -n 1 $bin_nvme --backend $backend --device-id 0 --setup $case --build-only $ntasks" >>$SFILE
fi
echo "" >> $SFILE

echo "# actual run" >>$SFILE
echo "srun -N $nodes -n $ntasks $bin_nvme --backend $backend --device-id 0 --setup $case" >>$SFILE

sbatch $SFILE

# clean-up

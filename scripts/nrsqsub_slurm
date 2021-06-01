#!/bin/bash -l
#
#SBATCH --job-name="nekRS"
##SBATCH --constraint=gpu
#SBATCH --exclusive
#SBATCH --gres=gpu:a100-sxm4-40gb:8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=32

module load cmake
module load gcc
module load cuda
module load openmpi

ulimit -s unlimited 
export NEKRS_HOME=${HOME}/.local/nekrs
export OCCA_CXX=g++
export OCCA_CXXFLAGS="-O2 -ftree-vectorize -funroll-loops -march=native -mtune=native"

export OMPI_MCA_pml=ucx
export OMPI_MCA_btl="^vader,tcp,openib,smcuda"
export OMPI_MCA_osc=ucx
export OMPI_MCA_io=romio321

export UCX_TLS=rc,sm,cuda
export UCX_RNDV_SCHEME=put_zcopy
export UCX_RNDV_THRESH=1024

export NEKRS_GPU_MPI=1 

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

# precompile
date
mpirun -np 1 ${NEKRS_HOME}/bin/nekrs --setup $1 --backend CUDA --build-only $SLURM_NTASKS

# actual run
date
CMD=.lhelper
echo "#!/bin/bash" >$CMD
echo "GPUS=(0 1 2 3 4 5 6 7)" >>$CMD
echo "CPUS=(3 3 1 1 7 7 5 5)" >>$CMD 
echo "NICS=(mlx5_0 mlx5_1 mlx5_2 mlx5_3 mlx5_6 mlx5_7 mlx5_8 mlx5_9)" >>$CMD
echo "export UCX_NET_DEVICES=\${NICS[\${OMPI_COMM_WORLD_LOCAL_RANK}]}:1" >>$CMD
echo "export CUDA_VISIBLE_DEVICES=\${GPUS[\${OMPI_COMM_WORLD_LOCAL_RANK}]}" >>$CMD
echo "#echo LOCAL_RANK=\${OMPI_COMM_WORLD_LOCAL_RANK} GPU=\$CUDA_VISIBLE_DEVICES NIC=\$UCX_NET_DEVICES CPU=\${CPUS[\${OMPI_COMM_WORLD_LOCAL_RANK}]}" >>$CMD
echo "numactl --cpunodebind=\${CPUS[\${OMPI_COMM_WORLD_LOCAL_RANK}]} --membind=\${CPUS[\${OMPI_COMM_WORLD_LOCAL_RANK}]} ${NEKRS_HOME}/bin/nekrs --backend CUDA --device-id 0 --setup $1" >>$CMD
chmod 755 $CMD

mpirun -np $SLURM_NTASKS ./$CMD

# clean-up
rm -rf $CMD $ROMIO_HINTS

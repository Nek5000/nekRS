#!/bin/bash   

#Nvidia DGX A100
GPUS=(0 1 2 3 4 5 6 7)                                             
CPUS=(3 3 1 1 7 7 5 5)  
NICS=(mlx5_0 mlx5_1 mlx5_2 mlx5_3 mlx5_6 mlx5_7 mlx5_8 mlx5_9)

#do not edit below here
APP="$*"

if [ ! -z ${OMPI_COMM_WORLD_RANK} ]; then
    RANK=${OMPI_COMM_WORLD_RANK}
    LOCAL_RANK=${OMPI_COMM_WORLD_LOCAL_RANK}
elif [ ! -z ${MV2_COMM_WORLD_RANK} ]; then
    RANK=${MV2_COMM_WORLD_RANK}
    LOCAL_RANK=${MV2_COMM_WORLD_LOCAL_RANK}
else    
    RANK=${SLURM_PROCID}
    LOCAL_RANK=${SLURM_LOCALID}
fi 
GPU=${GPUS[$LOCAL_RANK]}
NIC=${NICS[$LOCAL_RANK]}
CPU=${CPUS[$LOCAL_RANK]}

export CUDA_VISIBLE_DEVICES=$GPU
#export HIP_VISIBLE_DEVICES=$GPU

export UCX_NET_DEVICES=$NIC:1
export UCX_TLS=rc,sm,cuda
export UCX_RNDV_SCHEME=put_zcopy
export UCX_RNDV_THRESH=1024
export UCX_MEMTYPE_CACHE=n

export OMPI_MCA_pml=ucx
export OMPI_MCA_btl="^vader,tcp,openib,smcuda"
export OMPI_MCA_osc=ucx

export NEKRS_GPU_MPI=0
if [[ $UCX_TLS == *"cuda"* ]]; then
  export NEKRS_GPU_MPI=1
fi

ulimit -s unlimited 2>/dev/null

COMMAND="numactl --cpunodebind=$CPU --membind=$CPU $APP --device-id 0"
echo "RANK=$RANK, LOCAL_RANK=$LOCAL_RANK, GPU=$GPU, NIC=$NIC, CPU=$CPU, NEKRS_GPU_MPI=$NEKRS_GPU_MPI"

$COMMAND

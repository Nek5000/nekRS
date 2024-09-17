#!/bin/bash
#SBATCH -A CSC383_crusher
#SBATCH -J gs-julia-1MPI-1GPU
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH -t 0:02:00
#SBATCH -p batch
#SBATCH -N 1

date

GS_DIR=/gpfs/alpine/proj-shared/csc383/wgodoy/ADIOS2/examples/simulations/GrayScott.jl
GS_EXE=$GS_DIR/gray-scott.jl

srun -n 1 --gpus=1 julia --project=$GS_DIR $GS_EXE settings-files.json

# launch this file with sbatch `$ sbatch job_crusher.sh`

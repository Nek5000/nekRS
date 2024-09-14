#!/bin/bash
#    Begin LSF directives
#BSUB -P csc383
#BSUB -W 00:02
#BSUB -nnodes 1
#BSUB -J gs-julia
#BSUB -o output.%J
#BSUB -e output.%J 
#BSUB -N godoywf@ornl.gov
#    End BSUB directives and begin shell commands

date
GS_DIR=/gpfs/alpine/proj-shared/csc383/wgodoy/ADIOS2/examples/simulations/GrayScott.jl
GS_EXE=$GS_DIR/gray-scott.jl

jsrun -n 1 -g 1 julia --project=$GS_DIR $GS_EXE settings-files.json

# launch this file with bsub `$ bsub job_summit.sh`

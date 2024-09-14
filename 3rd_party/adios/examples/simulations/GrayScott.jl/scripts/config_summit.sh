#!/bin/bash

# Replace these 3 entries
PROJ_DIR=/gpfs/alpine/proj-shared/csc383
export JULIA_DEPOT_PATH=$PROJ_DIR/etc/summit/julia_depot
GS_DIR=$PROJ_DIR/wgodoy/ADIOS2/examples/simulations/GrayScott.jl

# remove existing generated Manifest.toml
rm -f $GS_DIR/Manifest.toml
rm -f $GS_DIR/LocalPreferences.toml

# good practice to avoid conflicts with existing default modules 
# needed to avoid seg fault with MPI
module purge

# load required modules
module load spectrum-mpi
module load gcc/12.1.0 # needed by julia libraries
module load cuda/11.0.3 # failure with 11.5.2
module load adios2/2.8.1
# module load julia/1.8.2 not working with CUDA.jl as it's missing libLLVM-13jl.so
export PATH=$PROJ_DIR/opt/summit/julia-1.9.0-beta3/bin:$PATH

# Required to enable underlying ADIOS2 library from loaded module
export JULIA_ADIOS2_PATH=$OLCF_ADIOS2_ROOT

# MPIPreferences to use spectrum-mpi
julia --project=$GS_DIR -e 'using Pkg; Pkg.add("MPIPreferences")'
julia --project=$GS_DIR -e 'using MPIPreferences; MPIPreferences.use_system_binary(; library_names=["libmpi_ibm"], mpiexec="jsrun")'

# Instantiate the project by installing packages in Project.toml
julia --project=$GS_DIR -e 'using Pkg; Pkg.instantiate()'

# Adds to LocalPreferences.toml to use underlying system CUDA since CUDA.jl v4.0.0
# https://cuda.juliagpu.org/stable/installation/overview/#Using-a-local-CUDA
julia --project=$GS_DIR -e 'using CUDA; CUDA.set_runtime_version!("local")'

# Adds a custom branch in case the development version is needed (for devs to test new features)
julia --project=$GS_DIR -e 'using Pkg; Pkg.add(url="https://github.com/eschnett/ADIOS2.jl.git", rev="main")'

# Build the new ADIOS2
julia --project=$GS_DIR -e 'using Pkg; Pkg.build()'
julia --project=$GS_DIR -e 'using Pkg; Pkg.precompile()'

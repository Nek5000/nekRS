#!/bin/bash
#
# This script is for configuring adios on the authors' machines
# You can study it to figure out how to configure adios on your system
#

# Usage: runconf  [ON|OFF]  [lean]
#    ON|OFF   build parallel version with MPI | build serial version
#    lean     build without external dependencies

USE_MPI=${1:-ON}
echo USE_MPI=$USE_MPI

BUILD_FULL=true
if [ x"$2" == x"lean" ]; then
    BUILD_FULL=false
    echo "Build ADIOS without extra packages"
fi

SRCDIR=`dirname ${BASH_SOURCE[0]}`

if [ `hostname | cut -c 1-5` == "titan" ]; then

    #########
    # Titan #
    #########
    TARGET=`module list 2>&1 | grep "PrgEnv"- | sed "s/^.*PrgEnv-\([a-z]*\).*/\1/"`
    if [ -z "$TARGET" ]; then
        echo "Cannot determine Programming environment. Exit"
        exit 1
    fi

    echo "Configure on Titan for $TARGET env. MPI = ${USE_MPI}  use external packages = ${BUILD_FULL}"

    module load cmake3/3.9.0
    module unload craype-interlagos craype-interlagos-cu 
    module unload craype-barcelona craype-abudhabi craype-abudhabi-cu 
    module unload craype-ivybridge craype-sandybridge craype-haswell 
    module unload craype-mc8 craype-mc12
    module load craype-istanbul
    unset EXTERNAL_OPTS
    if [ $BUILD_FULL == "true" ]; then
        # ADIOS 1.x
        #module load adios/1.13.0
        export ADIOS_DIR=/ccs/proj/e2e/pnorbert/ADIOS/xk6.$TARGET
        # HDF5
        if [ "x$USE_MPI" == "xON" ]; then
            module load cray-hdf5-parallel
        else
            module load cray-hdf5
        fi
        EXTERNAL_OPTS+="-DADIOS2_USE_HDF5=ON  -DHDF5_ROOT=${HDF5_DIR} "
        # ZFP
        #      -DZFP_ROOT=/opt/zfp 
        # EVPATH
        #      -DEVPath_ROOT=/opt/chaos 
        EXTERNAL_OPTS+="-DADIOS2_USE_Lustre=ON "
    else
        EXTERNAL_OPTS+="-DADIOS2_USE_HDF5=OFF -DADIOS2_USE_SST=OFF "
        EXTERNAL_OPTS+="-DADIOS2_USE_ZeroMQ=OFF -DADIOS2_USE_BZip2=OFF -DADIOS2_USE_ZFP=OFF -DADIOS2_USE_SZ=OFF"
    fi
    
    export CXXFLAGS="-fPIC -g -O2"
    export CFLAGS="-fPIC -g -O2"
    cmake -DCMAKE_INSTALL_PREFIX=${PWD}/install \
          -DADIOS2_USE_MPI=${USE_MPI} \
          -DADIOS2_USE_Fortran=ON \
          -DADIOS2_USE_Python=OFF \
          -DCMAKE_BUILD_TYPE=RelWithDebInfo \
          -DBUILD_TESTING=OFF \
          -DADIOS2_BUILD_EXAMPLES=ON \
          ${EXTERNAL_OPTS} \
          ${SRCDIR}

elif [ `hostname -f | cut -c 1-9` == "summitdev" ]; then
    #####################
    #  SUMMITDEV IBM    #
    #####################
    MLIST=`module list 2>&1 ` #| grep "PrgEnv"- | sed "s/^.*PrgEnv-\([a-z]*\).*/\1/"`
    TARGET=`echo $MLIST | grep " xl\/" | sed "s/^.* xl\/.*/xl/"`
    if [ -z $TARGET ]; then
        TARGET=`echo $MLIST | grep " pgi\/" | sed "s/^.* pgi\/.*/pgi/"`
        if [ -z $TARGET ]; then
            TARGET=`echo $MLIST | grep " gnu\/" | sed "s/^.* gnu\/.*/gnu/"`
            if [ -z $TARGET ]; then
                TARGET=`echo $MLIST | grep " clang\/" | sed "s/^.* clang\/.*/clang/"`
                if [ -z $TARGET ]; then
                    echo "Cannot determine Programming environment. Exit"
                    exit 1
                fi
            fi
        fi
    fi
    echo "Configure on Summitdev for $TARGET env. MPI = ${USE_MPI}  use external packages = ${BUILD_FULL}"

    module load cmake/3.9.2

    unset EXTERNAL_OPTS
    if [ $BUILD_FULL == "true" ]; then
        # ADIOS 1.x
        EXTERNAL_OPTS+="-DADIOS2_USE_HDF5=OFF "
        # ZFP
        #      -DZFP_ROOT=/opt/zfp 
        # EVPATH
        EXTERNAL_OPTS+="-DEVPath_ROOT=/ccs/proj/csc143/SST_XGC_coupling/summitdev/sw.gcc/chaos "
    else
        EXTERNAL_OPTS+="-DADIOS2_USE_HDF5=OFF"
    fi

    COMPILER_OPTS="-DCMAKE_CXX_COMPILER=${OMPI_CXX} -DCMAKE_C_COMPILER=${OMPI_CC} -DCMAKE_Fortran_COMPILER=${OMPI_FC}"
    echo "Compiler settings: ${COMPILER_OPTS}"
    
    cmake -DCMAKE_INSTALL_PREFIX=${PWD}/install \
          ${COMPILER_OPTS} \
          -DBUILD_SHARED_LIBS=ON \
          -DADIOS2_USE_MPI=${USE_MPI} \
          -DADIOS2_USE_Fortran=ON \
          -DADIOS2_USE_Python=OFF \
          -DCMAKE_BUILD_TYPE=Debug \
          -DBUILD_TESTING=ON \
          -DADIOS2_BUILD_EXAMPLES=ON \
          -DADIOS2_USE_DataMan=OFF \
          ${EXTERNAL_OPTS} \
          ${SRCDIR}
    
elif [ `hostname -f | cut -c 1-4` == "rhea" ]; then
    ############
    #  RHEA    #
    ############
    TARGET=`module list 2>&1 | grep "PE"- | sed "s/^.*PE-\([a-z]*\).*/\1/"`
    echo "Configure on Rhea for $TARGET env. MPI = ${USE_MPI}  use external packages = ${BUILD_FULL}"

    module load cmake/3.6.1

    unset EXTERNAL_OPTS
    if [ $BUILD_FULL == "true" ]; then
        # ADIOS 1.x
        EXTERNAL_OPTS+="-DADIOS2_USE_HDF5=OFF "
        # ZFP
        #      -DZFP_ROOT=/opt/zfp 
    else
        EXTERNAL_OPTS+="-DADIOS2_USE_HDF5=OFF"
    fi

    #COMPILER_OPTS="-DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_Fortran_COMPILER=ifort "
    COMPILER_OPTS="-DCMAKE_CXX_COMPILER=pgCC -DCMAKE_C_COMPILER=pgcc -DCMAKE_Fortran_COMPILER=pgf90 "
    #COMPILER_OPTS+="-DMPI_C_COMPILER=mpicc -DMPI_Fortran_COMPILER=mpif90"
    echo "Compiler settings: ${COMPILER_OPTS}"
    
    cmake -DCMAKE_INSTALL_PREFIX=${PWD}/install \
          ${COMPILER_OPTS} \
          -DBUILD_SHARED_LIBS=ON \
          -DADIOS2_USE_MPI=${USE_MPI} \
          -DADIOS2_USE_Fortran=ON \
          -DADIOS2_USE_Python=OFF \
          -DADIOS2_USE_SST=OFF \
          -DCMAKE_BUILD_TYPE=Debug \
          -DBUILD_TESTING=ON \
          -DADIOS2_BUILD_EXAMPLES=ON \
          -DADIOS2_USE_DataMan=OFF \
          ${EXTERNAL_OPTS} \
          ${SRCDIR}
    
else
    echo "Could not determine what machine is this."
    echo "This script is for configuring adios on the authors' machines."
    echo "You can study it to figure out how to configure adios on your system."
fi

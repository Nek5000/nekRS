This is the Matlab bindings to ADIOS2. It needs to be built manually and separately, after installing ADIOS2.

For having as little dependencies as possible, we suggest you build a lean package of ADIOS, no MPI, no Fortran, no staging, nothing extra and install it somewhere. (see the end of this README for an example build)

For the rest of this document, let's assume ADIOS is installed to a directory and environment variable ADIOS_INSTALL_DIR points to this directory.

Then edit the Makefile in the Matlab/ directory to point to that installation. And then run 

$ make matlab

The build requires that you have the MEX compiler ('mex') in the path. 

Finally, you need 
1) to set MATLABPATH to this Matlab build directory (or wherever you copy the .m and .mexa64 files) and 
2) if necessary, to force Matlab to use the same stdc++ library as what ADIOS was built with (in contrast to the Matlab executable). You can do this at startup of Matlab:

$ cd $ADIOS_INSTALL_DIR/lib
$ ldd libadios2_c.so
        libstdc++.so.6 => /usr/lib64/libstdc++.so.6 (0x00002adb84068000)

$ export LD_LIBRARY_PATH=$ADIOS_INSTALL_DIR/lib:$LD_LIBRARY_PATH
$ export MATLABPATH=<your Matlab/ directory>
$ LD_PRELOAD=/usr/lib64/libstdc++.so.6     matlab -nojvm -nodesktop -nosplash


In the Matlab/test directory there is a sample file and test script. The test function should complete without errors (provided your current directory is Matlab/test)

>>> TestADIOSRead 


For help with the Matlab bindings, use

>>> help adios
>>> help adiosopen
>>> help adiosread
>>> help adiosclose



Here is how to build the leanest ADIOS package:

$ cd <your adios source>
$ git pull
$ mkdir build.lean
$ cd build.lean
$ cmake -DCMAKE_INSTALL_PREFIX=/opt/adios/lean -DADIOS2_USE_MPI=OFF -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DCMAKE_BUILD_TYPE=RelWithDebInfo -DADIOS2_ENABLE_PIC=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DADIOS2_USE_SST=OFF -DADIOS2_USE_BZip2=OFF -DADIOS2_USE_DataMan=OFF  ..


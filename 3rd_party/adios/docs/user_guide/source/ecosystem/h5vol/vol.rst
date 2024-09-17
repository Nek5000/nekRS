**********
Disclaimer
**********

.. note::

   The Virtual Object Layer (VOL) is a feature introduced in recent release of  HDF5 1.12 (https://hdf5.wiki/index.php/New_Features_in_HDF5_Release_1.12).

   So please do make sure your HDF5 version supports the latest VOL.

Once the ADIOS VOL is compiled, There are two ways to apply it: 

* externally (through dynamic library, no code change) 
* internally (through client code). 

********
External
********

- Set the following environment parameters:

.. code-block:: c++

   HDF5_VOL_CONNECTOR=ADIOS2_VOL
   HDF5_PLUGIN_PATH=/replace/with/your/adios2_vol/lib/path/


Without code change, ADIOS2 VOL will be loaded at runtime by HDF5, to access ADIOS files without changing user code.

********
Internal
********

- include adios header 
- call the function to set VOL  when H5F is initiated
- call the function to unset VOL when H5F is closed

.. code-block:: c++

   // other includes
   #include <adios2/h5vol/H5Vol.h> // ADD THIS LINE TO INCLUDE ADIOS VOL

   hid_t  pid = H5Pcreate(H5P_FILE_ACCESS);
   // other declarations
   hid_t fid = H5Fopen(filename, mode, pid);

   H5VL_ADIOS2_set(pid); // ADD THIS LINE TO USE ADIOS VOL

   H5Fclose(fid);

   H5VL_ADIOS2_unset();  // ADD THIS LINE TO EXIT ADIOS VOL

..  To choose what ADIOS2 Engine to use, set env variable: ADIOS2_ENGINE (default is BP5)


**Note:** The following features are not supported in this VOL:

  * hyperslab support
  * HDF5   parameters are not passed down. e.g. compression/decompression
  * ADIOS2 parameters is not setup. 
  * user defined types
  * change of variable extent is not supported in ADIOS2. 







     

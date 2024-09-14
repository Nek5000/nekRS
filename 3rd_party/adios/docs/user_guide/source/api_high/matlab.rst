**********************
Matlab simple bindings
**********************

The ADIOS Matlab API supports reading data from ADIOS BP files with a 
simplified API that consists of three functions:

   * ``ADIOSOPEN``     returns a structure with information on an ADIOS BP File (variables and attributes).
   * ``ADIOSREAD``     reads in a variable from the file. It expects the info structure returned by ``ADIOSOPEN``.
   * ``ADIOSCLOSE``    closes the file.

Organization of an ADIOS BP file
--------------------------------

An ADIOS BP file contains a set of variables and attributes. Each variable in the group has a path, which defines a logical hierarchy of the variables within the file. 

Time dimension of a variable
----------------------------
Variables can be written several times from a program, if they have a time dimension. The reader exposes the variables with an extra dimension, i.e. a 2D variable written over time is seen as a 3D variable. In MATLAB, the extra dimension is the last dimension (the slowest changing dimension). Since the reader allows reading an arbitrary slice of a variable, data for one timestep can be read in with slicing.

Min/max of arrays
-----------------
The ADIOS BP format stores the min/max values in each variable.  The info structure therefore contains these min/max values. There is practically no overhead to provide this information (along with the values of all attributes) even for file sizes of several terabytes.


In the Matlab console use help for these functions

.. code-block:: matlabsession

   >>> help adiosopen
   >>> help adiosread
   >>> help adiosclose

``ADIOSOPEN``
-------------

``FILE = adiosopen(PATH)`` 
   Open a file for reading pointed by ``PATH`` and return an information structure (``FILE``). 

The returned FILE structure contains the following information

   .. code-block:: matlab

      Name              File path

      Handlers          Object handlers to pass on to ADIOS functions 
        FileHandler        uint64 file handler
        GroupHandler       uint64 IO group object handler
        ADIOSHandler       uint64 ADIOS object handler

      Variables         Structure array of variables
           Name            Path of variable
           Type            Matlab type class of data
           Dims            Array of dimensions
           StepsStart      First step's index for this variable in file, always at least 1
           StepsCount      Number of steps for this variable in file, always at least 1
           GlobalMin       Global minimum  of the variable (1-by-1 mxArray)
           GlobalMax       Global maximum of the variable
           
      Attribute         Structure array of attributes
           Name            Path of attribute
           Type            Matlab type class of data
           Value           Attribute value


``ADIOSREAD``
-------------

Read data from a BP file opened with ``adiosopen``. 
Provide the structure returned by ``adiosopen`` as the first input argument, 
and the path to a variable.
Inspect ``file.Variables`` and ``file.Attributes`` for the list of variables 
and attributes available in a file.

``data = adiosread(file, VARPATH)`` 

Read the entire variable ``VARPATH`` from a BP file. ``file`` is the output of ``ADIOSOPEN``. 
``VARPATH`` is a string to a variable or attribute. 
If an N-dimensional array variable has multiple steps in the file 
this function reads all steps and returns an N+1 dimensional array 
where the last dimension equals the number of steps.

``data = adiosread(file, INDEX)`` 

Read the entire variable from a BP file.
``INDEX`` points to a variable in the ``file.Variables`` array. 


``data = adiosread(..., START, COUNT, STEPSTART, STEPCOUNT)``

Read a portion of a variable. 

   .. code-block:: matlab
   
      START and COUNT:
      A slice is defined as two arrays of N integers, where N is the 
      number of dimensions of the variable, describing the
      "start" and "count" values. The "start" values start from 1.
          E.g. [1 5], [10 2] reads the first 10 values in the first dimension
      and 2 values from the 5th position in the second dimension resulting in
      a 10-by-2 array. 
          You can use negative numbers to index from the end of the array
      as in python. -1 refers to the last element of the array, -2 the one
      before and so on. 
          E.g. [-1], [1] reads in the last value of a 1D array. 
               [1], [-1] reads in the complete 1D array.

      STEPSTART and STEPCOUNT:
      Similarly, the number of steps from a specific step can be read instead
      of all data. Steps start from 1. Negative index can be used as well.
          E.g. -1, 1  will read in the last step from the file
               n, -1  will read all steps from 'n' to the last one

      
``ADIOSCLOSE``
--------------

``adiosclose(file)`` 
    Close file and free internal data structures. ``file`` is the structure returned by ``adiosopen``.
        
               

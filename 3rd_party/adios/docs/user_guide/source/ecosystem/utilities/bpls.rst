**********************
bpls : Inspecting Data
**********************

The ``bpls`` utility is for examining and pretty-printing the content of ADIOS output files (BP and HDF5 files).
By default, it lists the variables in the file including the type, name, and dimensionality. 

Let's assume we run the Heat Transfer tutorial example and produce the output by

.. code-block:: bash

    $ mpirun -n 12 ./heatSimulation  sim.bp  3 4   5 4   3 1
    Process decomposition  : 3 x 4
    Array size per process : 5 x 4
    Number of output steps : 3
    Iterations per step    : 1

    $ mpirun -n 3 ./heatAnalysis sim.bp a.bp 3 1

    $ bpls a.bp
      double   T     3*{15, 16}
      double   dT    3*{15, 16}

In our example, we have two arrays, ``T`` and ``dT``.
Both are 2-dimensional ``double`` arrays, their global size is ``15x16`` and the file contains ``3 output steps`` of these arrays.

.. note::

    bpls is written in C++ and therefore sees the order of the dimensions in *row major*. If the data was written from Fortran in column-major order, you will see the dimension order flipped when listing with bpls, just as a code written in C++ or python would see the data. 

 
Here is the description of the most used options
(use ``bpls -h`` to print help on all options for this utility).


* ``-l``

  Print the min/max of the arrays and the values of scalar values
  
  .. code-block:: bash

    $ bpls -l a.bp
      double   T     3*{15, 16} = 0 / 200
      double   dT    3*{15, 16} = -53.1922 / 49.7861


* ``-a`` ``-A``

  List the attributes along with the variables. ``-A`` will print the attributes only.

  .. code-block:: bash

    $ bpls a.bp -la
      double   T               3*{15, 16} = 0 / 200
      string   T/description   attr   = "Temperature from simulation"
      string   T/unit          attr   = "C"
      double   dT              3*{15, 16} = -53.1922 / 49.7861
      string   dT/description  attr   = "Temperature difference between two steps calculated in analysis"


* ``pattern``, ``-e`` 

  Select which variables/attributes to list or dump. By default the pattern(s) are interpreted as shell file patterns.

  .. code-block:: bash

    $ bpls a.bp -la T*
      double   T               3*{15, 16} = 0 / 200
      
  Multiple patterns can be defined in the command line. 

  .. code-block:: bash

    $ bpls a.bp -la T/* dT/* 
      string   T/description   attr   = "Temperature from simulation"
      string   T/unit          attr   = "C"
      string   dT/description  attr   = "Temperature difference between two steps calculated in analysis"

  If the -e option is given (all) the pattern(s) will be interpreted as regular expressions. 

  .. code-block:: bash

    $ bpls a.bp -la T.* -e
      double   T               3*{15, 16} = 0 / 200
      string   T/description   attr   = "Temperature from simulation"
      string   T/unit          attr   = "C"

* ``-D``

  Print the decomposition of a variable. In the BP file, the data blocks written by different writers are stored separately and have their own size info and min/max statistics. This option is useful at code development to check if the output file is written the way intended.


  .. code-block:: bash

    $ bpls a.bp -l T -D
      double   T               3*{15, 16} = 0 / 200
        step 0: 
          block 0: [ 0: 4,  0:15] = 3.54199e-14 / 200
          block 1: [ 5: 9,  0:15] = 58.3642 / 200
          block 2: [10:14,  0:15] = 0 / 200
        step 1: 
          block 0: [ 0: 4,  0:15] = 31.4891 / 153.432
          block 1: [ 5: 9,  0:15] = 68.2107 / 180.184
          block 2: [10:14,  0:15] = 31.4891 / 161.699
        step 2: 
          block 0: [ 0: 4,  0:15] = 48.0431 / 135.225
          block 1: [ 5: 9,  0:15] = 74.064 / 170.002
          block 2: [10:14,  0:15] = 48.0431 / 147.87

  In this case we find 3 blocks per output step and 3 output steps. We can see that the variable ``T`` was decomposed in the first (slow) dimension. In the above example, the ``T`` variable in the simulation output (``sim.bp``) had 12 blocks per step, but the analysis code was running on 3 processes, effectively reorganizing the data into fewer larger blocks.


* ``-d``

  Dump the data content of a variable. For pretty-printing, one should use the additional ``-n`` and ``-f`` options. For selecting only a subset of a variable, one should use the ``-s`` and ``-c`` options.

  By default, six values are printed per line and using C style ``-g`` prints for floating point values. 

  .. code-block:: bash

    $ bpls a.bp -d T
      double   T     3*{15, 16}
        (0, 0, 0)    124.925 124.296 139.024 95.2078 144.864 191.485
        (0, 0, 6)    139.024 140.814 124.925 109.035 110.825 58.3642
        (0, 0,12)    104.985 154.641 110.825 125.553 66.5603 65.9316
        ...
        (2,14, 4)    105.918 116.842 111.249 102.044 93.3121 84.5802
        (2,14,10)    75.3746 69.782 80.706 93.5492 94.7595 95.0709



  For pretty-printing, use the additional ``-n`` and ``-f`` options. 

  .. code-block:: bash

    $ bpls a.bp -d T -n 16 -f "%3.0f" 
      double   T     3*{15, 16}
        (0, 0, 0)    125 124 139  95 145 191 139 141 125 109 111  58 105 155 111 126
        (0, 1, 0)     67  66  81  37  86 133  81  82  67  51  52   0  47  96  52  67
        (0, 2, 0)    133 133 148 104 153 200 148 149 133 118 119  67 114 163 119 134
        ...
        (2,13, 0)     98  98  96  96 115 132 124 109  97  86  71  63  79  98  97  95
        (2,14, 0)     96  96  93  93 106 117 111 102  93  85  75  70  81  94  95  95


  For selecting a subset of a variable, use the ``-s`` and ``-c`` options. These options are N+1 dimensional for N-dimensional arrays with more than one steps. The first element of the options are used to select the starting step and the number of steps to print. 

  The following example dumps a ``4x4`` small subset from the center of the array, one step from the second (middle) step: 

  .. code-block:: bash

    $ bpls a.bp -d T -s "1,6,7" -c "1,4,4" -n 4
      double   T     3*{15, 16}
        slice (1:1, 6:9, 7:10)
        (1,6, 7)    144.09 131.737 119.383 106.787
        (1,7, 7)    145.794 133.44 121.086 108.49
        (1,8, 7)    145.794 133.44 121.086 108.49
        (1,9, 7)    144.09 131.737 119.383 106.787

* ``-y`` ``--noindex``

  Data can be dumped in a format that is easier to import later into other tools, like Excel. The leading array indexes can be omitted by using this option. Non-data lines, like the variable and slice info, are printed with a starting ``;``.

  .. code-block:: bash

    $ bpls a.bp -d T -s "1,6,7" -c "1,4,4" -n 4 --noindex
      ; double   T     3*{15, 16}
      ;   slice (1:1, 6:9, 7:10)
      144.09 131.737 119.383 106.787
      145.794 133.44 121.086 108.49
      145.794 133.44 121.086 108.49
      144.09 131.737 119.383 106.787


.. note::

  HDF5 files can also be dumped with bpls if ADIOS was built with HDF5 support. Note that the HDF5 files do not contain min/max information for the arrays and therefore bpls always prints 0 for them:


  .. code-block:: bash

    $ bpls -l a.h5
      double   T     3*{15, 16} = 0 / 0
      double   dT    3*{15, 16} = 0 / 0


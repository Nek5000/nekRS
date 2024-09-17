*****************
adios_reorganize
*****************

``adios_reorganize`` and ``adios_reorganize_mpi`` are generic ADIOS tools
to read in ADIOS streams and output the same data into another ADIOS stream.
The two tools are for serial and MPI environments, respectively.
They can be used for

* converting files between ADIOS BP and HDF5 format
* using separate compute nodes to stage I/O from/to disk to/from a large scale application
* reorganizing the data blocks for a different number of blocks

Let's assume we run the Heat Transfer tutorial example and produce the output by

.. code-block:: bash

    $ mpirun -n 12 ./heatSimulation  sim.bp  3 4   5 4   3 1
    Process decomposition  : 3 x 4
    Array size per process : 5 x 4
    Number of output steps : 3
    Iterations per step    : 1

    $ bpls sim.bp
      double   T     3*{15, 16}

In our example, we have an array, ``T``. It is a 2-dimensional ``double`` array, its global size is ``15x16`` and the file contains ``3 output steps`` of this array. The array is composed of 12 separate blocks coming from the 12 producers in the application. 

* Convert BP file to HDF5 file

   If ADIOS is built with HDF5 support, this tool can be used to convert between the two file formats.

   .. code-block:: bash
   
       $ mpirun -n 2 adios_reorganize_mpi sim.bp sim.h5 BPFile "" HDF5 "" 2 1
      
       $ bpls sim.h5
         double   T     3*{15, 16}

       $ h5ls -r sim.h5 
       /                        Group
       /Step0                   Group
       /Step0/T                 Dataset {15, 16}
       /Step1                   Group
       /Step1/T                 Dataset {15, 16}
       /Step2                   Group
       /Step2/T                 Dataset {15, 16}


       
* Stage I/O through extra compute nodes

    If writing data to disk is a bottleneck to the application, it may be worth to use extra nodes for receiving the data quickly from the application and then write to disk while the application continues computing. Similarly, data can be staged in from disk into extra nodes and make it available for fast read-in for an application. One can use one of the staging engines in ADIOS to perform this data staging (SST, SSC, DataMan).
    
    Assuming that the heatSimulation is using SST instead of file I/O in a run (set in its ``adios2.xml`` configuration file), staging to disk can be done this way:
    
    .. code-block:: bash
   
        Make sure adios2.xml sets SST for the simulation:
            <io name="SimulationOutput">
                <engine type="SST">
                </engine>
            </io>


        $ mpirun -n 12 ./heatSimulation  sim.bp  3 4   5 4   3 1 : \
                 -n 2 adios_reorganize_mpi sim.bp staged.bp SST "" BPFile "" 2 1
        
        $ bpls staged.bp
          double   T     3*{15, 16}

    
    Data is staged to the extra 2 cores and those will write the data to disk while the heatSimulation calculates the next step. Note, that this staging can only be useful if the tool can write all data to disk before the application produces the next output step. Otherwise, it will still block the application for I/O. 
    
    
* Reorganizing the data blocks in file for a different number of blocks

    In the above example, the application writes the array from 12 processes, but then ``adios_reorganize_mpi`` reads the global arrays on 2 processes. The output file on disk will therefore contain the array in 2 blocks. This reorganization of the array may be useful if reading is too slow for a dataset created by many-many processes. One may want to reorganize a file written by tens or hundreds of thousands of processes if one wants to read the content more than one time and the read time proves to be a bottleneck in one's work flow.
    
    .. code-block:: bash
    
        $ mpirun -n 12 ./heatSimulation  sim.bp  3 4   5 4   3 1
        $ bpls sim.bp -D
          double   T     3*{15, 16}
              step 0: 
                block  0: [ 0: 4,  0: 3]
                block  1: [ 5: 9,  0: 3]
                block  2: [10:14,  0: 3]
                block  3: [ 0: 4,  4: 7]
                block  4: [ 5: 9,  4: 7]
                block  5: [10:14,  4: 7]
                block  6: [ 0: 4,  8:11]
                block  7: [ 5: 9,  8:11]
                block  8: [10:14,  8:11]
                block  9: [ 0: 4, 12:15]
                block 10: [ 5: 9, 12:15]
                block 11: [10:14, 12:15]
              step 1: 
                block  0: [ 0: 4,  0: 3]
                block  1: [ 5: 9,  0: 3]
                block  2: [10:14,  0: 3]
                block  3: [ 0: 4,  4: 7]
                block  4: [ 5: 9,  4: 7]
                block  5: [10:14,  4: 7]
                block  6: [ 0: 4,  8:11]
                block  7: [ 5: 9,  8:11]
                block  8: [10:14,  8:11]
                block  9: [ 0: 4, 12:15]
                block 10: [ 5: 9, 12:15]
                block 11: [10:14, 12:15]
              step 2: 
                block  0: [ 0: 4,  0: 3]
                block  1: [ 5: 9,  0: 3]
                block  2: [10:14,  0: 3]
                block  3: [ 0: 4,  4: 7]
                block  4: [ 5: 9,  4: 7]
                block  5: [10:14,  4: 7]
                block  6: [ 0: 4,  8:11]
                block  7: [ 5: 9,  8:11]
                block  8: [10:14,  8:11]
                block  9: [ 0: 4, 12:15]
                block 10: [ 5: 9, 12:15]
                block 11: [10:14, 12:15]

          
        $ mpirun -n 2 adios_reorganize_mpi sim.bp reorg.bp BPFile "" BPFile "" 2 1
        $ bpls reorg.bp -D
          double   T     3*{15, 16}
              step 0: 
                block 0: [ 0: 6,  0:15]
                block 1: [ 7:14,  0:15]
              step 1: 
                block 0: [ 0: 6,  0:15]
                block 1: [ 7:14,  0:15]
              step 2: 
                block 0: [ 0: 6,  0:15]
                block 1: [ 7:14,  0:15]



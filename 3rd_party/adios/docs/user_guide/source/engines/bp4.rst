***
BP4 
***

The BP4 Engine writes and reads files in ADIOS2 native binary-pack (bp version 4) format. 
This was a new format for ADIOS 2.5 and improved on the metadata operations of the older BP3 format. 
Compared to the older format, BP4 provides three main advantages:

  * Fast and safe **appending** of multiple output steps into the same file. Better performance than writing new files each step. 
    Existing steps cannot be corrupted by appending new steps. 
  * **Streaming** through files (i.e. online processing). Consumer apps can read existing steps while the Producer is still writing new steps.
    Reader's loop can block (with timeout) and wait for new steps to arrive. Same reader code can read the entire data in post or in situ.
    No restrictions on the Producer.  
  * **Burst buffer support** for writing data. It can write the output to a local file system on each compute node and drain the data to the parallel file system in a separate asynchronous thread. 
    Streaming read from the target file system are still supported when data goes through the burst buffer. Appending to an existing file on the target file system is NOT supported currently.

BP4 files have the following structure given a "name" string passed as the first argument of ``IO::Open``:

.. code-block:: c++

   io.SetEngine("BP4");
   adios2::Engine bpFile = io.Open("name", adios2::Mode::Write);

will generate:

.. code-block:: bash

   % BP4 datasets are always a directory
   name.bp/

   % data and metadata files
   name.bp/
           data.0
           data.1
           ...
           data.M
           md.0
           md.idx

.. note::

   BP4 file names are compatible with the Unix (``/``) and Windows (``\\``) file system naming convention for directories and files.


This engine allows the user to fine tune the buffering operations through the following optional parameters:

1. **Profile**: turns ON/OFF profiling information right after a run

2. **ProfileUnits**: set profile units according to the required measurement scale for intensive operations

3. **Threads**: number of threads provided from the application for buffering, use this for very large variables in data size

4. **InitialBufferSize**: initial memory provided for buffering (minimum is 16Kb)

5. **BufferGrowthFactor**: exponential growth factor for initial buffer > 1, default = 1.05.

6. **MaxBufferSize**: maximum allowable buffer size (must be larger than 16Kb). If too large adios2 will throw an exception.

7. **FlushStepsCount**: users can select how often to produce the more expensive collective metadata file in terms of steps: default is 1. Increase to reduce adios2 collective operations footprint, with the trade-off of reducing checkpoint frequency. Buffer size will increase until first steps count if ``MaxBufferSize`` is not set.

8. **NumAggregators** (or **SubStreams**): Users can select how many sub-files (``M``) are produced during a run, ranges between 1 and the number of mpi processes from ``MPI_Size`` (``N``), adios2 will internally aggregate data buffers (``N-to-M``) to output the required number of sub-files. Default is 0, which will let adios2 to group processes per shared-memory-access (i.e. one per compute node) and use one process per node as an aggregator. If NumAggregators is larger than the number of processes then it will be set to the number of processes.

9. **AggregatorRatio**: An alternative option to NumAggregators to pick every Nth process as aggregator. An integer divider of the number of processes is required, otherwise a runtime exception is thrown. 

10. **OpenTimeoutSecs**: (Streaming mode) Reader may want to wait for the creation of the file in ``io.Open()``. By default the Open() function returns with an error if file is not found.

11. **BeginStepPollingFrequencySecs**: (Streaming mode) Reader can set how frequently to check the file (and file system) for new steps. Default is 1 seconds which may be stressful for the file system and unnecessary for the application.

12. **StatsLevel**: Turn on/off calculating statistics for every variable (Min/Max). Default is On. It has some cost to generate this metadata so it can be turned off if there is no need for this information.

13. **StatsBlockSize**: Calculate Min/Max for a given size of each process output. Default is one Min/Max per writer. More fine-grained min/max can be useful for querying the data. 

14. **NodeLocal** or **Node-Local**: For distributed file system. Every writer process must make sure the .bp/ directory is created on the local file system. Required when writing to local disk/SSD/NVMe in a cluster. Note: the BurstBuffer* parameters are newer and should be used for using the local storage as temporary instead of this parameter.

15. **BurstBufferPath**: Redirect output file to another location and drain it to the original target location in an asynchronous thread. It requires to be able to launch one thread per aggregator (see SubStreams) on the system. This feature can be used on machines that have local NVMe/SSDs on each node to accelerate the output writing speed. On Summit at OLCF, use "/mnt/bb/<username>" for the path where <username> is your user account name. Temporary files on the accelerated storage will be automatically deleted after the application closes the output and ADIOS drains all data to the file system, unless draining is turned off (see the next parameter). Note: at this time, this feature cannot be used to append data to an existing dataset on the target system. 

16. **BurstBufferDrain**: To write only to the accelerated storage but to not drain it to the target file system, set this flag to false. Data will NOT be deleted from the accelerated storage on close. By default, setting the BurstBufferPath will turn on draining. 

17. **BurstBufferVerbose**: Verbose level 1 will cause each draining thread to print a one line report at the end (to standard output) about where it has spent its time and the number of bytes moved. Verbose level 2 will cause each thread to print a line for each draining operation (file creation, copy block, write block from memory, etc). 

18. **StreamReader**: By default the BP4 engine parses all available metadata in Open(). An application may turn this flag on to parse a limited number of steps at once, and update metadata when those steps have been processed. If the flag is ON, reading only works in streaming mode (using BeginStep/EndStep); file reading mode will not work as there will be zero steps processed in Open().

============================== ===================== ===========================================================
 **Key**                       **Value Format**      **Default** and Examples
============================== ===================== ===========================================================
 Profile                        string On/Off         **On**, Off
 ProfileUnits                   string                **Microseconds**, Milliseconds, Seconds, Minutes, Hours
 Threads                        integer > 1           **1**, 2, 3, 4, 16, 32, 64
 InitialBufferSize              float+units >= 16Kb   **16Kb**, 10Mb, 0.5Gb
 MaxBufferSize                  float+units >= 16Kb   **at EndStep**, 10Mb, 0.5Gb
 BufferGrowthFactor             float > 1             **1.05**, 1.01, 1.5, 2
 FlushStepsCount                integer > 1           **1**, 5, 1000, 50000
 NumAggregators                 integer >= 1          **0 (one file per compute node)**, ``MPI_Size``/2, ... , 2, (N-to-1) 1
 AggregatorRatio                integer >= 1          not used unless set, ``MPI_Size``/N must be an integer value
 OpenTimeoutSecs                float                 **0**, ``10.0``, ``5``
 BeginStepPollingFrequencySecs  float                 **1**, ``10.0`` 
 StatsLevel                     integer, 0 or 1       **1**, ``0``
 StatsBlockSize                 integer > 0           **a very big number**, ``1073741824`` for blocks with 1M elements
 NodeLocal                      string On/Off         **Off**, On
 Node-Local                     string On/Off         **Off**, On
 BurstBufferPath                string                **""**, /mnt/bb/norbert, /ssd
 BurstBufferDrain               string On/Off         **On**, Off
 BurstBufferVerbose             integer, 0-2          **0**, ``1``, ``2`` 
 StreamReader                   string On/Off         On, **Off**
============================== ===================== ===========================================================


Only file transport types are supported. Optional parameters for ``IO::AddTransport`` or in runtime config file transport field:

**Transport type: File**

============= ================= ================================================
 **Key**       **Value Format**  **Default** and Examples
============= ================= ================================================
 Library           string        **POSIX** (UNIX), **FStream** (Windows), stdio, IME
============= ================= ================================================

The IME transport directly reads and writes files stored on DDN's IME burst
buffer using the IME native API. To use the IME transport, IME must be
avaiable on the target system and ADIOS2 needs to be configured with
``ADIOS2_USE_IME``. By default, data written to the IME is automatically
flushed to the parallel filesystem at every ``EndStep()`` call. You can
disable this automatic flush by setting the transport parameter ``SyncToPFS``
to ``OFF``.

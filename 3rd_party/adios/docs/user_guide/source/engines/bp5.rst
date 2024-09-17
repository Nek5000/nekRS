***
BP5
***

The BP5 Engine writes and reads files in ADIOS2 native binary-pack (bp version 5) format. 
This was a new format for ADIOS 2.8, improving on the metadata operations and the memory consumption 
of the older BP4/BP3 formats. BP5 is the default file format as of
ADIOS 2.9.  As compared to the older format, BP5 provides three main advantages:

  * **Lower memory** consumption. Deferred Puts will use user buffer for I/O wherever possible thus saving on a memory copy. 
    Aggregation uses a fixed-size shared-memory segment on each compute node instead of using MPI to send data from one process to another. 
    Memory consumption can get close to half of BP4 in some cases. 
  * **Faster metadata** management improves write/read performance where hundreds or more variables are added to the output. 
  * Improved functionality around **appending** many output steps into the same file. Better performance than writing new files each step. 
    Restart can append to an existing series by truncating unwanted steps. Readers can filter out unwanted steps to only see and process a 
    limited set of steps. Just like as in BP4, existing steps cannot be corrupted by appending new steps.

In 2.8 BP5 was a brand new file format and engine. It still does **NOT** support some functionality of BP4:

  * **Burst buffer support** for writing data.

BP5 files have the following structure given a "name" string passed as the first argument of ``IO::Open``:

.. code-block:: c++

   io.SetEngine("BP5");
   adios2::Engine bpFile = io.Open("name", adios2::Mode::Write);

will generate:

.. code-block:: bash

   % BP5 datasets are always a directory
   name.bp/

   % data and metadata files
   name.bp/
           data.0
           data.1
           ...
           data.M
           md.0
           md.idx
           mmd.0 

.. note::

   BP5 file names are compatible with the Unix (``/``) and Windows (``\\``) file system naming convention for directories and files.

.. note::

   BP5 has an ``mmd.0`` file in the directory, which BP4 does not have.


This engine allows the user to fine tune the buffering operations through the following optional parameters:

1. Streaming through file

   1. **OpenTimeoutSecs**: (Streaming mode) Reader may want to wait for the creation of the file in ``io.Open()``. By default the Open() function returns with an error if file is not found.

   #. **BeginStepPollingFrequencySecs**: (Streaming mode) Reader can set how frequently to check the file (and file system) for new steps. Default is 1 seconds which may be stressful for the file system and unnecessary for the application.

#. Aggregation

   #. **AggregationType**: *TwoLevelShm*, *EveryoneWritesSerial* and *EveryoneWrites* are three aggregation strategies. See :ref:`Aggregation in BP5`. The default is *TwoLevelShm*.
 
   #. **NumAggregators**: The number of processes that will ever write data directly to storage. The default is set to the number of compute nodes the application is running on (i.e. one process per compute node). TwoLevelShm will select a fixed number of processes *per compute-node* to get close to the intention of the user but does not guarantee the exact number of aggregators.

   #. **AggregatorRatio**: An alternative option to NumAggregators to pick every nth process as aggregator. The number of aggregators will be automatically kept to be within 1 and total number of processes no matter what bad number is supplied here. Moreover, TwoLevelShm will select an fixed number of processes *per compute-node* to get close to the intention of this ratio but does not guarantee the exact number of aggregators.

   #. **NumSubFiles**: The number of data files to write to in the *.bp/* directory. Only used by *TwoLevelShm* aggregator, where the number of files can be smaller then the number of aggregators. The default is set to *NumAggregators*. 

   #. **StripeSize**: The data blocks of different processes are aligned to this size (default is 4096 bytes) in the files. Its purpose is to avoid multiple processes to write to the same file system block and potentially slow down the write.  

   #. **MaxShmSize**: Upper limit for how much shared memory an aggregator process in *TwoLevelShm* can allocate. For optimum performance, this should be at least *2xM +1KB* where *M* is the maximum size any process writes in a single step. However, there is no point in allowing for more than 4GB. The default is 4GB.


#. Buffering

   #. **BufferVType**: *chunk* or *malloc*, default is chunking. Chunking maintains the buffer as a list of memory blocks, either ADIOS-owned for sync-ed Puts and small Puts, and user-owned pointers of deferred Puts. Malloc maintains a single memory block and extends it (reallocates) whenever more data is buffered. Chunking incurs extra cost in I/O by having to write data in chunks (multiple write system calls), which can be helped by increasing *BufferChunkSize* and *MinDeferredSize*. Malloc incurs extra cost by reallocating memory whenever more data is buffered (by Put()), which can be helped by increasing *InitialBufferSize*. 

   #. **BufferChunkSize**: (for *chunk* buffer type) The size of each memory buffer chunk, default is 128MB but it is worth increasing up to 2147381248 (a bit less than 2GB) if possible for maximum write performance.

   #. **MinDeferredSize**: (for *chunk* buffer type) Small user variables are always buffered, default is 4MB. 

   #. **InitialBufferSize**: (for *malloc* buffer type) initial memory provided for buffering (default and minimum is 16Kb). To avoid reallocations, it is worth increasing this size to the expected maximum total size of data any process would write in any step (not counting deferred Puts). 

   #. **GrowthFactor**: (for *malloc* buffer type) exponential growth factor for initial buffer > 1, default = 1.05.
      
#. Managing steps

   #. **AppendAfterSteps**: BP5 enables overwriting some existing steps by opening in *adios2::Mode::Append* mode and specifying how many existing steps to keep. Default value is MAX_INT, so it always appends after the last step. -1 would achieve the same thing. If you have 10 steps in the file,

      - value 0 means starting from the beginning, truncating all existing data
      - value 1 means appending after the first step, so overwrite 2,3...10
      - value 10 means appending after all existing steps
      - value >10 means the same, append after all existing steps (gaps in steps are impossible)
      - -1 means appending after the last step, i.e. same as 10 or higher
      - -2 means removing the last step, i.e. starting from the 10th
      - -11 (and <-11) means truncating all existing data
  
   #. **SelectSteps**: BP5 reading allows for only seeing selected steps. This is a string of space-separated list of range definitions in
      the form of "start:end:step". Indexing starts from 0. If 'end' is 'n' or 'N', then it is an unlimited range expression. Range definitions are adding up. Note that in the reading functions, counting the steps is *always* *0* to *s-1* where *s* steps are presented, so even after applying this selection, the selected steps are presented as *0* to *s-1*. Examples:

      - "0 6 3 2" selects four steps indexed 0,2,3 and 6 (presented in reading as 0,1,2,3)
      - "1:5" selects 5 consecutive steps, skipping step 0, and starting from 1
      - "2:n" selects all steps from step 2
      - "0:n:2" selects every other steps from the beginning (0,2,4,6...)
      - "0:n:3  10:n:5" selects every third step from the beginning and additionally every fifth steps from step 10.

#. Asynchronous writing I/O

   #. **AsyncOpen**: *true/false* Call the open function asynchronously. It decreases I/O overhead when creating lots of subfiles (*NumAggregators* is large) and one calls *io.Open()* well ahead of the first write step. Only implemented for writing. Default is *true*.

   #. **AsyncWrite**: *true/false* Perform data writing operations asynchronously after *EndStep()*. Default is *false*. If the application calls *EnterComputationBlock()/ExitComputationBlock()* to indicate phases where no communication is happening, ADIOS will try to perform all data writing during those phases, otherwise it will write immediately and eagerly after *EndStep()*. 
   
#. Direct I/O. Experimental, see discussion on `GitHub <https://github.com/ornladios/ADIOS2/issues/3029>`_.
 
   #. **DirectIO**: Turn on O_DIRECT when using POSIX transport. Do not use this on parallel file systems. 

   #. **DirectIOAlignOffset**: Alignment for file offsets. Default is 512 which is usually 

   #. **DirectIOAlignBuffer**: Alignment for memory pointers. Default is to be same as *DirectIOAlignOffset*. 

#. Miscellaneous

   #. **StatsLevel**: 1 turns on *Min/Max* calculation for every variable, 0 turns this off. Default is 1. It has some cost to generate this metadata so it can be turned off if there is no need for this information.

   #. **MaxOpenFilesAtOnce**: Specify how many subfiles a process can keep open at once. Default is unlimited. If a dataset contains more subfiles than how many open file descriptors the system allows (see *ulimit -n*) then one can either try to raise that system limit (set it with *ulimit -n*), or set this parameter to force the reader to close some subfiles to stay within the limits.
   
   #. **Threads**: Read side: Specify how many threads one process can use to speed up reading. The default value is *0*, to let the engine estimate the number of threads based on how many processes are running on the compute node and how many hardware threads are available on the compute node but it will use maximum 16 threads. Value *1* forces the engine to read everything within the main thread of the process. Other values specify the exact number of threads the engine can use. Although multithreaded reading works in a single *Get(adios2::Mode::Sync)* call if the read selection spans multiple data blocks in the file, the best parallelization is achieved by using deferred mode and reading everything in *PerformGets()/EndStep()*.   

   #. **FlattenSteps**: This is a writer-side parameter specifies that the
      reader should interpret multiple writer-created timesteps as a
      single timestep, essentially flattening all Put()s into a single step.

   #. **IgnoreFlattenSteps**: This is a reader-side parameter that
      tells the reader to ignore any FlattenSteps parameter supplied
      to the writer.

============================== ===================== ===========================================================
 **Key**                       **Value Format**      **Default** and Examples
============================== ===================== ===========================================================
 OpenTimeoutSecs                float                 **0** for *ReadRandomAccess* mode, **3600** for *Read* mode, ``10.0``, ``5``
 BeginStepPollingFrequencySecs  float                 **1**, 10.0 
 AggregationType                string                **TwoLevelShm**, EveryoneWritesSerial, EveryoneWrites
 NumAggregators                 integer >= 1          **0 (one file per compute node)**
 AggregatorRatio                integer >= 1          not used unless set
 NumSubFiles                    integer >= 1          **=NumAggregators**, only used when *AggregationType=TwoLevelShm*
 StripeSize                     integer+units         **4KB**
 MaxShmSize                     integer+units         **4294762496**
 BufferVType                    string                **chunk**, malloc
 BufferChunkSize                integer+units         **128MB**, worth increasing up to min(2GB, datasize/process/step)
 MinDeferredSize                integer+units         **4MB**
 InitialBufferSize              float+units >= 16Kb   **16Kb**, 10Mb, 0.5Gb
 GrowthFactor                   float > 1             **1.05**, 1.01, 1.5, 2
 AppendAfterSteps               integer >= 0          **INT_MAX**
 SelectSteps                    string                "0 6 3 2", "1:5", "0:n:3  10:n:5"
 AsyncOpen                      string On/Off         **On**, Off, true, false
 AsyncWrite                     string On/Off         **Off**, On, true, false
 DirectIO                       string On/Off         **Off**, On, true, false
 DirectIOAlignOffset            integer >= 0          **512**
 DirectIOAlignBuffer            integer >= 0          set to DirectIOAlignOffset if unset
 StatsLevel                     integer, 0 or 1       **1**, 0
 MaxOpenFilesAtOnce             integer >= 0          **UINT_MAX**, 1024, 1
 Threads                        integer >= 0          **0**, 1, 32
 FlattenSteps                   boolean               **off**, on, true, false
 IgnoreFlattenSteps             boolean               **off**, on, true, false
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

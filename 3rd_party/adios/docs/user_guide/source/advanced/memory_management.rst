###################
 Memory Management
###################

BP4 buffering
-------------

BP4 has a simple buffering mechanism to provide ultimate performance at the cost of high memory usage: all user data (passed in `Put()` calls) is buffered in one contiguous memory allocation and writing/aggregation is done with this large buffer in `EndStep()`. Aggregation in BP4 uses MPI to send this buffer to the aggregator and hence maintaining two such large buffers. Basically, if an application writes `N` bytes of data in a step, then BP4 needs approximately `2xN` bytes extra memory for buffering. 

A potential performance problem is that BP4 needs to extend the buffer occasionally to fit more incoming data (more `Put()` calls). At large sizes the reallocation may need to move the buffer into a different place in memory, which requires copying the entire existing buffer. When there are GBs of data already buffered, this copy will noticably decrease the overall observed write performance. This situation can be avoided if one can guess a usable upper limit to how much data each process is going to write, and telling this to the BP4 engine through the **InitialBufferSize** parameter before `Open()`.

Another potential problem is that reallocation may fail at some point, well before the limits of memory, since it needs a single contiguous allocation be available.

BP5 buffering
-------------

BP5 is designed to use less memory than BP4. The buffer it manages is a list of large chunks. The advantages of the list of chunks is that no reallocation of existing buffer is needed, and that BP5 can potentially allocate more buffer than BP4 since it requests many smaller chunks instead of a large contiguous buffer. In general, chunks should be as big as the system/application can afford, up to **2147381248** bytes (almost but less than 2GB, the actual size limit POSIX write() calls have). Each chunk will result in a separate write call, hence minimizing the number of chunks is preferred. The current default is set to 128MB, so please increase this on large computers if you can and if you write more than that amount of data per process, using the parameter **BufferChunkSize**. 

Second, BP5 can add a large user variable as a chunk to this list without copying it at all and use it directly to write (or send to aggregator). `Put(..., adios2::Mode::Deferred)` will handle the user data directly, unless its size is below a threshold (see parameter **MinDeferredSize**). 

.. note::
    Do not call `PerformPuts()` when using BP5, because this call forces copying all user data into the internal buffer before writing, eliminating all benefits of zero-copy that BP5 provides when operating with large buffers.  Instead, consider using Put() with the Sync option if you want to force ADIOS to copy data immediately.  Alternatively, BP5 offers PerformDataWrite(), an collective operation that actually moves data to storage, potentially freeing up buffer and application memory.

Third, BP5 is using a shared memory segment on each compute node for aggregation, instead of MPI. The best settings for the shared memory is 4GB (see parameter **MaxShmSize**), enough place for two chunks with the POSIX write limit. More is useless but can be smaller if a system/application cannot allow this much space for aggregation (but there will be more write calls to disk as a result).

Span object in internal buffer
------------------------------

Another option to decrease memory consumption is to pre-allocate space in the BP4/BP5 buffer and then prepare output variables directly in that space. This will avoid a copy and the need for doubling memory for temporary variables that are only created for output purposes. This Span feature is only available in C++. 
See the `Span()` function in :ref:`C++11 Engine class`  ../api_full/api_full.html#engine-class

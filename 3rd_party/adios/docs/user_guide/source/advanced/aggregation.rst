#############
 Aggregation 
#############

The basic problem of large-scale I/O is that the N-to-1 and N-to-N (process-to-file) patterns do not scale and one must set the number of files in an output to the capability of the file system, not the size of the application. Hence, *N* processes need to write to *M* files to 

    1) utilize the bandwidth of the file system and to 
    2) minimize the cost of multiple process writing to a single file, while
    3) not overwhelming the file system with too many files.


Aggregation in BP5
-------------------

There are two implementations of aggregation in BP5, none of them is the same as the one in BP4. The aggregation setup in ADIOS2 consist of: a) *NumAggregators*, which processes do write to disk (others will send data to them), and b) *NumSubFiles*, how many files they will write. 

**EveryoneWritesSerial** is a simple aggregation strategy. Every process is writing its own data to disk, to one particular file only, and the processes are serialized over each particular file. In this aggregator, *NumAggregators* = *NumSubFiles* (= *M*). This approach should scale well with application size. On Summit's GPFS though we observe that a single writer per compute node is better than multiple process writing to the file system, hence this aggregation method performs poorly there.

**EveryoneWrites** is the same strategy as the previous except that every process immediately writes its own data to its designated file. Since it basically implements an N-to-N write pattern, this method does not scale, so only use it up to a moderate number of processes (1-4 process * number of file system servers). At small scale, as long as the file system can deal with the on-rush of the write requests, this method can provide the fastest I/O. 

**TwoLevelShm** has a subset of processes that actually write to disk (*NumAggregators*). There must be at least one process per compute node, which creates a shared-memory segment for other processes on the node to send their data. The aggregator process basically serializes the writing of data from this subset of processes (itself and the processes that send data to it). TwoLevelShm performs similarly to EveryoneWritesSerial on Lustre, and is the only good option on Summit's GPFS. 

The number of files (*NumSubFiles*) can be smaller than *NumAggregators*, and then multiple aggregators will write to one file concurrently. Such a setup becomes useful when the number of nodes is many times more than the number of file servers.

TwoLevelShm works best if each process's output data fits into the shared-memory segment, which holds two pages. Since POSIX writes are limited to about 2GB, the best setup is to use 4GB shared-memory size by each aggregator. This is the default size, but you can use the *MaxShmSize* parameter to set this lower if necessary. At runtime, BP5 will only allocate twice the maximum size of the largest data size any process has, but up to MaxShmSize. If the data from two processes does not fit into the shared-memory segment, BP5 will need to perfom multiple iterations of copy and disk-write, which is generally slower than writing large data blocks at once.  

The **default setup** is *TwoLevelShm*, where *NumAggregators* is the number of compute nodes the application is running on, and the number of files is the same. This setup is good for Summit's GPFS and good for Lustre at large scale. However, the default setup leaves potential performance on the table when running applications at smaller scale, where the one process per node setup cannot utilize the full bandwidth of a large parallel file system. 

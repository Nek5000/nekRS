******************************************
DataMan for Wide Area Network Data Staging
******************************************

The DataMan engine is designed for data staging over the wide area network.
It is supposed to be used in cases where a few writers send data to a few readers
over long distance.

DataMan supports compression operators such as ZFP lossy compression and BZip2 lossless compression.
Please refer to the operator section for usage.

The DataMan engine takes the following parameters:

1. ``IPAddress``: No default value. The IP address of the host where the writer application runs.
   This parameter is compulsory in wide area network data staging.

2. ``Port``: Default **50001**. The port number on the writer host that will be used for data transfers.

3. ``Timeout``: Default **5**. Timeout in seconds to wait for every send / receive operation.
   Packages not sent or received within this time are considered lost.

4. ``RendezvousReaderCount``: Default **1**. This integer value specifies the number of readers for which the writer should wait before the writer-side Open() returns.
   By default, an early-starting writer will wait for the reader to start, or vice versa.
   A number >1 will cause the writer to wait for more readers, and a value of 0 will allow the writer to proceed without any readers present.
   This value is interpreted by DataMan Writer engines only.

5. ``Threading``: Default **true** for reader, **false** for writer. Whether to use threads for send and receive operations.
   Enabling threading will cause extra overhead for managing threads and buffer queues, but will improve the continuity of data steps for readers, and help overlap data transfers with computations for writers.

6. ``TransportMode``: Default **fast**. Only DataMan writers take this parameter.
   Readers are automatically synchronized at runtime to match writers' transport mode.
   The fast mode is optimized for latency-critical applications.
   It enforces readers to only receive the latest step.
   Therefore, in cases where writers are faster than readers, readers will skip some data steps.
   The reliable mode ensures that all steps are received by readers, by sacrificing performance compared to the fast mode.

7. ``MaxStepBufferSize``: Default **128000000**. In order to bring down the latency in wide area network staging use cases, DataMan uses a fixed receiver buffer size.
   This saves an extra communication operation to sync the buffer size for each step, before sending actual data.
   The default buffer size is 128 MB, which is sufficient for most use cases.
   However, in case 128 MB is not enough, this parameter must be set correctly, otherwise DataMan will fail.


=============================== ================== ================================================
 **Key**                         **Value Format**   **Default** and Examples
=============================== ================== ================================================
 IPAddress                       string             **N/A**, 22.195.18.29
 Port                            integer            **50001**, 22000, 33000
 Timeout                         integer            **5**, 10, 30
 RendezvousReaderCount           integer            **1**, 0, 3
 Threading                       bool               **true** for reader, **false** for writer
 TransportMode                   string             **fast**, reliable
 MaxStepBufferSize               integer            **128000000**, 512000000, 1024000000
=============================== ================== ================================================



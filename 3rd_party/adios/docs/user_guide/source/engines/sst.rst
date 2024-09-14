*********************************
SST Sustainable Staging Transport
*********************************

In ADIOS2, the Sustainable Staging Transport (SST) is an engine that allows
direct connection of data producers and consumers via the ADIOS2 write/read
APIs.  This is a classic streaming data architecture where the data passed
to ADIOS on the write side (via Put() deferred and sync, and similar calls)
is made directly available to a reader (via Get(), deferred and sync, and
similar calls).

SST is designed for use in HPC environments and can take advantage of RDMA
network interconnects to speed the transfer of data between communicating
HPC applications; however, it is also capable of operating in a Wide Area
Networking environment over standard sockets.  SST supports full MxN data
distribution, where the number of reader ranks can differ from the number of
writer ranks.  SST also allows multiple reader cohorts to get access to a writer's
data simultaneously.

To use this engine, you can either specify it in your xml config file, with
tag ``<engine type=SST>`` or, set it in client code. For example, here is
how to create an SST reader:

.. code-block:: c++

 adios2::IO sstIO = adios.DeclareIO("SomeName");
 sstIO.SetEngine("SST");
 adios2::Engine sstReader = sstIO.Open(filename, adios2::Mode::Read);

and a sample code for SST writer is:

.. code-block:: c++

 adios2::IO sstIO = adios.DeclareIO("SomeName");
 sstIO.SetEngine("SST");
 adios2::Engine sstWriter = sstIO.Open(filename, adios2::Mode::Write);

The general goal of ADIOS2 is to ease the conversion of a file-based
application to instead use a non-file streaming interconnect, for example,
data producers such as computational physics codes and consumers such as
analysis applications.  However, there are some uses of ADIOS2 APIs that
work perfectly well with the ADIOS2 file engines, but which will not work or
will perform badly with streaming.  For example, SST is based upon the *"step"* concept and
ADIOS2 applications that use SST must call ``BeginStep()`` and ``EndStep()``.  On
the writer side, the ``Put()`` calls between ``BeginStep`` and ``EndStep`` are the unit
of communication and represent the data that will be available between the
corresponding ``Begin``/``EndStep`` calls on the reader.

Also, it is recommended that SST-based applications not use the ADIOS2
Get() sync method unless there is only one data item to be read per step.
This is because SST implements MxN data transfer (and avoids having to
deliver all data to every reader), by queueing data on the writer ranks
until it is known which reader rank requires it.  Normally this data fetch
stage is initiated by ``PerformGets()`` or ``EndStep()``, both of which fulfill any
pending ``Get()`` deferred operations.  However, unlike ``Get()`` deferred, the
semantics of ``Get()`` sync require the requested data to be fetched from the
writers before the call can return.   If there are multiple calls to
``Get()`` sync per step, each one may require a communication with many writers,
something that would have only had to happen once if ``Get()`` differed were used
instead.  Thus the use of ``Get()`` sync is likely to incur a substantial
performance penalty.

On the writer side, depending upon the chosen data marshaling option there
may be some (relatively small) performance differences between ``Put()`` sync and
``Put()`` deferred, but they are unlikely to be as substantial as between
``Get()`` sync and ``Get()`` deferred.

Note that SST readers and writers do not necessarily move in lockstep, but
depending upon the queue length parameters and queueing policies specified,
differing reader and writer speeds may cause one or the other side to wait
for data to be produced or consumed, or data may be dropped if allowed by
the queueing policy.  However, steps themselves are atomic and no step will
be partially dropped, delivered to a subset of ranks, or otherwise divided.

The SST engine allows the user to customize the streaming operations through
the following optional parameters:

1. ``RendezvousReaderCount``: Default **1**.  This integer value specifies
the number of readers for which the writer should wait before the
writer-side Open() returns.   The default of 1 implements an ADIOS1/flexpath
style "rendezvous", in which an early-starting reader will wait for the
writer to start, or vice versa.  A number >1 will cause the writer to wait
for more readers and a value of 0 will allow the writer to proceed without
any readers present.  This value is interpreted by SST Writer engines only.

2. ``RegistrationMethod``:  Default **"File"**.  By default, SST reader and
writer engines communicate network contact information via files in a shared
filesystem.  Specifically, the ``"filename"`` parameter in the ``Open()`` call is
interpreted as a path which the writer uses as the name of a file to which
contact information is written, and from which a reader will attempt to read
contact information.  As with other file-based engines, file creation and
access is subject to the usual considerations (directory components are
interpreted, but must exist and be traversable, writer must be able to
create the file and the reader must be able to read it).  Generally the file
so created will exist only for as long as the writer keeps the stream
Open(), but abnormal process termination may leave "stale" files in those
locations.  These stray ".sst" files should be deleted to avoid confusing
future readers.  SST also offers a **"Screen"** registration method in which
writers and readers send their contact information to, and read it from,
stdout and stdin respectively.  The "screen" registration method doesn't
support batch mode operations in any way, but may be useful when manually
starting jobs on machines in a WAN environment that don't share a
filesystem. A future release of SST will also support a **"Cloud"**
registration method where contact information is registered to and retrieved
from a network-based third-party server so that both the shared filesystem
and interactivity can be avoided. This value is interpreted by both SST
Writer and Reader engines.

3. ``QueueLimit``:  Default **0**.  This integer value specifies the number
of steps which the writer will allow to be queued before taking specific
action (such as discarding data or waiting for readers to consume the
data).  The default value of 0 is interpreted as no limit.  This value is
interpreted by SST Writer engines only.

4. ``QueueFullPolicy``: Default **"Block"**.  This value controls what
policy is invoked if a non-zero **QueueLimit** has been specified and
new data would cause the queue limit to be reached.  Essentially, the
**"Block"** option ensures data will not be discarded and if the queue
fills up the writer will block on **EndStep** until the data has been
read. If there is one active reader, **EndStep** will block until data
has been consumed off the front of the queue to make room for newly
arriving data.  If there is more than one active reader, it is only
removed from the queue when it has been read by all readers, so the
slowest reader will dictate progress.  **NOTE THAT THE NO READERS
SITUATION IS A SPECIAL CASE**: If there are no active readers, new
timesteps are considered to have completed their active queueing
immediately upon submission.  They may be retained in the "reserve
queue" if the ReserveQueueLimit is non-zero.  However, if that
ReserveQueueLimit parameter is zero, timesteps submitted when there
are no active readers will be immediately discarded.

Besides **"Block"**, the other
acceptable value for **QueueFullPolicy** is **"Discard"**.  When
**"Discard"** is specified, and an **EndStep** operation would add
more than the allowed number of steps to the queue, some step is
discarded.  If there are no current readers connected to the stream,
the *oldest* data in the queue is discarded.  If there are current
readers, then the *newest* data (I.E. the just-created step) is
discarded.  (The differential treatment is because SST sends metadata
for each step to the readers as soon as the step is accepted and
cannot reliably prevent that use of that data without a costly
all-to-all synchronization operation.  Discarding the *newest* data
instead is less satisfying, but has a similar long-term effect upon
the set of steps delivered to the readers.)  This value is interpreted
by SST Writer engines only.

5. ``ReserveQueueLimit``:  Default **0**.  This integer value specifies the
number of steps which the writer will keep in the queue for the benefit
of late-arriving readers.  This may consist of timesteps that have
already been consumed by any readers, as well as timesteps that have not
yet been consumed.  In some sense this is target queue minimum size,
while QueueLimit is a maximum size.  This value is interpreted by SST
Writer engines only. 

6. ``DataTransport``: Default **varies**.  This string value specifies
the underlying network communication mechanism to use for exchanging
data in SST.  Generally this is chosen by SST based upon what is
available on the current platform.  However, specifying this engine
parameter allows overriding SST's choice.  Current allowed values are
**"UCX"**, **"MPI"**, **"RDMA"**, and **"WAN"**.  (**ib** and **fabric** are accepted as
equivalent to **RDMA** and **evpath** is equivalent to **WAN**.)
Generally both the reader and writer should be using the same network
transport, and the network transport chosen may be dictated by the
situation.  For example, the RDMA transport generally operates only
between applications running on the same high-performance interconnect
(e.g. on the same HPC machine).  If communication is desired between
applications running on different interconnects, the Wide Area Network
(WAN) option should be chosen.  This value is interpreted by both SST
Writer and Reader engines.

7. ``WANDataTransport``: Default **sockets**.  If the SST
**DataTransport** parameter is **"WAN**, this string value specifies
the EVPath-level data transport to use for exchanging data.  The value
must be a data transport known to EVPath, such as **"sockets"**,
**"enet"**, or **"ib"**.  Generally both the reader and writer should
be using the same EVPath-level data transport.  This value is
interpreted by both SST Writer and Reader engines.

8. ``ControlTransport``: Default **tcp**.  This string value specifies
the underlying network communication mechanism to use for performing
control operations in SST.  SST can be configured to standard TCP
sockets, which are very reliable and efficient, but which are limited
in their scalability.  Alternatively, SST can use a reliable UDP
protocol, that is more scalable, but as of ADIOS2 Release 2.4.0 still
suffers from some reliability problems.  (**sockets** is accepted as
equivalent to **tcp** and **udp**, **rudp**, and **enet** are
equivalent to **scalable**.  Generally both the reader and writer
should be using the same control transport.  This value is interpreted
by both SST Writer and Reader engines.

9. ``NetworkInterface``: Default **NULL**.  In situations in which
there are multiple possible network interfaces available to SST, this
string value specifies which should be used to generate SST's contact
information for writers.  Generally this should *NOT* be specified
except for narrow sets of circumstances.  It has no effect if
specified on Reader engines.  If specified, the string value should
correspond to a name of a network interface, such as are listed by
commands like "netstat -i".  For example, on most Unix systems,
setting the NetworkInterface parameter to "lo" (or possibly "lo0")
will result in SST generating contact information that uses the
network address associated with the loopback interface (127.0.0.1).
This value is interpreted by only by the SST Writer engine.

10. ``ControlInterface``: Default **NULL**.  This value is similar to the
NetworkInterface parameter, but only applies to the SST layer which does
messaging for control (open, close, flow and timestep management, but not
actual data transfer).  Generally the NetworkInterface parameter can be used
to control this, but that also aplies to the Data Plane.  Use
ControlInterface in the event of conflicting specifications.

11. ``DataInterface``: Default **NULL**.  This value is similar to the
NetworkInterface parameter, but only applies to the SST layer which does
messaging for data transfer, not control (open, close, flow and timestep
management).  Generally the NetworkInterface parameter can be used to
control this, but that also aplies to the Control Plane.  Use DataInterface
in the event of conflicting specifications.  In the case of the RDMA data
plane, this parameter controls the libfabric interface choice.

12. ``FirstTimestepPrecious``: Default **FALSE**.
FirstTimestepPrecious is a boolean parameter that affects the queueing
of the first timestep presented to the SST Writer engine. If
FirstTimestepPrecious is **TRUE**, then the first timestep is
effectively never removed from the output queue and will be presented
as a first timestep to any reader that joins at a later time.  This
can be used to convey run parameters or other information that every
reader may need despite joining later in a data stream.  Note that
this queued first timestep does count against the QueueLimit parameter
above, so if a QueueLimit is specified, it should be a value larger
than 1.  Further note while specifying this parameter guarantees that
the preserved first timestep will be made available to new readers,
other reader-side operations (like requesting the LatestAvailable
timestep in Engine parameters) might still cause the timestep to be skipped.
This value is interpreted by only by the SST Writer engine.

13. ``AlwaysProvideLatestTimestep``: Default **FALSE**.
AlwaysProvideLatestTimestep is a boolean parameter that affects what
of the available timesteps will be provided to the reader engine.  If
AlwaysProvideLatestTimestep is **TRUE**, then if there are multiple
timesteps available to the reader, older timesteps will be skipped and
the reader will see only the newest available upon BeginStep.
This value is interpreted by only by the SST Reader engine.

14. ``OpenTimeoutSecs``: Default **60**.  OpenTimeoutSecs is an integer
parameter that specifies the number of seconds SST is to wait for a peer
connection on Open().  Currently this is only implemented on the Reader side
of SST, and is a timeout for locating the contact information file created
by Writer-side Open, not for completing the entire Open() handshake.
Currently value is interpreted by only by the SST Reader engine.

15. ``SpeculativePreloadMode``: Default **AUTO**.  In some
circumstances, SST eagerly sends all data from writers to every
readers without first waiting for read requests.  Generally this
improves performance if every reader needs all the data, but can be
very detrimental otherwise.  The value **AUTO** for this engine
parameter instructs SST to apply its own heuristic for determining if
data should be eagerly sent.  The value **OFF** disables this feature
and the value **ON** causes eager sending regardless of heuristic.
Currently SST's heuristic is simple.  If the size of the reader cohort
is less than or equal to the value of the ``SpecAutoNodeThreshold``
engine parameter (Default value 1), eager sending is initiated.
Currently value is interpreted by only by the SST Reader engine.

16.  ``SpecAutoNodeThreshold``:  Default **1**.  If the size of the
reader cohort is less than or equal to this value *and* the
``SpeculativePreloadMode`` parameter is **AUTO**, SST will initiate
eager data sending of all data from each writer to all readers.
Currently value is interpreted by only by the SST Reader engine.


17. ``StepDistributionMode``: Default **"AllToAll"**.  This value
controls how steps are distributed, particularly when there are
multiple readers.  By default, the value is **"AllToAll"**, which
means that all timesteps are to be delivered to all readers (subject
to discard rules, etc.).  In other distribution modes, this is not the
case.  For example, in **"RoundRobin"**, each step is delivered
only to a single reader, determined in a round-robin fashion based
upon the number or readers who have opened the stream at the time the
step is submitted.  In **"OnDemand"** each step is delivered to a
single reader, but only upon request (with a request being initiated
by the reader doing BeginStep()).  Normal reader-side rules (like
BeginStep timeouts) and writer-side rules (like queue limit behavior) apply.

+-----------------------------+---------------------+----------------------------------------------------+
| **Key**                     | **Value Format**    | **Default** and Examples                           |
+-----------------------------+---------------------+----------------------------------------------------+
| RendezvousReaderCount       | integer             | **1**                                              |
| RegistrationMethod          | string              | **File**, Screen                                   |
| QueueLimit                  | integer             | **0** (no queue limits)                            |
| QueueFullPolicy             | string              | **Block**, Discard                                 |
| ReserveQueueLimit           | integer             | **0** (no queue limits)                            |
| DataTransport               | string              | **default varies by platform**, UCX, MPI, RDMA, WAN|
| WANDataTransport            | string              | **sockets**, enet, ib                              |
| ControlTransport            | string              | **TCP**, Scalable                                  |
| MarshalMethod               | string              | **BP5**, BP, FFS                                   |
| NetworkInterface            | string              | **NULL**                                           |
| ControlInterface            | string              | **NULL**                                           |
| DataInterface               | string              | **NULL**                                           |
| FirstTimestepPrecious       | boolean             | **FALSE**, true, no, yes                           |
| AlwaysProvideLatestTimestep | boolean             | **FALSE**, true, no, yes                           |
| OpenTimeoutSecs             | integer             | **60**                                             |
| SpeculativePreloadMode      | string              | **AUTO**, ON, OFF                                  |
| SpecAutoNodeThreshold       | integer             | **1**                                              |
+-----------------------------+---------------------+----------------------------------------------------+

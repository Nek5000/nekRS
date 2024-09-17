**********
DataSpaces
**********

The DataSpaces engine for ADIOS2 is experimental. DataSpaces is an asynchronous I/O transfer method within ADIOS that enables 
low-overhead, high-throughput data extraction from a running simulation. 
DataSpaces is designed for use in HPC environments and can take advantage of RDMA
network interconnects to speed the transfer of data between communicating
HPC applications.  DataSpaces supports full MxN data distribution, where the number 
of reader ranks can differ from the number of
writer ranks. In addition, this engine supports multiple reader and writer applications, which
must be distinguished by unique values of ``AppID`` for different applications. It can be set
in the xml config file with tag ``<parameter key="AppID" value="2"/>``. The value should be unique 
for each applications or clients.

To use this engine, you can either specify it in your xml config file, with
tag ``<engine type=DATASPACES>`` or, set it in client code. For example, here is
how to create an DataSpaces reader:

.. code-block:: c++

    adios2::IO dspacesIO = adios.DeclareIO("SomeName");
    dspacesIO.SetEngine("DATASPACES");
    adios2::Engine dspacesReader = dspacesIO.Open(filename, adios2::Mode::Read);

and a sample code for DataSpaces writer is:

.. code-block:: c++

    adios2::IO dspacesIO = adios.DeclareIO("SomeName");
    dspacesIO.SetEngine("DATASPACES");
    adios2::Engine dspacesWriter = dspacesIO.Open(filename, adios2::Mode::Write);

To make use of the DataSpaces engine, an application job needs to also run the dataspaces_server
component together with the application. The server should be configured and started 
before the application as a separate job in the system. For example:

``aprun -n $SPROC ./dataspaces_server -s $SPROC &> log.server &``


The variable ``$SPROC`` represents the number of server instances to run. The ``&`` character 
at the end of the line would place the ``aprun`` command in the background, and will 
allow the job script to continue and run the other applications. The server processes 
produce a configuration file, i.e., ``conf.0`` that is used by the application  
to connect to the servers. This file contains identifying information of the 
master server, which coordinates the client registration 
and discovery process. The job script should wait for the servers to start-up and 
produce the ``conf.0`` file before starting the client application processes.

The server also needs a user configuration read from a text file called ``dataspaces.conf``. 
How many output timesteps of the same dataset (called versions) should be kept in the server's memory 
and served to readers should be specified in the file. If this file does not exist in the current directory, 
the server will assume default values (only 1 timestep stored).
.. code-block::

    ## Config file for DataSpaces
    max_versions = 5
    lock_type = 3

The dataspaces_server module is a stand-alone service that runs independently of a simulation 
on a set of dedicated nodes in the staging area. It transfers data from the application through RDMA,  
and can save it to local storage system, e.g., the Lustre file system, stream it to 
remote sites, e.g., auxilliary clusters, or serve it directly from the staging area to 
other applications. One instance of the dataspaces_server can service multiple applications 
in parallel. Further, the server can run in cooperative mode (i.e., multiple 
instances of the server cooperate to service the application in parallel and to balance 
load). The dataspaces_server receives notification messages from the transport method, schedules 
the requests, and initiates the data transfers  in parallel. The 
server schedules and prioritizes the data transfers while the simulation is computing 
in order to overlap data transfers with computations, to maximize data throughput, 
and to minimize the overhead on the application.

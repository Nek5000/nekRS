*********************************************
sst_conn_tool : SST network connectivity tool
*********************************************

The ``sst_conn_tool`` utility exposes some aspects of SST network
connectivity parameters and activity in order to allow debugging of
SST connections. 

In its simplest usage, it just lets you test an SST connection
(between two runs of the program) and tells you the network
information its trying, I.E. what IP address and port it determined to
use for listening, and if it’s connecting somewhere what those
parameters are.  For example, you’d first run sst_conn_tool in one
window and its output would look like this: 

.. code-block:: bash

    bash-3.2$ bin/sst_conn_tool

            Sst writer is listening for TCP/IP connection at IP 192.168.1.17, port 26051

            Sst connection tool waiting for connection…

To try to connect from another window, you run sst_conn_tool with the -c or —connect option:

.. code-block:: bash

    bash-3.2$ bin/sst_conn_tool -c

            Sst reader at IP 192.168.1.17, listening UDP port 26050


            Attempting TCP/IP connection to writer at IP 192.168.1.17, port 26051

    Connection success, all is well!
    bash-3.2$ 

Here, it has found the contact information file, tried and succeeded in making the connection and has indicated that all is well.  In the first window, we get a similar message about the success of the connection.

In the event that there is trouble with the connection, there is a
“-i” or “—info” option that will provide additional information about
the network configuration options.  For example:

.. code-block:: bash

    bash-3.2$ bin/sst_conn_tool -i

            ADIOS2_IP_CONFIG best guess hostname is "sandy.local"
            ADIOS2_IP_CONFIG Possible interface lo0 : IPV4 127.0.0.1
            ADIOS2_IP_CONFIG Possible interface en0 : IPV4 192.168.1.56
            ADIOS2_IP_CONFIG Possible interface en5 : IPV4 192.168.1.17
            ADIOS2_IP_CONFIG best guess IP is "192.168.1.17"
            ADIOS2_IP_CONFIG default port range is "any"

            The following environment variables can impact ADIOS2_IP_CONFIG operation:
                    ADIOS2_IP               -  Publish the specified IP address for contact
                    ADIOS2_HOSTNAME         -  Publish the specified hostname for contact
                    ADIOS2_USE_HOSTNAME     -  Publish a hostname preferentially over IP address
                    ADIOS2_INTERFACE        -  Use the IP address associated with the specified network interface
                    ADIOS2_PORT_RANGE       -  Use a port within the specified range "low:high", 
					       or specify "any" to let the OS choose

            Sst writer is listening for TCP/IP connection at IP 192.168.1.17, port 26048

            Sst connection tool waiting for connection...

Full options for sst_conn_tool:

Operational Modes:


* ``-l`` ``--listen``
    Display connection parameters and wait for an SST connection (default)

* ``-c`` ``--connect``
  Attempt a connection to an already-running instance of sst_conn_tool


Additional Options:

* ``-i`` ``--info``
   Display additional networking information for this host

* ``-f`` ``--file``
   Use file-based contact info sharing (default).  The SST contact
   file is created in the current directory

* ``-s`` ``--screen``
   Use screen-based contact info sharing, SST contact info is displayed/entered via terminal

* ``-h`` ``--help``
    Display this message usage and options




*********
ADIOS
*********

The ``adios2::ADIOS`` component is the initial contact point between an application and the ADIOS2 library.
Applications can be classified as MPI and non-MPI based.
We start by focusing on MPI applications as their non-MPI equivalent just removes the MPI communicator.

.. code-block:: c++

    /** ADIOS class factory of IO class objects */
    adios2::ADIOS adios("config.xml", MPI_COMM_WORLD);

This component is created by passing :

    1. **Runtime config file** (optional): ADIOS2 xml runtime config file, see :ref:`Runtime Configuration Files`.

    2. **MPI communicator** : which determines the scope of the ADIOS library components in an application.

``adios2::ADIOS`` objects can be created in MPI and non-MPI (serial) mode.
Optionally, a runtime configuration file can be passed to the constructor indicating the full file path, name and extension.

**Constructors for MPI applications**

.. code-block:: c++

    /** Constructors */

    // version that accepts an optional runtime adios2 config file
    adios2::ADIOS(const std::string configFile,
                  MPI_COMM mpiComm = MPI_COMM_SELF);

    adios2::ADIOS(MPI_COMM mpiComm = MPI_COMM_SELF);

    /** Examples */
    adios2::ADIOS adios(MPI_COMM_WORLD);
    adios2::ADIOS adios("config.xml", MPI_COMM_WORLD);

**Constructors for non-MPI (serial) applications**

.. code-block:: c++

    /** Constructors */
    adios2::ADIOS(const std::string configFile);

    adios2::ADIOS();

    /** Examples */
    adios2::ADIOS adios("config.xml");
    adios2::ADIOS adios; // Do not use () for empty constructor.


**Factory of IO components**: Multiple IO components (IO tasks) can be created from within the scope of an ADIOS object by calling the ``DeclareIO`` function:

.. code-block:: c++

    /** Signature */
    adios2::IO ADIOS::DeclareIO(const std::string ioName);

    /** Examples */
    adios2::IO bpWriter = adios.DeclareIO("BPWriter");
    adios2::IO bpReader = adios.DeclareIO("BPReader");

This function returns a reference to an existing IO class object that lives inside the ADIOS object that created it.
The ``ioName`` string must be unique; declaring two IO objects with the same name will throw an exception.
IO names are used to identify IO components in the runtime configuration file, :ref:`Runtime Configuration Files`.

As shown in the diagram below, each resulting IO object is self-managed and independent, thus providing an adaptable way to perform different kinds of I/O operations.
Users must be careful not to create conflicts between system level unique I/O identifiers: file names, IP address and port, MPI Send/Receive message rank and tag, etc.

.. blockdiag::

    blockdiag {
        default_fontsize = 18;
        default_shape = roundedbox;
        default_linecolor = blue;
        span_width = 150;

        ADIOS -> IO_1, B, IO_N[label = "DeclareIO", fontsize = 13];
        B[shape = "dots"];
        ADIOS -> B[style = "none"];
    }

.. tip::

    The ADIOS component is the only one whose memory is owned by the application.
    Thus applications must decide on its scope.
    Any other component of the ADIOS2 API refers to a component that lives inside the ADIOS component(e.g. IO, Operator) or indirectly in the IO component(Variable, Engine)

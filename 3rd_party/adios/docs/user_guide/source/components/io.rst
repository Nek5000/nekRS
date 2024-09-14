**
IO
**

The ``IO`` component is the connection between how applications set up their input/output options by selecting an ``Engine`` and its specific parameters, subscribing variables to data, and setting supported transport modes to a particular ``Engine``.
Think of ``IO`` as a control panel for all the user-defined parameters that applications would like to fine tune.
None of the ``IO`` operations are heavyweight until the ``Open`` function that generates an ``Engine`` is called.
Its API allows

* generation of ``Variable`` and ``Attribute`` components containing information about the data in the input output process
* setting ``Engine``-specific parameters and adding supported modes of transport
* generation of ``Engine`` objects to execute the actual IO tasks.

.. note::
   If two different engine types are needed (*e.g.* ``BPFile``, ``SST``), you must define two ``IO`` objects.
   Also, at reading always define separate IOs to avoid ``Variable`` name clashes.


.. blockdiag::

   diagram {
      default_fontsize = 17;
      default_shape = roundedbox;
      default_linecolor = blue;
      span_width = 220;

      IO -> Var_1, B, Var_N [label = "DefineVariable<T>",fontsize = 13];
      B [shape = "dots"];
      IO -> B [style = "none"];

      IO -> Att_1, C, Att_N [label = "DefineAttribute<T>",fontsize = 13];
      C [shape = "dots"];
      IO -> C [style = "none"];

      IO -> Transport_1, D, Transport_N [label = "AddTransport",fontsize = 13];
      D [shape = "dots"];
      IO -> D [style = "none"];

      IO -> Engine_1, E, Engine_N [label = "Open",fontsize = 13];
      E [shape = "dots"];
      IO -> E [style = "none"];
   }


Setting a Particular Engine and its Parameters
----------------------------------------------

Engines execute the heavy operations in ADIOS2.
Each ``IO`` may select a type of ``Engine`` through the ``SetEngine`` function.
If ``SetEngine`` is not called, then the ``BPFile`` engine is used.

.. code-block:: c++

   /** Signature */
   void adios2::IO::SetEngine( const std::string engineType );

   /** Example */
   bpIO.SetEngine("BPFile");

Each ``Engine`` allows the user to fine tune execution of buffering and output tasks via parameters passed to the ``IO`` object.
These parameters are then propagated to the ``Engine``.
For a list of parameters allowed by each engine see :ref:`Available Engines`.

.. note::

   ``adios2::Params`` is an alias to ``std::map<std::string,std::string>`` to pass parameters as key-value string pairs, which can be initialized with curly-brace initializer lists.

.. code-block:: c++

    /** Signature */
    /** Passing several parameters at once */
    void SetParameters(const adios2:Params& parameters);
    /** Passing one parameter key-value pair at a time */
    void SetParameter(const std::string key, const std::string value);

    /** Examples */
    io.SetParameters( { {"Threads", "4"},
                        {"ProfilingUnits", "Milliseconds"},
                        {"MaxBufferSize","2Gb"},
                        {"BufferGrowthFactor", "1.5" }
                        {"FlushStepsCount", "5" }
                      } );
    io.SetParameter( "Threads", "4" );


Adding Supported Transports with Parameters
-------------------------------------------

The ``AddTransport`` function allows the user to specify how data is moved through the system, *e.g.* RDMA, wide-area networks, or files.
It returns an ``unsigned int`` handler for each transport that can be used with the ``Engine::Close`` function at different times.
``AddTransport`` must provide library specific settings that the low-level system library interface allows.

.. code-block:: c++

    /** Signature */
    unsigned int AddTransport( const std::string transportType,
                               const adios2::Params& parameters );

    /** Examples */
    const unsigned int file1 = io.AddTransport( "File",
                                                { {"Library", "fstream"},
                                                  {"Name","file1.bp" }
                                                } );

    const unsigned int file2 = io.AddTransport( "File",
                                                { {"Library", "POSIX"},
                                                  {"Name","file2.bp" }
                                                } );

    const unsigned int wan = io.AddTransport( "WAN",
                                              { {"Library", "Zmq"},
                                                {"IP","127.0.0.1" },
                                                {"Port","80"}
                                              } );


Defining, Inquiring and Removing Variables and Attributes
---------------------------------------------------------

The template functions ``DefineVariable<T>`` allows subscribing to data into ADIOS2 by returning a reference to a ``Variable`` class object whose scope is the same as the ``IO`` object that created it.
The user must provide a unique name, the dimensions: MPI global: shape, MPI local: start and offset, optionally a flag indicating that dimensions are known to be constant, and a data pointer if defined in the application.
Note: data is not passed at this stage.
This is done by the ``Engine`` functions ``Put`` and ``Get`` for Variables.
See the :ref:`Variable` section for supported types and shapes.

.. tip::
   ``adios2::Dims`` is an alias to ``std::vector<std::size_t>``, while ``adios2::ConstantDims`` is an alias to bool ``true``. Use them for code clarity.

.. code-block:: c++

    /** Signature */
    adios2::Variable<T>
        DefineVariable<T>(const std::string name,
                          const adios2::Dims &shape = {}, // Shape of global object
                          const adios2::Dims &start = {}, // Where to begin writing
                          const adios2::Dims &count = {}, // Where to end writing
                          const bool constantDims = false);

    /** Example */
    /** global array of floats with constant dimensions */
    adios2::Variable<float> varFloats =
        io.DefineVariable<float>("bpFloats",
                                 {size * Nx},
                                 {rank * Nx},
                                 {Nx},
                                 adios2::ConstantDims);

Attributes are extra-information associated with the current ``IO`` object.
The function ``DefineAttribute<T>`` allows for defining single value and array attributes.
Keep in mind that Attributes apply to all Engines created by the ``IO`` object and, unlike Variables which are passed to each ``Engine`` explicitly, their definition contains their actual data.

.. code-block:: c++

    /** Signatures */

    /** Single value */
    adios2::Attribute<T> DefineAttribute(const std::string &name,
                                  const T &value);

    /** Arrays */
    adios2::Attribute<T> DefineAttribute(const std::string &name,
                                  const T *array,
                                  const size_t elements);

In situations in which a variable and attribute has been previously defined:
1) a variable/attribute reference goes out of scope, or 2) when reading from an incoming stream, the ``IO`` can inquire about the status of variables and attributes.
If the inquired variable/attribute is not found, then the overloaded ``bool()`` operator of returns ``false``.

.. code-block:: c++

    /** Signature */
    adios2::Variable<T> InquireVariable<T>(const std::string &name) noexcept;
    adios2::Attribute<T> InquireAttribute<T>(const std::string &name) noexcept;

    /** Example */
    adios2::Variable<float> varPressure = io.InquireVariable<float>("pressure");
    if( varPressure ) // it exists
    {
      ...
    }


.. note::
   ``adios2::Variable`` overloads ``operator bool()`` so that we can check for invalid states (e.g. variables haven't arrived in a stream, weren't previously defined, or weren't written in a file).

.. caution::

   Since ``InquireVariable`` and ``InquireAttribute`` are template functions, both the name and type must match the data you are looking for.


Opening an Engine
-----------------

The ``IO::Open`` function creates a new derived object of the abstract ``Engine`` class and returns a reference handler to the user.
A particular ``Engine`` type is set to the current ``IO`` component with the ``IO::SetEngine`` function.
Engine polymorphism is handled internally by the ``IO`` class, which allows subclassing future derived ``Engine`` types without changing the basic API.

``Engine`` objects are created in various modes.
The available modes are ``adios2::Mode::Read``, ``adios2::Mode::ReadRandomAccess``, ``adios2::Mode::Write``, ``adios2::Mode::Append``, ``adios2::Mode::Sync``, ``adios2::Mode::Deferred``, and ``adios2::Mode::Undefined``.


.. code-block:: c++

    /** Signatures */
    /** Provide a new MPI communicator other than from ADIOS->IO->Engine */
    adios2::Engine adios2::IO::Open(const std::string &name,
                                    const adios2::Mode mode,
                                    MPI_Comm mpiComm );

    /** Reuse the MPI communicator from ADIOS->IO->Engine \n or non-MPI serial mode */
    adios2::Engine adios2::IO::Open(const std::string &name,
                                    const adios2::Mode mode);


    /** Examples */

    /** Engine derived class, spawned to start Write operations */
    adios2::Engine bpWriter = io.Open("myVector.bp", adios2::Mode::Write);

    /** Engine derived class, spawned to start Read operations on rank 0 */
    if( rank == 0 )
    {
        adios2::Engine bpReader = io.Open("myVector.bp",
                                           adios2::Mode::Read,
                                           MPI_COMM_SELF);
    }

.. caution::

   Always pass ``MPI_COMM_SELF`` if an ``Engine`` lives in only one MPI process.
   ``Open`` and ``Close`` are collective operations.

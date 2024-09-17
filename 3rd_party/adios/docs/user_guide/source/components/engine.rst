******
Engine
******

.. _sec:basics_interface_components_engine:

The Engine abstraction component serves as the base interface to the actual IO systems executing the heavy-load tasks performed when producing and consuming data.

Engine functionality works around two concepts:

1. Variables are published (``Put``) and consumed (``Get``) in "steps" in either "File" random-access (all steps are available) or "Streaming" (steps are available as they are produced in a step-by-step fashion).
2. Variables are published (``Put``) and consumed (``Get``) using a "sync" or "deferred" (lazy evaluation) policy.

.. caution::

   The ADIOS2 "step" is a logical abstraction that means different things depending on the application context.
   Examples: "time step", "iteration step", "inner loop step", or "interpolation step", "variable section", etc.
   It only indicates how the variables were passed into ADIOS2 (e.g. I/O steps) without the user having to index this information on their own.

.. tip::
   
   Publishing and consuming data is a round-trip in ADIOS2.
   ``Put`` and ``Get`` APIs for write/append and read modes aim to be "symmetric", reusing functions, objects, and semantics as much as possible.

The rest of the section explains the important concepts.

BeginStep
---------

   Begins a logical step and return the status (via an enum) of the stream to be read/written.
   In streaming engines ``BeginStep`` is where the receiver tries to acquire a new step in the reading process.
   The full signature allows for a mode and timeout parameters.
   See :ref:`Supported Engines` for more information on what engine allows.
   A simplified signature allows each engine to pick reasonable defaults.

.. code-block:: c++

   // Full signature
   StepStatus BeginStep(const StepMode mode,
                        const float timeoutSeconds = -1.f); 

   // Simplified signature
   StepStatus BeginStep();

EndStep
-------
        
   Ends logical step, flush to transports depending on IO parameters and engine default behavior.


.. tip::
   
   To write portable code for a step-by-step access across ADIOS2 engines (file and streaming engines) use ``BeginStep`` and ``EndStep``.

.. danger:: 
   
   Accessing random steps in read mode (e.g. ``Variable<T>::SetStepSelection`` in file engines) will create a conflict with ``BeginStep`` and ``EndStep`` and will throw an exception.
   In file engines, data is either consumed in a random-access or step-by-step mode, but not both.


Close
-----

   Close current engine and underlying transports.
   An ``Engine`` object can't be used after this call.


Put: modes and memory contracts
-------------------------------

``Put`` publishes data in ADIOS2.
It is unavailable unless the ``Engine`` is created in ``Write`` or ``Append`` mode.

The most common signature is the one that passes a ``Variable<T>``
object for the metadata, a ``const`` piece of contiguous memory for
the data, and a mode for either ``Deferred`` (data may be collected at
Put() or not until EndStep/PerformPuts/Close) or ``Sync`` (data is reusable immediately).
This is the most common use case in applications.

1. Deferred (default) or Sync mode, data is contiguous memory 

   .. code-block:: c++

      void Put(Variable<T> variable, const T* data, const adios2::Mode = adios2::Mode::Deferred);

ADIOS2 Engines also provide direct access to their buffer memory.
``Variable<T>::Span`` is based on a subset of the upcoming `C++20 std::span <https://en.cppreference.com/w/cpp/container/span>`_, which is a non-owning reference to a block of contiguous memory.
Spans act as a 1D container meant to be filled out by the application.
They provide the standard API of an STL container, providing ``begin()`` and ``end()`` iterators, ``operator[]`` and ``at()``, as well as ``data()`` and ``size()``.

``Variable<T>::Span`` is helpful in situations in which temporaries are needed to create contiguous pieces of memory from non-contiguous pieces (*e.g.* tables, arrays without ghost-cells), or just to save memory as the returned ``Variable<T>::Span`` can be used for computation, thus avoiding an extra copy from user memory into the ADIOS2 buffer.
``Variable<T>::Span`` combines a hybrid ``Sync`` and ``Deferred`` mode, in which the initial value and memory allocations are ``Sync``, while data population and metadata collection are done at EndStep/PerformPuts/Close.
Memory contracts are explained later in this chapter followed by examples.

The following ``Variable<T>::Span`` signatures are available:

2. Return a span setting a default ``T()`` value into a default buffer
 
   .. code-block:: c++
   
      Variable<T>::Span Put(Variable<T> variable);
      
3. Return a span setting an initial fill value into a certain buffer.
If span is not returned then the ``fillValue`` is fixed for that block.

   .. code-block:: c++

      Variable<T>::Span Put(Variable<T> variable, const size_t bufferID, const T fillValue);


In summary, the following are the current Put signatures for publishing data in ADIOS 2:

1. ``Deferred`` (default) or ``Sync`` mode, data is contiguous memory put in an ADIOS2 buffer.

   .. code-block:: c++

      void Put(Variable<T> variable, const T* data, const adios2::Mode = adios2::Mode::Deferred);
   
2. Return a span setting a default ``T()`` value into a default ADIOS2 buffer.
If span is not returned then the default ``T()`` is fixed for that block (e.g. zeros).
 
   .. code-block:: c++
   
      Variable<T>::Span Put(Variable<T> variable);
   
3. Return a span setting an initial fill value into a certain buffer.
If span is not returned then the ``fillValue`` is fixed for that block.

   .. code-block:: c++

      Variable<T>::Span Put(Variable<T> variable, const size_t bufferID, const T fillValue);


The following table summarizes the memory contracts required by ADIOS2 engines between ``Put`` signatures and the data memory coming from an application:

+----------+-------------+----------------------------------------------------+
| Put      | Data Memory | Contract                                           |
+----------+-------------+----------------------------------------------------+
|          | Pointer     | do not modify until PerformPuts/EndStep/Close      |
| Deferred |             |                                                    |
|          | Contents    | consumed at Put or PerformPuts/EndStep/Close       |
+----------+-------------+----------------------------------------------------+
|          | Pointer     | modify after Put                                   |
| Sync     |             |                                                    |
|          | Contents    | consumed at Put                                    |
+----------+-------------+----------------------------------------------------+
|          | Pointer     | modified by new Spans, updated span iterators/data |
| Span     |             |                                                    |
|          | Contents    | consumed at PerformPuts/EndStep/Close              |
+----------+-------------+----------------------------------------------------+


.. note::

   In Fortran (array) and Python (numpy array) avoid operations that modify the internal structure of an array (size) to preserve the address. 
   
   
Each ``Engine`` will give a concrete meaning to  each functions signatures, but all of them must follow the same memory contracts to the "data pointer": the memory address itself, and the "data contents": memory bits (values).
   
1. **Put in Deferred or lazy evaluation mode (default)**: this is the preferred mode as it allows ``Put`` calls to be "grouped" before potential data transport at the first encounter of ``PerformPuts``, ``EndStep`` or ``Close``.
   
     .. code-block:: c++
         
         Put(variable, data);
         Put(variable, data, adios2::Mode::Deferred);
         

   Deferred memory contracts: 
      
   - "data pointer" do not modify (e.g. resize) until first call to ``PerformPuts``, ``EndStep`` or ``Close``.
      
   - "data contents" may be consumed immediately or at first call to
     ``PerformPuts``, ``EndStep`` or ``Close``.  Do not modify data contents after Put.


   Usage:

      .. code-block:: c++
         
         // recommended use: 
         // set "data pointer" and "data contents"
         // before Put
         data[0] = 10;  
         
         // Puts data pointer into adios2 engine
         // associated with current variable metadata
         engine.Put(variable, data);
         
         // Modifying data after Put(Deferred) may result in different
	 // results with different engines
         // Any resize of data after Put(Deferred) may result in
	 // memory corruption or segmentation faults
         data[1] = 10; 
         
         // "data contents" must not have been changed
         // "data pointer" must be the same as in Put
         engine.EndStep();   
         //engine.PerformPuts();  
         //engine.Close();
         
         // now data pointer can be reused or modified
        
   .. tip::

      It's recommended practice to set all data contents before ``Put`` in deferred mode to minimize the risk of modifying the data pointer (not just the contents) before PerformPuts/EndStep/Close.


2.  **Put in Sync mode**: this is the special case, data pointer becomes reusable right after ``Put``.
Only use it if absolutely necessary (*e.g.* memory bound application or out of scope data, temporary).
   
      .. code-block:: c++
         
         Put(variable, *data, adios2::Mode::Sync);
         

   Sync memory contracts:
      
   - "data pointer" and "data contents" can be modified after this call.
   
   
   Usage:

      .. code-block:: c++
         
         // set "data pointer" and "data contents"
         // before Put in Sync mode
         data[0] = 10;  
         
         // Puts data pointer into adios2 engine
         // associated with current variable metadata
         engine.Put(variable, data, adios2::Mode::Sync);
         
         // data pointer and contents can be reused
         // in application 
   
   
3. **Put returning a Span**: signature that allows access to adios2 internal buffer. 

   Use cases: 
   
   -  population from non-contiguous memory structures
   -  memory-bound applications 


   Limitations:
   
   -  does not allow operations (compression)
   -  must keep engine and variables within scope of span usage 
     


   Span memory contracts: 
      
   - "data pointer" provided by the engine and returned by ``span.data()``, might change with the generation of a new span. It follows iterator invalidation rules from std::vector. Use `span.data()` or iterators, `span.begin()`, `span.end()` to keep an updated data pointer.
      
   - span "data contents" are published at the first call to ``PerformPuts``, ``EndStep`` or ``Close``


   Usage:

       .. code-block:: c++
         
         // return a span into a block of memory
         // set memory to default T()
         adios2::Variable<int32_t>::Span span1 = Put(var1);
         
         // just like with std::vector::data()
         // iterator invalidation rules
         // dataPtr might become invalid
         // always use span1.data() directly
         T* dataPtr = span1.data();
         
         // set memory value to -1 in buffer 0
         adios2::Variable<float>::Span span2 = Put(var2, 0, -1);

         // not returning a span just sets a constant value 
         Put(var3);
         Put(var4, 0, 2);
         
         // fill span1
         span1[0] = 0;
         span1[1] = 1;
         span1[2] = 2;
         
         // fill span2
         span2[1] = 1;
         span2[2] = 2;
         
         // here collect all spans
         // they become invalid
         engine.EndStep();
         //engine.PerformPuts();  
         //engine.Close();
         
         // var1 = { 0, 1, 2 };
         // var2 = { -1., 1., 2.};
         // var3 = { 0, 0, 0};
         // var4 = { 2, 2, 2};


The ``data`` fed to the ``Put`` function is assumed to be allocated on the Host (default mode). In order to use data allocated on the device, the memory space of the variable needs to be set to Cuda.

     .. code-block:: c++

         variable.SetMemorySpace(adios2::MemorySpace::CUDA);
         engine.Put(variable, gpuData, mode);

.. note::

   Only CUDA allocated buffers are supported for device data.
   Only the BP4 and BP5 engines are capable of receiving device allocated buffers.


PerformPuts
-----------

   Executes all pending ``Put`` calls in deferred mode and collects
   span data.  Specifically this call copies Put(Deferred) data into
   internal ADIOS buffers, as if Put(Sync) had been used instead.

.. note::

   This call allows the reuse of user buffers, but may negatively
   impact performance on some engines.


PerformDataWrite
----------------

   If supported by the engine, moves data from prior ``Put`` calls to disk

.. note::

   Currently only supported by the BP5 file engine.



Get: modes and memory contracts
-------------------------------

``Get`` is the function for consuming data in ADIOS2.
It is available when an Engine is created using ``Read`` mode at ``IO::Open``.
ADIOS2 ``Put`` and ``Get`` semantics are as symmetric as possible considering that they are opposite operations (*e.g.* ``Put`` passes ``const T*``, while ``Get`` populates a non-const ``T*``).

The ``Get`` signatures are described below.

1. ``Deferred`` (default) or ``Sync`` mode, data is contiguous pre-allocated memory:

   .. code-block:: c++

      Get(Variable<T> variable, const T* data, const adios2::Mode = adios2::Mode::Deferred);


2. In this signature, ``dataV`` is automatically resized by ADIOS2 based on the ``Variable`` selection:

   .. code-block:: c++

      Get(Variable<T> variable, std::vector<T>& dataV, const adios2::Mode = adios2::Mode::Deferred);


The following table summarizes the memory contracts required by ADIOS2 engines between ``Get`` signatures and the pre-allocated (except when using C++11 ``std::vector``) data memory coming from an application:

+----------+-------------+-----------------------------------------------+
| Get      | Data Memory | Contract                                      |
+----------+-------------+-----------------------------------------------+
|          | Pointer     | do not modify until PerformGets/EndStep/Close |
| Deferred |             |                                               |
|          | Contents    | populated at Get or PerformGets/EndStep/Close |
+----------+-------------+-----------------------------------------------+
|          | Pointer     | modify after Get                              |
| Sync     |             |                                               |
|          | Contents    | populated at Get                              |
+----------+-------------+-----------------------------------------------+


1. **Get in Deferred or lazy evaluation mode (default)**: this is the preferred mode as it allows ``Get`` calls to be "grouped" before potential data transport at the first encounter of ``PerformPuts``, ``EndStep`` or ``Close``.
   
     .. code-block:: c++
         
         Get(variable, data);
         Get(variable, data, adios2::Mode::Deferred);
         

   Deferred memory contracts: 
      
   - "data pointer": do not modify (e.g. resize) until first call to ``PerformPuts``, ``EndStep`` or ``Close``.
      
   - "data contents": populated at ``Put``, or at first call to ``PerformPuts``, ``EndStep`` or ``Close``.

   Usage:`

      .. code-block:: c++

         std::vector<double> data;

         // resize memory to expected size 
         data.resize(varBlockSize);
         // valid if all memory is populated 
         // data.reserve(varBlockSize);

         // Gets data pointer to adios2 engine
         // associated with current variable metadata
         engine.Get(variable, data.data() );

         // optionally pass data std::vector 
         // leave resize to adios2
         //engine.Get(variable, data);

         // "data pointer" must be the same as in Get
         engine.EndStep();   
         // "data contents" are now ready
         //engine.PerformPuts();  
         //engine.Close();

         // now data pointer can be reused or modified



2.  **Put in Sync mode**: this is the special case, data pointer becomes reusable right after Put.
Only use it if absolutely necessary (*e.g.* memory bound application or out of scope data, temporary).
   
      .. code-block:: c++
         
         Get(variable, *data, adios2::Mode::Sync);
         

   Sync memory contracts:
      
   - "data pointer" and "data contents" can be modified after this call.
   
   
   Usage:

      .. code-block:: c++
         
         .. code-block:: c++
         
         std::vector<double> data;
         
         // resize memory to expected size 
         data.resize(varBlockSize);
         // valid if all memory is populated 
         // data.reserve(varBlockSize);
         
         // Gets data pointer to adios2 engine
         // associated with current variable metadata
         engine.Get(variable, data.data() );
         
         // "data contents" are ready
         // "data pointer" can be reused by the application

.. note::

   ``Get`` doesn't support returning spans.


PerformGets
-----------

   Executes all pending ``Get`` calls in deferred mode.


Engine usage example
--------------------

The following example illustrates the basic API usage in write mode for data generated at each application step:

.. code-block:: c++

   adios2::Engine engine = io.Open("file.bp", adios2::Mode::Write);

   for( size_t i = 0; i < steps; ++i )
   {
      // ... Application *data generation

      engine.BeginStep(); //next "logical" step for this application

      engine.Put(varT, dataT, adios2::Mode::Sync);
      // dataT memory already consumed by engine
      // Application can modify dataT address and contents
      
      // deferred functions return immediately (lazy evaluation),
      // dataU, dataV and dataW pointers and contents must not be modified
      // until PerformPuts, EndStep or Close.
      // 1st batch
      engine.Put(varU, dataU);
      engine.Put(varV, dataV);
      
      // in this case adios2::Mode::Deferred is redundant,
      // as this is the default option
      engine.Put(varW, dataW, adios2::Mode::Deferred);

      // effectively dataU, dataV, dataW are "deferred"
      // possibly until the first call to PerformPuts, EndStep or Close.
      // Application MUST NOT modify the data pointer (e.g. resize
      // memory) or change data contents.
      engine.PerformPuts();

      // dataU, dataV, dataW pointers/values can now be reused
      
      // ... Application modifies dataU, dataV, dataW 

      //2nd batch
      dataU[0] = 10
      dataV[0] = 10
      dataW[0] = 10 
      engine.Put(varU, dataU);
      engine.Put(varV, dataV);
      engine.Put(varW, dataW);
      // Application MUST NOT modify dataU, dataV and dataW pointers (e.g. resize),
      // Contents should also not be modified after Put() and before
      // PerformPuts() because ADIOS may access the data immediately
      // or not until PerformPuts(), depending upon the engine
      engine.PerformPuts();
      
      // dataU, dataV, dataW pointers/values can now be reused
      
      // Puts a varP block of zeros
      adios2::Variable<double>::Span spanP = Put<double>(varP);
      
      // Not recommended mixing static pointers, 
      // span follows 
      // the same pointer/iterator invalidation  
      // rules as std::vector
      T* p = spanP.data();

      // Puts a varMu block of 1e-6
      adios2::Variable<double>::Span spanMu = Put<double>(varMu, 0, 1e-6);
      
      // p might be invalidated 
      // by a new span, use spanP.data() again
      foo(spanP.data());

      // Puts a varRho block with a constant value of 1.225
      Put<double>(varMu, 0, 1.225);
      
      // it's preferable to start modifying spans 
      // after all of them are created
      foo(spanP.data());
      bar(spanMu.begin(), spanMu.end()); 
      
      
      engine.EndStep();
      // spanP, spanMu are consumed by the library
      // end of current logical step,
      // default behavior: transport data
   }

   engine.Close();
   // engine is unreachable and all data should be transported
   ...

.. tip::

   Prefer default ``Deferred`` (lazy evaluation) functions as they have the potential to group several variables with the trade-off of not being able to reuse the pointers memory space until ``EndStep``, ``PerformPuts``, ``PerformGets``, or ``Close``.
   Only use ``Sync`` if you really have to (*e.g.* reuse memory space from pointer).
   ADIOS2 prefers a step-based IO in which everything is known ahead of time when writing an entire step.


.. danger::
   The default behavior of ADIOS2 ``Put`` and ``Get`` calls IS NOT synchronized, but rather deferred.
   It's actually the opposite of ``MPI_Put`` and more like ``MPI_rPut``.
   Do not assume the data pointer is usable after a ``Put`` and ``Get``, before ``EndStep``, ``Close`` or the corresponding ``PerformPuts``/``PerformGets``.
   Avoid using temporaries, r-values, and out-of-scope variables in ``Deferred`` mode.
   Use ``adios2::Mode::Sync`` in these cases.


Available Engines
-----------------

A particular engine is set within the ``IO`` object that creates it with the ``IO::SetEngine`` function in a case insensitive manner.
If the ``SetEngine`` function is not invoked the default engine is the ``BPFile``.

+-------------------------+---------+---------------------------------------------+
| Application             | Engine  | Description                                 |
+-------------------------+---------+---------------------------------------------+
| File                    | BP5     | DEFAULT write/read ADIOS2 native bp files   |
|                         |         |                                             |
|                         | HDF5    | write/read interoperability with HDF5 files |
+-------------------------+---------+---------------------------------------------+
| Wide-Area-Network (WAN) | DataMan | write/read TCP/IP streams                   |
+-------------------------+---------+---------------------------------------------+
| Staging                 | SST     | write/read to a "staging" area: *e.g.* RDMA |
+-------------------------+---------+---------------------------------------------+


``Engine`` polymorphism has two goals:

1. Each ``Engine`` implements an orthogonal IO scenario targeting a use case (e.g. Files, WAN, InSitu MPI, etc) using a simple, unified API.

2. Allow developers to build their own custom system solution based on their particular requirements in the own playground space.
Reusable toolkit objects are available inside ADIOS2 for common tasks: bp buffering, transport management, transports, etc.

A class that extends ``Engine`` must be thought of as a solution to a range of IO applications.
Each engine must provide a list of supported parameters, set in the IO object creating this engine using ``IO::SetParameters``, and supported transports (and their parameters) in ``IO::AddTransport``.
Each Engine's particular options are documented in :ref:`Supported Engines`.

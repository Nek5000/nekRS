********
Variable
********

An ``adios2::Variable`` is the link between a piece of data coming from an application and its metadata.
This component handles all application variables classified by data type and shape.

Each ``IO`` holds a set of Variables, and each ``Variable`` is identified with a unique name.
They are created using the reference from ``IO::DefineVariable<T>`` or retrieved using the pointer from
``IO::InquireVariable<T>`` functions in :ref:`IO`.

Data Types
----------

Only primitive types are supported in ADIOS2.
Fixed-width types from `<cinttypes> and <cstdint> <https://en.cppreference.com/w/cpp/types/integer>`_  should be
preferred when writing portable code. ADIOS2 maps primitive types to equivalent fixed-width types
(e.g. ``int`` -> ``int32_t``). In C++, acceptable types ``T`` in ``Variable<T>`` along with their preferred fix-width
equivalent in 64-bit platforms are given below:

.. code-block:: c++

   Data types Variables supported by ADIOS2 Variable<T>

   std::string (only used for global and local values, not arrays)
   char                      -> int8_t or uint8_t depending on compiler flags
   signed char               -> int8_t 
   unsigned char             -> uint8_t
   short                     -> int16_t
   unsigned short            -> uint16_t
   int                       -> int32_t
   unsigned int              -> uint32_t 
   long int                  -> int32_t or int64_t (Linux)
   long long int             -> int64_t 
   unsigned long int         -> uint32_t or uint64_t (Linux)
   unsigned long long int    -> uint64_t  
   float                     -> always 32-bit = 4 bytes  
   double                    -> always 64-bit = 8 bytes
   long double               -> platform dependent
   std::complex<float>       -> always  64-bit = 8 bytes = 2 * float
   std::complex<double>      -> always 128-bit = 16 bytes = 2 * double

.. tip::

   It's recommended to be consistent when using types for portability.
   If data is defined as a fixed-width integer, define variables in ADIOS2 using a fixed-width type, *e.g.*  for ``int32_t`` data types use ``DefineVariable<int32_t>``.

.. note::

   C, Fortran APIs: the enum and parameter adios2_type_XXX only provides fixed-width types.
   
.. note::

   Python APIs: use the equivalent fixed-width types from numpy.
   If ``dtype`` is not specified, ADIOS2 handles numpy defaults just fine as long as primitive types are passed.

Shapes
------

ADIOS2 is designed for MPI applications.
Thus different application data shapes must be supported depending on their scope within a particular MPI communicator.
The shape is defined at creation from the ``IO`` object by providing the dimensions: shape, start, count in the
``IO::DefineVariable<T>``. The supported shapes are described below.


1. **Global Single Value**:
Only a name is required for their definition.
These variables are helpful for storing global information, preferably managed by only one MPI process, that may or may
not change over steps: *e.g.* total number of particles, collective norm, number of nodes/cells, etc.

   .. code-block:: c++

      if( rank == 0 )
      {
         adios2::Variable<uint32_t> varNodes = io.DefineVariable<uint32_t>("Nodes");
         adios2::Variable<std::string> varFlag = io.DefineVariable<std::string>("Nodes flag");
         // ...
         engine.Put( varNodes, nodes );
         engine.Put( varFlag, "increased" );
         // ...
      }

   .. note::

      Variables of type ``string`` are defined just like global single values.
      Multidimensional strings are supported for fixed size strings through variables of type ``char``.


2. **Global Array**:
This is the most common shape used for storing data that lives in several MPI processes.
The image below illustrates the definitions of the dimension components in a global array: shape, start, and count.

   .. image:: https://i.imgur.com/MKwNe5e.png
   
   .. warning::

      Be aware of data ordering in your language of choice (row-major or column-major) as depicted in the figure above.
      Data decomposition is done by the application, not by ADIOS2.

   Start and Count local dimensions can be later modified with the ``Variable::SetSelection`` function if it is not a constant dimensions variable.


3. **Local Value**:
Values that are local to the MPI process.
They are defined by passing the ``adios2::LocalValueDim`` enum as follows:

   .. code-block:: c++

      adios2::Variable<int32_t> varProcessID =
            io.DefineVariable<int32_t>("ProcessID", {adios2::LocalValueDim})
      //...
      engine.Put<int32_t>(varProcessID, rank);

These values become visible on the reader as a single merged 1-D
Global Array whose size is determined by the number of writer ranks.

4. **Local Array**:
Arrays that are local to the MPI process.
These are commonly used to write checkpoint-restart data.
Reading, however, needs to be handled differently: each process' array has to be read separately, using ``SetSelection`` per rank.
The size of each process selection should be discovered by the reading application by inquiring per-block size information of the variable, and allocate memory accordingly.

  .. image:: https://i.imgur.com/XLh2TUG.png


.. note::

   Constants are not handled separately from step-varying values in ADIOS2.
   Simply write them only once from one rank.

5. **Joined Array**:
Joined arrays are a variation of the Local Array described above.
Where LocalArrays are only available to the reader via their block
number, JoinedArrays are merged into a single global array whose
global dimensions are determined by the sum of the contributions of
each writer rank.   Specifically:  JoinedArrays are N-dimensional
arrays where one (and only one) specific dimension is the Joined
dimension.  (The other dimensions must be constant and the same across
all contributions.)  When defining a Joined variable, one specifies a
shape parameter that give the dimensionality of the array with the
special constant ``adios2::JoinedDim`` in the dimension to be joined.
Unlike a Global Array definition, the start parameter must be an empty
Dims value.
For example, the definition below defines a 2-D Joined array where the
first dimension is the one along which blocks will be joined and the
2nd dimension is 5.  Here this rank is contributing two rows to this array.

.. code-block:: c++

  auto var = outIO.DefineVariable<double>("table", {adios2::JoinedDim, 5}, {}, {2, 5});

If each of N writer ranks were to declare a variable like this and do
a single Put() in a timestep, the reader-side GlobalArray would have
shape {2*N, 5} and all normal reader-side GlobalArray operations would
be applicable to it.  


.. note::

   JoinedArrays are currently only supported by BP4 and BP5 engines,
   as well as the SST engine with BP5 marshalling.

Global Array Capabilities and Limitations
-----------------------------------------

ADIOS2 is focusing on writing and reading N-dimensional, distributed, global arrays of primitive types. The basic idea
is that, usually, a simulation has such a data structure in memory (distributed across multiple processes) and wants to
dump its content regularly as it progresses. ADIOS2 was designed to:

1. to do this writing and reading as fast as possible
2. to enable reading any subsection of the array

.. image:: https://imgur.com/6nX67yq.png
   :width: 400

The figure above shows a parallel application of 12 processes producing a 2D array. Each process has a 2D array locally
and the output is created by placing them into a 4x3 pattern. A reading application's individual process then can read
any subsection of the entire global array. In the figure, a 6 process application decomposes the array in a 3x2 pattern
and each process reads a 2D array whose content comes from multiple producer processes.

The figure hopefully helps to understand the basic concept but it can be also misleading if it suggests limitations that
are not there. Global Array is simply a boundary in N-dimensional space where processes can place their blocks of data.
In the global space:

1. one process can place multiple blocks

  .. image:: https://imgur.com/Pb1s03h.png
     :width: 400

2. does NOT need to be fully covered by the blocks

  .. image:: https://imgur.com/qJBXYcQ.png
     :width: 400

  * at reading, unfilled positions will not change the allocated memory

3. blocks can overlap

  .. image:: https://imgur.com/GA59lZ2.png
     :width: 300

  * the reader will get values in an overlapping position from one of the block but there is no control over from which
    block

4. each process can put a different size of block, or put multiple blocks of different sizes

5. some process may not contribute anything to the global array

Over multiple output steps

1. the processes CAN change the size (and number) of blocks in the array

  * E.g. atom table: global size is fixed but atoms wander around processes, so their block size is changing

    .. image:: https://imgur.com/DorjG2q.png
       :width: 400

2. the global dimensions CAN change over output steps

  * but then you cannot read multiple steps at once
  * E.g. particle table size changes due to particles disappearing or appearing

    .. image:: https://imgur.com/nkuHeVX.png
       :width: 400


Limitations of the ADIOS global array concept

1. Indexing starts from 0
2. Cyclic data patterns are not supported; only blocks can be written or read
3. If Some blocks may fully or partially fall outside of the global boundary, the reader will not be able to read those
   parts

.. note::

   Technically, the content of the individual blocks is kept in the BP format (but not in HDF5 format) and in staging.
   If you really, really want to retrieve all the blocks, you need to handle this array as a Local Array and read the
   blocks one by one.

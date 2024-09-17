*************************
adios2 classes for Python
*************************

``Stream`` is a high-level class that can perform most of the ADIOS functionality. ``FileReader`` is just a convenience class and is the same as Stream with "rra" (ReadRandomAccess) mode. FileReaders do not work with for loops as opposed to Streams that work step-by-step, rather one can access any step of any variable at will. The other classes, ``Adios``, ``IO``, ``Engine``, ``Variable``, ``Attribute`` and ``Operator`` correspond to the C++ classes. One needs to use them to extend the capabilities of the ``Stream`` class (e.g. using an external XML file for runtime configuration, changing the engine for the run, setting up a compression operator for an output variable, etc.)

.. autoclass:: adios2::Stream
    :members:

.. autoclass:: adios2::FileReader
    :members:

.. autoclass:: adios2::Adios
    :members:

.. autoclass:: adios2::IO
    :members:

.. autoclass:: adios2::Engine
    :members:

.. autoclass:: adios2::Variable
    :members:

.. autoclass:: adios2::Attribute
    :members:

.. autoclass:: adios2::Operator
    :members:

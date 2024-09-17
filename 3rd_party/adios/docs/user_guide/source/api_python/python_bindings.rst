**********************
Python bindings to C++
**********************

.. note::

   The bindings to the C++ functions is the basis of the native Python API described before. It is still accessible to users who used the "Full Python API" pre-2.10. In order to make old scripts working with 2.10 and later versions, change the import line in the python script.

.. code-block:: python

    import adios2.bindings as adios2

The full Python APIs follows very closely the full C++11 API interface. All of its functionality is now in the native API as well, so its use is discouraged for future scripts.

Examples using the Python bindings in the ADIOS2 repository
-----------------------------------------------------------

- Simple file-based examples
    - examples/hello/helloWorld/hello-world-bindings.py
    - examples/hello/bpReader/bpReaderHeatMap2D-bindings.py
    - examples/hello/bpWriter/bpWriter-bindings.py

- Staging examples using staging engines SST and DataMan
    - examples/hello/sstWriter/sstWriter-bindings.py
    - examples/hello/sstReader/sstReader-bindings.py


ADIOS class
--------------
.. autoclass:: adios2.bindings::ADIOS
    :members:

IO class
--------------
.. autoclass:: adios2.bindings::IO
    :members:

Variable class
--------------
.. autoclass:: adios2.bindings::Variable
    :members:

Attribute class
---------------
.. autoclass:: adios2.bindings::Attribute
    :members:

Engine class
--------------
.. autoclass:: adios2.bindings::Engine
    :members:

Operator class
--------------
.. autoclass:: adios2.bindings::Operator
    :members:

Query class
--------------
.. autoclass:: adios2.bindings::Query
    :members:

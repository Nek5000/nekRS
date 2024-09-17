****
Null 
****

The ``Null`` Engine performs no internal work and no I/O. It was created for testing applications that have ADIOS2 output in it by turning off the I/O easily. The runtime difference between a run with the Null engine and another engine tells us the IO overhead of that particular output with that particular engine.  

.. code-block:: c++

    adios2::IO io = adios.DeclareIO("Output");
    io.SetEngine("Null");

or using the XML config file:

.. code-block:: xml

    <adios-config>
        <io name="Output">
            <engine type="Null">
            </engine>
        </io>
    </adios-config>

Although there is a reading engine as well, which will not fail, any variable/attribute inquiry returns `nullptr` and any subsequent Get() calls will throw an exception in C++/Python or return an error in C/Fortran. 

Note that there is also a `Null transport` that can be used by a BP engine instead of the default `File transport`. In that case, the BP engine will perform all internal work including buffering and aggregation but no data will be output at all. A run like this can be used to assess the overhead of the internal work of the BP engine. 

.. code-block:: c++

    adios2::IO io = adios.DeclareIO("Output");
    io.SetEngine("BP5");
    io.AddTransport("Null", {});

or using the XML config file

.. code-block:: xml

    <adios-config>
        <io name="Output">
            <engine type="BP5">
            </engine>
            <transport type="Null">
            </transport>
        </io>
    </adios-config>
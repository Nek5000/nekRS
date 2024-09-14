********************
Inline for zero-copy
********************

The ``Inline`` engine provides in-process communication between the writer and reader, avoiding the copy of data buffers.

This engine is focused on the `N` → `N` case: `N` writers share a process with `N` readers, and the analysis happens 'inline' without writing the data to a file or copying to another buffer.
It is similar to the streaming SST engine, since analysis must happen per step.

To use this engine, you can either add ``<engine type=Inline>`` to your XML config file, or set it in your application code:

.. code-block:: c++

    adios2::IO io = adios.DeclareIO("ioName");
    io.SetEngine("Inline");
    adios2::Engine inlineWriter = io.Open("inline_write", adios2::Mode::Write);
    adios2::Engine inlineReader = io.Open("inline_read", adios2::Mode::Read);

Notice that unlike other engines, the reader and writer share an IO instance.
Both the writer and reader must be opened before either tries to call ``BeginStep()``/``PerformPuts()``/``PerformGets()``.
There must be exactly one writer, and exactly one reader.

For successful operation, the writer will perform a step, then the reader will perform a step in the same process.
When the reader starts its step, the only data it has available is that written by the writer in its process.
The reader then can retrieve whatever data was written by the writer by using the double-pointer ``Get`` call:

.. code-block:: c++

    void Engine::Get<T>(Variable<T>, T**) const;

This version of ``Get`` is only used for the inline engine.
See the example below for details.

.. note::
    Since the inline engine does not copy any data, the writer should avoid changing the data before the reader has read it.

Typical access pattern:

.. code-block:: c++

    // ... Application data generation

    inlineWriter.BeginStep();
    inlineWriter.Put(var, in_data); // always use deferred mode
    inlineWriter.EndStep();
    // Unlike other engines, data should not be reused at this point (since ADIOS
    // does not copy the data), though ADIOS cannot enforce this.
    // Must wait until reader is finished using the data.

    inlineReader.BeginStep();
    double* out_data;
    inlineReader.Get(var, &data);
    // Now in_data == out_data.
    inlineReader.EndStep();

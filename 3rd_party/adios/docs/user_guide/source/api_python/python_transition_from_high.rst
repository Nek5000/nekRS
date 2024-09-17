**********************************
Transition from old API to new API
**********************************

A python script using the high-level API of 2.9 and earlier needs to be modified to work with 2.10 and later.

- adios2.open() is replaced with adios2.Stream(), and does not have 4th and 5th optional arguments for external xml and IO name.
- the ``for in file`` is replaced with ``for _ in file.steps()`` but it works for both writing (by specifying the number of output steps) and reading (for the number of available steps in a stream/file).

.. code-block:: python

    # OLD API
    import adios2

    # NEW API
    from adios2 import Adios, Stream

    # NEW API: this still works
    import adios2


    # OLD API
    fr = adios2.open(args.instream, "r", mpi.comm_app,"adios2.xml", "SimulationOutput")

    # NEW API
    adios = Adios("adios2.xml", mpi.comm_app)
    io = adios.declare_io("SimulationOutput")
    fr = Stream(io, args.instream, "r", mpi.comm_app)


    # OLD API
    for fr_step in fr:
        fr_step....

    # NEW API 1
    for _ in fr.steps():
        fr....

    # NEW API 2
    for fr_step in fr.steps():
        fr_step....

#################
Supported Engines
#################

This section provides a description of the :ref:`Available Engines` in ADIOS2 and their specific parameters to allow extra-control from the user. Parameters are passed in key-value pairs for:

1. Engine specific parameters

2. Engine supported transports and parameters

Parameters are passed at:

1. Compile time ``IO::SetParameters`` (``adios2_set_parameter`` in C, Fortran)

2. Compile time ``IO::AddTransport`` (``adios2_set_transport_parameter`` in C, Fortran)

3. :ref:`Runtime Configuration Files` in the :ref:`ADIOS` component.

.. include:: bp5.rst
.. include:: bp4.rst
.. include:: bp3.rst
.. include:: hdf5.rst
.. include:: sst.rst
.. include:: ssc.rst
.. include:: dataman.rst
.. include:: dataspaces.rst
.. include:: inline.rst
.. include:: null.rst
.. include:: plugin.rst

###################
Supported Operators
###################

The Operator abstraction allows ADIOS2 to act upon the user application data,
either from a ``adios2::Variable`` or a set of Variables in an ``adios2::IO``
object.  Current supported operations are:

1. Data compression/decompression, lossy and lossless.

This section provides a description of the supported operators in ADIOS2
and their specific parameters to allow extra-control from the user. Parameters
are passed in key-value pairs for:

1. Operator general supported parameters.

2. Operator specific supported parameters.

Parameters are passed at:

1. Compile time: using the second parameter of the method
   ``ADIOS2::DefineOperator``

2. :ref:`Runtime Configuration Files` in the :ref:`ADIOS` component.

.. include:: CompressorZFP.rst
.. include:: plugin.rst
.. include:: encryption.rst

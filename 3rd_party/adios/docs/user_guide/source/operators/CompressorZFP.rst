*************
CompressorZFP
*************

The ``CompressorZFP`` Operator is compressor that uses a lossy but optionally
error-bounded compression to achieve high compression ratios.

ZFP provides compressed-array classes that support high throughput read and
write random access to individual array elements. ZFP also supports serial and
parallel (OpenMP and CUDA) compression of whole arrays, e.g., for applications
that read and write large data sets to and from disk.

ADIOS2 provides a ``CompressorZFP`` operator that lets you compress an
decompress variables. Below there is an example of how to invoke
``CompressorZFP`` operator:

.. code-block:: c++

    adios2::IO io = adios.DeclareIO("Output");
    auto ZFPOp    = adios.DefineOperator("CompressorZFP", adios2::ops::LossyZFP);

    auto var_r32 = io.DefineVariable<float>("r32", shape, start, count);
    var_r32.AddOperation(ZFPOp, {{adios2::ops::zfp::key::rate, rate}});

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CompressorZFP Specific parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``CompressorZFP`` operator accepts the following operator specific
parameters:

+-------------------+---------------------------------------------+
| ``CompressorZFP`` available parameters                          |
+===================+=============================================+
| ``accuracy``      | Fixed absolute error tolerance              |
+-------------------+---------------------------------------------+
| ``rate``          | Fixed number of bits in a compression unit  |
+-------------------+---------------------------------------------+
| ``precision``     | Fixed number of uncompressed bits per value |
+-------------------+---------------------------------------------+
| ``backend``       | Backend device: ``cuda`` ``omp`` ``serial`` |
+-------------------+---------------------------------------------+

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CompressorZFP Execution Policy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``CompressorZFP`` can run in multiple backend devices: GPUs (CUDA), OpenMP, and
in the host CPU. By default ``CompressorZFP`` will choose its backend following
the above order upon the availability of the device adapter.

Exceptionally, if its corresponding ADIOS2 variable contains a CUDA memory
address, this is a CUDA buffer, the CUDA backend will be called if available.

In any case, the user can manually set the backend using the ZFPOperator
specific parameter ``backend``.

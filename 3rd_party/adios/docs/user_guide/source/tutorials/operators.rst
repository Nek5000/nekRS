Operators
=========

.. _sec:tutorials_basics_operators:

In the previous tutorial we learned how to write and read attributes.

For this example to work, you would need to have the SZ compression library installed, which ADIOS automatically detects.
The easiest way to install SZ is with Spack, and you can do that as follows:

.. code-block:: bash

   git clone https://github.com/spack/spack.git ~/spack
   cd ~/spack
   . share/spack/setup-env.sh
   spack install sz
   spack load sz

In this tutorial we will learn how to use operators. Operators are used for Data compression/decompression, lossy and
lossless. They act upon the user application data, either from a variable or a set of variables in a IO object.

Additionally, we will explore how to simply write variables across multiple steps.

So, let's dig in!

Start editing the skeleton file `ADIOS2/examples/hello/bpOperatorSZWriter/bpOperatorSZWriter_tutorialSkeleton.cpp <https://github.com/ornladios/ADIOS2/blob/master/examples/hello/bpOperatorSZWriter/bpOperatorSZWriter_tutorialSkeleton.cpp>`_.

1. In an MPI application first we need to always initialize MPI. We do that with the following lines:

.. code-block:: cpp


   int rank, size;
   int rank, size;
   int provided;

   // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
   MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

2. This application has command line arguments for the size of the data, and the compression accuracy,
   which we can read as follows:

.. code-block:: cpp

   const std::size_t Nx = static_cast<std::size_t>(std::stoull(argv[1]));
   const double accuracy = std::stod(argv[2]);

3. Now we need to create some application variables which will be used to define ADIOS2 variables.

.. code-block:: cpp

   std::vector<float> myFloats(Nx);
   std::vector<double> myDoubles(Nx);
   std::iota(myFloats.begin(), myFloats.end(), 0.);
   std::iota(myDoubles.begin(), myDoubles.end(), 0.);

4. Now we need to create an ADIOS2 instance and IO object.

.. code-block:: cpp

   adios2::ADIOS adios(MPI_COMM_WORLD);
   adios2::IO bpIO = adios.DeclareIO("BPFile_SZ");

5. Now we need to define the variables we want to write.

.. code-block:: cpp

   adios2::Variable<float> bpFloats = bpIO.DefineVariable<float>(
       "bpFloats", {size * Nx}, {rank * Nx}, {Nx}, adios2::ConstantDims);
   adios2::Variable<double> bpDoubles = bpIO.DefineVariable<double>(
       "bpDoubles", {size * Nx}, {rank * Nx}, {Nx}, adios2::ConstantDims);

6. Now we need to define the compression operator we want to use. In this case we will use the SZ compressor.

.. code-block:: cpp

   adios2::Operator op = bpIO.DefineOperator("SZCompressor", "sz");
   varFloats.AddOperation(op, {{"accuracy", std::to_string(accuracy)}});
   varDoubles.AddOperation(op, {{"accuracy", std::to_string(accuracy)}});

.. note::

   ``DefineOperator()'`` s second parameter can be either zfp or sz. For more information regarding operators and their
   properties you can look at :ref:`Basics: Interface Components: Operator <sec:basics_interface_components_operator>`.

7. Let's also create an attribute to store the accuracy value.

.. code-block:: cpp

   adios2::Attribute<double> attribute = bpIO.DefineAttribute<double>("accuracy", accuracy);

8. Now we need to open the file for writing.

.. code-block:: cpp

   adios2::Engine bpWriter = bpIO.Open("SZexample.bp", adios2::Mode::Write);

9. Now we need to write the data. We will write the data for 3 steps, and edit them in between.

.. code-block:: cpp

   for (unsigned int step = 0; step < 3; ++step)
   {
       bpWriter.BeginStep();

       bpWriter.Put<double>(bpDoubles, myDoubles.data());
       bpWriter.Put<float>(bpFloats, myFloats.data());

       bpWriter.EndStep();

       // here you can modify myFloats, myDoubles per step
       std::transform(myFloats.begin(), myFloats.end(), myFloats.begin(),
                      [&](float v) -> float { return 2 * v; });
       std::transform(myDoubles.begin(), myDoubles.end(), myDoubles.begin(),
                      [&](double v) -> double { return 3 * v; });
   }

10. Now we need to close the file.

.. code-block:: cpp

   bpWriter.Close();

11. Finally we need to finalize MPI.

.. code-block:: cpp

   MPI_Finalize();

12. The final code should look as follows (excluding try/catch and optional usage of MPI), and it was derived from the
    example `ADIOS2/examples/hello/bpOperatorSZWriter/bpOperatorSZWriter.cpp <https://github.com/ornladios/ADIOS2/blob/master/examples/hello/bpOperatorSZWriter/bpOperatorSZWriter.cpp>`_.

.. literalinclude:: ../../../../examples/hello/bpOperatorSZWriter/bpOperatorSZWriter.cpp
   :language: cpp

13. You can compile and run it as follows:

.. code-block:: bash

   cd Path-To-ADIOS2/examples/hello/bpOperatorSZWriter
   mkdir build
   cd build
   cmake -DADIOS2_DIR=Path-To-ADIOS2/build/ ..
   cmake --build .
   mpirun -np 2 ./adios2_hello_bpOperatorSZWriter_mpi 20 0.000001

12. You can check the content of the output file "SZexample.bp" using *bpls* as follows:

.. code-block:: bash

   Path-To-ADIOS2/build/bin/bpls ./SZexample.bp

     double   bpDoubles  3*{40}
     float    bpFloats   3*{40}

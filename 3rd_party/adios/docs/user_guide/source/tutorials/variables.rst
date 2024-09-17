Variables
=========

.. _sec:tutorials_basics_variables:

In the previous tutorial we learned how to define a simple string variable, write it, and read it back.

In this tutorial we will go two steps further:

1. We will define variables which include arrays, and we will write them and read them back.
2. We will use MPI to write and read the above variables in parallel.

Let's start with the writing part.

Start editing the skeleton file `ADIOS2/examples/hello/bpWriter/bpWriter_tutorialSkeleton.cpp <https://github.com/ornladios/ADIOS2/blob/master/examples/hello/bpWriter/bpWriter_tutorialSkeleton.cpp>`_.

1. In an MPI application first we need to always initialize MPI. We do that with the following lines:

.. code-block:: cpp

   int rank, size;
   int provided;

   // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
   MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

2. Now we need to create some application variables which will be used to define ADIOS2 variables.

.. code-block:: cpp

   // Application variable
   std::vector<float> myFloats = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
   std::vector<int> myInts = {0, -1, -2, -3, -4, -5, -6, -7, -8, -9};
   const std::size_t Nx = myFloats.size();
   const std::string myString("Hello Variable String from rank " + std::to_string(rank));

3. Now we need to define an ADIOS2 instance and the ADIOS2 variables.

.. code-block:: cpp

   adios2::ADIOS adios(MPI_COMM_WORLD);
   adios2::IO bpIO = adios.DeclareIO("BPFile_N2N");

   adios2::Variable<float> bpFloats = bpIO.DefineVariable<float>(
       "bpFloats", {size * Nx}, {rank * Nx}, {Nx}, adios2::ConstantDims);

   adios2::Variable<int> bpInts = bpIO.DefineVariable<int>("bpInts", {size * Nx}, {rank * Nx},
                                                           {Nx}, adios2::ConstantDims);

   // For the sake of the tutorial we create an unused variable
   adios2::Variable<std::string> bpString = bpIO.DefineVariable<std::string>("bpString");

.. note::

   The above int/float variables are global arrays. The 1st argument of the ``DefineVariable`` function is the variable
   name, the 2nd are the global dimensions, the 3rd is the start index for a rank, the 4th are the rank/local
   dimensions, and the 5th is a boolean variable to indicate if the dimensions are constant or not over multiple steps,
   where ``adios2::ConstantDims == true`` We will explore other tutorials that don't use constant dimensions.

4. Now we need to open the ADIOS2 engine and write the variables.

.. code-block:: cpp

   adios2::Engine bpWriter = bpIO.Open("myVector_cpp.bp", adios2::Mode::Write);

   bpWriter.BeginStep();
   bpWriter.Put(bpFloats, myFloats.data());
   bpWriter.Put(bpInts, myInts.data());
   // bpWriter.Put(bpString, myString);
   bpWriter.EndStep();

   bpWriter.Close();

5. Finally we need to finalize MPI.

.. code-block:: cpp

   MPI_Finalize();

6. The final code should look as follows (excluding try/catch and the optional usage of MPI), and it was derived
   from the example `ADIOS2/examples/hello/bpWriter/bpWriter.cpp <https://github.com/ornladios/ADIOS2/blob/master/examples/hello/bpWriter/bpWriter.cpp>`_.

.. literalinclude:: ../../../../examples/hello/bpWriter/bpWriter.cpp
   :language: cpp

7. You can compile and run it as follows:

.. code-block:: bash

   cd Path-To-ADIOS2/examples/hello/bpWriter
   mkdir build
   cd build
   cmake -DADIOS2_DIR=Path-To-ADIOS2/build/ ..
   cmake --build .
   mpirun -np 2 ./adios2_hello_bpWriter_mpi

8. You can check the content of the output file "myVector_cpp.bp" using *bpls* as follows:

.. code-block:: bash

   Path-To-ADIOS2/build/bin/bpls ./myVector_cpp.bp

     float    bpFloats  {10}
     int32_t  bpInts    {10}

Now let's move to the reading part.

Start editing the skeleton file `ADIOS2/examples/hello/bpReader/bpReader_tutorialSkeleton.cpp <https://github.com/ornladios/ADIOS2/blob/master/examples/hello/bpReader/bpReader_tutorialSkeleton.cpp>`_.

9. In an MPI application first we need to always initialize MPI. We do that with the following line:

.. code-block:: cpp

   int rank, size;
   int provided;

   // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
   MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

10. Now we need to define an ADIOS2 instance and open the ADIOS2 engine.

.. code-block:: cpp

   adios2::ADIOS adios(MPI_COMM_WORLD);
   
   adios2::IO bpIO = adios.DeclareIO("BPFile_N2N");
   
   adios2::Engine bpReader = bpIO.Open("myVector_cpp.bp", adios2::Mode::Read);

11. Now we need to read the variables. In this case we know the variables that we need to inquire, so we can use the
    ``InquireVariable`` function immediately. But let's explore how to check the available variables in a file first,
    and then we will use the ``InquireVariable`` function.

.. code-block:: cpp

   bpReader.BeginStep();
   const std::map<std::string, adios2::Params> variables = bpIO.AvailableVariables();

   for (const auto &variablePair : variables)
   {
       std::cout << "Name: " << variablePair.first;
       for (const auto &parameter : variablePair.second)
       {
           std::cout << "\t" << parameter.first << ": " << parameter.second << "\n";
       }
   }

   adios2::Variable<float> bpFloats = bpIO.InquireVariable<float>("bpFloats");
   adios2::Variable<int> bpInts = bpIO.InquireVariable<int>("bpInts");

12. Now we need to read the variables from each rank. We will use the ``SetSelection`` to set the start index and rank
    dimensions, then ``Get`` function to read the variables, and print the contents from rank 0.

.. code-block:: cpp

   const std::size_t Nx = 10;
   if (bpFloats) // means found
   {
      std::vector<float> myFloats;

      // read only the chunk corresponding to our rank
      bpFloats.SetSelection({{Nx * rank}, {Nx}});
      bpReader.Get(bpFloats, myFloats, adios2::Mode::Sync);

      if (rank == 0)
      {
          std::cout << "MyFloats: \n";
          for (const auto number : myFloats)
          {
              std::cout << number << " ";
          }
          std::cout << "\n";
      }
   }

   if (bpInts) // means not found
   {
      std::vector<int> myInts;
      // read only the chunk corresponding to our rank
      bpInts.SetSelection({{Nx * rank}, {Nx}});

      bpReader.Get(bpInts, myInts, adios2::Mode::Sync);

      if (rank == 0)
      {
          std::cout << "myInts: \n";
          for (const auto number : myInts)
          {
              std::cout << number << " ";
          }
          std::cout << "\n";
      }
   }

.. note::

   While using the ``Get`` function, we used the third parameter named ``Mode``. The mode parameter can also be used
   for the ``Put`` function.

   For the ``Put`` function, there are three modes: ``Deferred`` (default), ``Sync``, and ``Span``. and for the ``Get``
   there are two modes: ``Deferred`` (default) and ``Sync``.

   1. The ``Deferred`` mode is the default mode, because it is the fastest mode, as it allows ``Put`` / ``Get`` to be
      grouped before potential data transport at the first encounter of ``PerformPuts`` / ``PerformGets``, ``EndStep``
      or ``Close``.

   2. The ``Sync`` mode forces ``Put`` / ``Get`` to be performed immediately so that the data are available immediately.

   3. The ``Span`` mode is special mode of ``Deferred`` that allows population from non-contiguous memory structures.

   For more information about the ``Mode`` parameter for both ``Put`` and ``Get`` functions, and when you should use
   each option see :ref:`Basics: Interface Components: Engine <sec:basics_interface_components_engine>`.

13. Now we need close the ADIOS2 engine.

.. code-block:: cpp

   bpReader.EndStep();
   bpReader.Close();

14. Finally we need to finalize MPI.

.. code-block:: cpp

   MPI_Finalize();

15. The final code should look as follows (excluding try/catch), and it was derived from the example
    `ADIOS2/examples/hello/bpWriter/bpWriter.cpp <https://github.com/ornladios/ADIOS2/blob/master/examples/hello/bpWriter/bpWriter.cpp>`_.

.. literalinclude:: ../../../../examples/hello/bpReader/bpReader.cpp
   :language: cpp

16. You can compile and run it as follows:

.. code-block:: bash

   cd Path-To-ADIOS2/examples/hello/bpReader
   mkdir build
   cd build
   cmake -DADIOS2_DIR=Path-To-ADIOS2/build/ ..
   cmake --build .
   mpirun -np 2 ./adios2_hello_bpReader_mpi

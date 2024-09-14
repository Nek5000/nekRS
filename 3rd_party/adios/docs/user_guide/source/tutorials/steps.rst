Steps
=====

.. _sec:tutorials_basics_steps:

In the previous tutorial, we introduced the concept of operators, and briefly touched upon the concept of steps.

In this tutorial, we will explore how to write data for multiple steps, and how to read them back.

So let's dig in!

Start editing the skeleton file `ADIOS2/examples/hello/bpStepsWriteRead/bpStepsWriteRead_tutorialSkeleton.cpp <https://github.com/ornladios/ADIOS2/blob/master/examples/hello/bpStepsWriteRead/bpStepsWriteRead_tutorialSkeleton.cpp>`_.

1. In an MPI application first we need to always initialize MPI. We do that with the following lines:

.. code-block:: cpp

   int rank, size;
   int rank, size;
   int provided;

   // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
   MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

2. This application has an optional command line argument for engine being used. If
   no argument is provided, the default engine is BPFile.

.. code-block:: cpp

    const std::string engine = argv[1] ? argv[1] : "BPFile";

3. We will define the number of steps and the size of the data that we will create.

.. code-block:: cpp

    const std::string filename = engine + "StepsWriteRead.bp";
    const unsigned int nSteps = 10;
    const unsigned int Nx = 60000;

4. Now we need to create an ADIOS2 instance.

.. code-block:: cpp

   adios2::ADIOS adios(MPI_COMM_WORLD);

5. Now we will populate the writer function with the following signature:

.. code-block::

   void writer(adios2::ADIOS &adios, const std::string &engine, const std::string &fname,
             const size_t Nx, unsigned int nSteps, int rank, int size)
   {
     ...
   }

6. Let's create some simulation data. We will create a 1D array of size Nx, and fill it with 0.block

.. code-block:: cpp

   std::vector<double> simData(Nx, 0.0);

7. Now we will create an IO object and set the engine type.block

.. code-block:: cpp

   adios2::IO bpIO = adios.DeclareIO("SimulationOutput");
   io.SetEngine(engine);

.. note::

   The beauty of ADIOS2 is that you write the same code for all engines. The only thing that changes is the engine name.
   The underlying engine handles all the intricacies of the engine's format, and the user enjoys the API's simplicity.

8. Now we will create a variable for the simulation data and the step.

.. code-block:: cpp

    const adios2::Dims shape{static_cast<size_t>(size * Nx)};
    const adios2::Dims start{static_cast<size_t>(rank * Nx)};
    const adios2::Dims count{Nx};
    auto bpFloats = bpIO.DefineVariable<float>("bpFloats", shape, start, count);

    auto bpStep = bpIO.DefineVariable<unsigned int>("bpStep");

9. Now we will open the file for writing.

.. code-block:: cpp

    adios2::Engine bpWriter = bpIO.Open(fname, adios2::Mode::Write);

10. Now we will write the data for each step.

.. code-block:: cpp

   for (unsigned int step = 0; step < nSteps; ++step)
   {
       const adios2::Box<adios2::Dims> sel({0}, {Nx});
       bpFloats.SetSelection(sel);

       bpWriter.BeginStep();
       bpWriter.Put(bpFloats, simData.data());
       bpWriter.Put(bpStep, step);
       bpWriter.EndStep();

       // Update values in the simulation data
       update_array(simData, 10);
   }

11. Now we will close the file.

.. code-block:: cpp

   bpWriter.Close();

12. Now we will populate the reader function with the following signature:

.. code-block:: cpp

   void reader(adios2::ADIOS &adios, const std::string &engine, const std::string &fname,
               const size_t Nx, unsigned int /*nSteps*/, int rank, int /*size*/)
   {
     ...
   }

13. Now we will create an IO object and set the engine type.

.. code-block:: cpp

   adios2::IO bpIO = adios.DeclareIO("SimulationOutput");
   io.SetEngine(engine);

14. Now we will open the file for reading.

.. code-block:: cpp

   adios2::Engine bpReader = bpIO.Open(fname, adios2::Mode::Read);

15. Now we will create a vector to store simData and a variable for the step.

.. code-block:: cpp

   std::vector<float> simData(Nx, 0);
   unsigned int inStep = 0;

16. Now we will read the data for each step.

.. code-block:: cpp

   for (unsigned int step = 0; bpReader.BeginStep() == adios2::StepStatus::OK; ++step)
   {
       auto bpFloats = bpIO.InquireVariable<float>("bpFloats");
       if (bpFloats)
       {
           const adios2::Box<adios2::Dims> sel({{Nx * rank}, {Nx}});
           bpFloats.SetSelection(sel);
           bpReader.Get(bpFloats, simData.data());
       }
       auto bpStep = bpIO.InquireVariable<unsigned int>("bpStep");
       if (bpStep)
       {
           bpReader.Get(bpStep, &inStep);
       }

       bpReader.EndStep();
   }

17. Now we will close the file.

.. code-block:: cpp

   bpReader.Close();

18. Now we will call the writer and reader functions:

.. code-block:: cpp

   writer(adios, engine, filename, Nx, nSteps, rank, size);
   reader(adios, engine, filename, Nx, nSteps, rank, size);

19. Finally we need to finalize MPI.

.. code-block:: cpp

   MPI_Finalize();

20. The final code should look as follows (excluding try/catch and optional usage of MPI), and it was derived from the
    example `ADIOS2/examples/hello/bpStepsWriteRead/bpStepsWriteRead.cpp <https://github.com/ornladios/ADIOS2/blob/master/examples/hello/bpStepsWriteRead/bpStepsWriteRead.cpp>`_.

.. literalinclude:: ../../../../examples/hello/bpStepsWriteRead/bpStepsWriteRead.cpp
   :language: cpp

21. You can compile and run it as follows:

.. code-block:: bash

   cd Path-To-ADIOS2/examples/hello/bpStepsWriteRead
   mkdir build
   cd build
   cmake -DADIOS2_DIR=Path-To-ADIOS2/build/ ..
   cmake --build .
   mpirun -np 2 ./adios2_hello_bpStepsWriteRead_mpi

22. You can check the content of the output file "BPFileStepsWriteRead.bp" using *bpls* as follows:

.. code-block:: bash

   Path-To-ADIOS2/build/bin/bpls ./BPFileStepsWriteRead.bp

     float     bpFloats  10*{120000}
     uint32_t  bpStep    10*scalar

Attributes
==========

.. _sec:tutorials_basics_attributes:

In the previous tutorial, we learned how to write/read variables.

In this tutorial, we will explore how to write/read attributes. Attributes are metadata related to the whole dataset or
to a specific variable. In this tutorial, we will only focus on attributes related to the whole dataset, but we will
explain how variable's attributes can be used too.

Start editing the skeleton file `ADIOS2/examples/hello/bpAttributeWriteRead/bpAttributeWriteRead_tutorialSkeleton.cpp <https://github.com/ornladios/ADIOS2/blob/master/examples/hello/bpAttributeWriteRead/bpAttributeWriteRead_tutorialSkeleton.cpp>`_.

1. In an MPI application first we need to always initialize MPI. We do that with the following lines:

.. code-block:: cpp

   int rank, size;
   int provided;

   // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
   MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

2.  Now we need to create a application variable which will be used to define an ADIOS2 variable.

.. code-block:: cpp

   std::vector<float> myFloats = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

3. Then, we need to create an ADIOS2 instance.

.. code-block:: cpp

   adios2::ADIOS adios(MPI_COMM_WORLD);

4. Then, we create the following writer function:

.. code-block:: cpp

   void writer(adios2::ADIOS &adios, int rank, int size, std::vector<float> &myFloats)
   {
      ...
   }

5. In this writer function, we define an IO object, and a float vector variable as follows:

.. code-block:: cpp

   adios2::IO bpIO = adios.DeclareIO("BPFile_N2N");

   const std::size_t Nx = myFloats.size();
   adios2::Variable<float> bpFloats = bpIO.DefineVariable<float>(
     "bpFloats", {size * Nx}, {rank * Nx}, {Nx}, adios2::ConstantDims);

6. Now, we will define various types of attributes as follows:

.. code-block:: cpp

   bpIO.DefineAttribute<std::string>("Single_String", "File generated with ADIOS2");

   std::vector<std::string> myStrings = {"one", "two", "three"};
   bpIO.DefineAttribute<std::string>("Array_of_Strings", myStrings.data(), myStrings.size());

   bpIO.DefineAttribute<double>("Attr_Double", 0.f);
   std::vector<double> myDoubles = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
   bpIO.DefineAttribute<double>("Array_of_Doubles", myDoubles.data(), myDoubles.size());

.. note::

   if we want to define an attribute for a specific variable, we can use one of the following API:

   .. code-block:: cpp

      template <class T>
      Attribute<T> DefineAttribute(const std::string &name, const T *data, const size_t size,
                                   const std::string &variableName = "", const std::string separator = "/",
                                   const bool allowModification = false);

      template <class T>
      Attribute<T> DefineAttribute(const std::string &name, const T &value,
                                   const std::string &variableName = "", const std::string separator = "/",
                                   const bool allowModification = false);

   As we can see, by default the attributes don't change over multiple steps, but we can change that by setting
   ``allowModification`` to ``true``.

7. Then, we open a file for writing:

.. code-block:: cpp

   adios2::Engine bpWriter = bpIO.Open("fileAttributes.bp", adios2::Mode::Write);

8. Now, we write the data and close the file:

.. code-block:: cpp

   bpWriter.BeginStep();
   bpWriter.Put<float>(bpFloats, myFloats.data());
   bpWriter.EndStep();
   bpWriter.Close();

9. Steps 1-8 are used for writing, we will define a reader function in the rest of the steps:

.. code-block:: cpp

   void reader(adios2::ADIOS &adios, int rank, int size)
   {
      ...
   }

10. In this reader function, we define an IO object, and open the file for reading:

.. code-block:: cpp

   adios2::IO bpIO = adios.DeclareIO("BPFile_N2N");
   adios2::Engine bpReader = bpIO.Open("fileAttributes.bp", adios2::Mode::Read);

11. Now, we check the available attributes as follows:

.. code-block:: cpp

   bpReader.BeginStep();
   const auto attributesInfo = bpIO.AvailableAttributes();

   for (const auto &attributeInfoPair : attributesInfo)
   {
     std::cout << "Attribute: " << attributeInfoPair.first;
     for (const auto &attributePair : attributeInfoPair.second)
     {
         std::cout << "\tKey: " << attributePair.first << "\tValue: " << attributePair.second
                   << "\n";
     }
     std::cout << "\n";
   }

12. Now we will inquire and get the attributes as follows:

.. code-block:: cpp

    adios2::Attribute<float> singleString = bpIO.InquireAttribute<float>("Single_String");
    if (singleString)
    {
        std::cout << singleString.Name() << ": " << singleString.Data()[0] << "\n";
    }
    adios2::Attribute<std::string> arrayOfStrings =
        bpIO.InquireAttribute<std::string>("Array_of_Strings");
    if (arrayOfStrings)
    {
        std::cout << arrayOfStrings.Name() << ": ";
        for (const auto &value : arrayOfStrings.Data())
        {
            std::cout << value << " ";
        }
        std::cout << "\n";
    }
    adios2::Attribute<double> attrDouble = bpIO.InquireAttribute<double>("Attr_Double");
    if (attrDouble)
    {
        std::cout << attrDouble.Name() << ": " << attrDouble.Data()[0] << "\n";
    }
    adios2::Attribute<double> arrayOfDoubles = bpIO.InquireAttribute<double>("Array_of_Doubles");
    if (arrayOfDoubles)
    {
        std::cout << arrayOfDoubles.Name() << ": ";
        for (const auto &value : arrayOfDoubles.Data())
        {
            std::cout << value << " ";
        }
        std::cout << "\n";
    }

13. Afterward, we will inquire and get the variable as follows:

.. code-block:: cpp

    adios2::Variable<float> bpFloats = bpIO.InquireVariable<float>("bpFloats");
    const std::size_t Nx = 10;
    std::vector<float> myFloats(Nx);
    if (bpFloats)
    {
        bpFloats.SetSelection({{Nx * rank}, {Nx}});
        bpReader.Get(bpFloats, myFloats.data());
    }
    bpReader.EndStep();

14. Finally, we close the file:

.. code-block:: cpp

   bpReader.Close();

15. In the main function, we call the writer and reader functions as follows:

.. code-block:: cpp

   writer(adios, rank, size, myFloats);
   reader(adios, rank, size);

16. Finally, we finalize MPI:

.. code-block:: cpp

   MPI_Finalize();

17. The final code should look as follows (excluding try/catch and optional usage MPI), and it was derived from the
    example `ADIOS2/examples/hello/bpAttributeWriteRead/bpAttributeWriteRead.cpp <https://github.com/ornladios/ADIOS2/blob/master/examples/hello/bpAttributeWriteRead/bpAttributeWriteRead.cpp>`_.

.. literalinclude:: ../../../../examples/hello/bpAttributeWriteRead/bpAttributeWriteRead.cpp
   :language: cpp

18. You can compile and run it as follows:

.. code-block:: bash

   cd Path-To-ADIOS2/examples/hello/bpAttributeWriteRead
   mkdir build
   cd build
   cmake -DADIOS2_DIR=Path-To-ADIOS2/build/ ..
   cmake --build .
   mpirun -np 2 ./adios2_hello_bpAttributeWriteRead_mpi

19. You can check the content of the output file "fileAttributes.bp" using *bpls* as follows:

.. code-block:: bash

   Path-To-ADIOS2/build/bin/bpls ./fileAttributes.bp

     float    bpFloats  {20}

******************
C++ High-Level API
******************

C++11 High-Level APIs are based on a single object `adios2::fstream`

.. caution::

   DO NOT place ``use namespace adios2`` in your C++ code.
   Use ``adios2::fstream`` directly to prevent conflicts with ``std::fstream``. 


C++11 Write example
-------------------

.. code-block:: c++

   #include <adios2.h>
   ...

   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   
   // Nx, Ny from application, std::size_t
   const adios2::Dims shape{Nx, Ny * static_cast<std::size_t>(size)};
   const adios2::Dims start{0, Ny * static_cast<std::size_t>(rank)};
   const adios2::Dims count{Nx, Ny};
   
   adios2::fstream oStream("cfd.bp", adios2::fstream::out, MPI_COMM_WORLD);

   // NSteps from application
   for (std::size_t step = 0; step < NSteps; ++step)
   {
       if(rank == 0 && step == 0) // global variable
       {
           oStream.write<int32_t>("size", size);
       }

       // physicalTime double, <double> is optional
       oStream.write<double>( "physicalTime", physicalTime );
       // T and P are std::vector<float>
       oStream.write( "temperature", T.data(), shape, start, count );
       // adios2::endl will advance the step after writing pressure
       oStream.write( "pressure", P.data(), shape, start, count, adios2::end_step );
       
   }
   
   // Calling close is mandatory!
   oStream.close(); 

C++11 Read "step-by-step" example
---------------------------------

.. code-block:: c++

   #include <adios2.h>
   ...
   
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   
   // Selection Window from application, std::size_t
   const adios2::Dims start{0, 0};
   const adios2::Dims count{SelX, SelY};
   
   if( rank == 0)
   {
      // if only one rank is active use MPI_COMM_SELF
      adios2::fstream iStream("cfd.bp", adios2::fstream::in, MPI_COMM_SELF);
   
      adios2::fstep iStep;
      while (adios2::getstep(iStream, iStep))
      {
          if( iStep.currentstep() == 0 )
          {
              const std::size_t sizeOriginal = iStep.read<std::size_t>("size");
          }
          const double physicalTime = iStream.read<double>( "physicalTime");
          const std::vector<float> temperature = iStream.read<float>( "temperature", start, count );
          const std::vector<float> pressure = iStream.read<float>( "pressure", start, count );
      }
      // Don't forget to call close!
      iStream.close(); 
   }
   
``adios2::fstream`` API documentation
-------------------------------------

.. doxygenclass:: adios2::fstream
   :project: CXX11
   :path: ../../bindings/CXX11/cxx11/fstream
   :members:

/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * inlineReaderWriter.cpp  example borrowed from bpTimeWriter, using
 * the inline engine. Writes a variable using the Advance function for time
 * aggregation. Time step is saved as an additional (global) single value
 * variable, just for tracking purposes.
 *
 *  Created on: Nov 16, 2018
 *      Author: Aron Helser aron.helser@kitware.com
 */

#include <algorithm> //std::for_each
#include <ios>       //std::ios_base::failure
#include <iostream>  //std::cout
#include <stdexcept> //std::invalid_argument std::exception
#include <vector>

#include <adios2.h>

void DoAnalysis(adios2::IO &inlineIO, adios2::Engine &inlineReader, int rank, unsigned int step)
{
    inlineReader.BeginStep();
    /////////////////////READ
    adios2::Variable<float> inlineFloats000 = inlineIO.InquireVariable<float>("inlineFloats000");

    adios2::Variable<std::string> inlineString =
        inlineIO.InquireVariable<std::string>("inlineString");

    if (inlineFloats000)
    {
        auto blocksInfo = inlineReader.BlocksInfo(inlineFloats000, step);

        std::cout << "Data StepsStart " << inlineFloats000.StepsStart() << " from rank " << rank
                  << ": ";
        for (auto &info : blocksInfo)
        {
            // bp file reader would see all blocks, inline only sees local
            // writer's block(s).
            size_t myBlock = info.BlockID;
            inlineFloats000.SetBlockSelection(myBlock);

            // info passed by reference
            // engine must remember data pointer (or info) to fill it out at
            // PerformGets()
            inlineReader.Get<float>(inlineFloats000, info, adios2::Mode::Deferred);
        }
        inlineReader.PerformGets();

        for (const auto &info : blocksInfo)
        {
            adios2::Dims count = info.Count;
            const float *vectData = info.Data();
            for (size_t i = 0; i < count[0]; ++i)
            {
                float datum = vectData[i];
                std::cout << datum << " ";
            }
            std::cout << "\n";
        }
    }
    else
    {
        std::cout << "Variable inlineFloats000 not found\n";
    }

    if (inlineString && rank == 0)
    {
        std::string myString;
        inlineReader.Get(inlineString, myString, adios2::Mode::Sync);
        std::cout << "inlineString: " << myString << "\n";
    }
    inlineReader.EndStep();
    // all deferred block info are now valid - need data pointers to be
    // valid, filled with data
}

int main(int argc, char *argv[])
{
    int rank = 0, size = 1;
#if ADIOS2_USE_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif

    // Application variable
    std::vector<float> myFloats = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    const std::size_t Nx = myFloats.size();

    try
    {
        // Inline uses single IO for write/read
        adios2::IO inlineIO = adios.DeclareIO("InlineReadWrite");
        /// WRITE
        inlineIO.SetEngine("Inline");
        inlineIO.SetParameter("verbose", "4");

        /** global array: name, { shape (total dimensions) }, { start
         * (local) },
         * { count (local) }, all are constant dimensions */
        const unsigned int variablesSize = 10;
        std::vector<adios2::Variable<float>> inlineFloats(variablesSize);

        adios2::Variable<std::string> inlineString =
            inlineIO.DefineVariable<std::string>("inlineString");

        for (unsigned int v = 0; v < variablesSize; ++v)
        {
            std::string namev("inlineFloats");
            if (v < 10)
            {
                namev += "00";
            }
            else if (v < 100)
            {
                namev += "0";
            }
            namev += std::to_string(v);

            inlineFloats[v] = inlineIO.DefineVariable<float>(namev, {size * Nx}, {rank * Nx}, {Nx},
                                                             adios2::ConstantDims);
        }

        /** global single value variable: name */
        adios2::Variable<unsigned int> inlineTimeStep =
            inlineIO.DefineVariable<unsigned int>("timeStep");

        adios2::Engine inlineWriter = inlineIO.Open("myWriteID", adios2::Mode::Write);

        adios2::Engine inlineReader = inlineIO.Open("myReadID", adios2::Mode::Read);

        for (unsigned int timeStep = 0; timeStep < 3; ++timeStep)
        {
            inlineWriter.BeginStep();
            if (rank == 0) // global single value, only saved by rank 0
            {
                inlineWriter.Put<unsigned int>(inlineTimeStep, timeStep);
            }

            // template type is optional, but recommended
            for (unsigned int v = 0; v < variablesSize; ++v)
            {
                // Note: Put is deferred, so all variables will see v == 9
                // and myFloats[0] == 9, 10, or 11
                myFloats[rank] = static_cast<float>(v + timeStep + rank);
                inlineWriter.Put(inlineFloats[v], myFloats.data());
            }

            const std::string myString("Hello from rank: " + std::to_string(rank) +
                                       " and timestep: " + std::to_string(timeStep));

            if (rank == 0)
            {
                inlineWriter.Put(inlineString, myString);
            }

            inlineWriter.EndStep();

            DoAnalysis(inlineIO, inlineReader, rank, timeStep);
        }
    }
    catch (std::exception const &e)
    {
        std::cout << "Caught exception from rank " << rank << "\n";
        std::cout << e.what() << "\n";
#if ADIOS2_USE_MPI
        return MPI_Abort(MPI_COMM_WORLD, 1);
#else
        return 1;
#endif
    }

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}

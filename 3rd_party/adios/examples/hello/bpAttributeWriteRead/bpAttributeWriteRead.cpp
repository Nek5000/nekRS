/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * bpAttributeWriteRead.cpp: Simple self-descriptive example of how to write/read attributes and
 * a variable to a BP File that lives in several MPI processes.
 *
 *  Created on: Feb 16, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include <ios>      //std::ios_base::failure
#include <iostream> //std::cout
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif
#include <stdexcept> //std::invalid_argument std::exception
#include <string>
#include <vector>

#include <adios2.h>

void writer(adios2::ADIOS &adios, int rank, int size, std::vector<float> &myFloats)
{
    /*** IO class object: settings and factory of Settings: Variables,
     * Parameters, Transports, and Execution: Engines */
    adios2::IO bpIO = adios.DeclareIO("BPFile_N2N");

    const std::size_t Nx = myFloats.size();

    /** global array : name, { shape (total) }, { start (local) }, { count
     * (local) }, all are constant dimensions */
    adios2::Variable<float> bpFloats = bpIO.DefineVariable<float>(
        "bpFloats", {size * Nx}, {rank * Nx}, {Nx}, adios2::ConstantDims);

    bpIO.DefineAttribute<std::string>("Single_String", "File generated with ADIOS2");

    std::vector<std::string> myStrings = {"one", "two", "three"};
    bpIO.DefineAttribute<std::string>("Array_of_Strings", myStrings.data(), myStrings.size());

    bpIO.DefineAttribute<double>("Attr_Double", 0.f);
    std::vector<double> myDoubles = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    bpIO.DefineAttribute<double>("Array_of_Doubles", myDoubles.data(), myDoubles.size());

    /** Engine derived class, spawned to start IO operations */
    adios2::Engine bpWriter = bpIO.Open("fileAttributes.bp", adios2::Mode::Write);

    bpWriter.BeginStep();

    /** Write variable for buffering */
    bpWriter.Put<float>(bpFloats, myFloats.data());

    bpWriter.EndStep();

    /** Create bp file, engine becomes unreachable after this*/
    bpWriter.Close();
}

void reader(adios2::ADIOS &adios, int rank, int /*size*/)
{
    adios2::IO bpIO = adios.DeclareIO("BPReader");

    adios2::Engine bpReader = bpIO.Open("fileAttributes.bp", adios2::Mode::Read);

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

    adios2::Variable<float> bpFloats = bpIO.InquireVariable<float>("bpFloats");
    const std::size_t Nx = 10;
    std::vector<float> myFloats(Nx);
    if (bpFloats)
    {
        bpFloats.SetSelection({{Nx * rank}, {Nx}});
        bpReader.Get(bpFloats, myFloats.data());
    }

    bpReader.EndStep();

    bpReader.Close();
}

int main(int argc, char *argv[])
{
    int rank, size;
#if ADIOS2_USE_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
    rank = 0;
    size = 1;
#endif

    /** Application variable */
    std::vector<float> myFloats = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    try
    {
        /** ADIOS class factory of IO class objects */
#if ADIOS2_USE_MPI
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif

        writer(adios, rank, size, myFloats);
        reader(adios, rank, size);
    }
    catch (std::invalid_argument &e)
    {
        std::cout << "Invalid argument exception, STOPPING PROGRAM from rank " << rank << "\n";
        std::cout << e.what() << "\n";
#if ADIOS2_USE_MPI
        std::cerr << "STOPPING PROGRAM from rank " << rank << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }
    catch (std::ios_base::failure &e)
    {
        std::cout << "IO System base failure exception, STOPPING PROGRAM from rank " << rank
                  << "\n";
        std::cout << e.what() << "\n";
#if ADIOS2_USE_MPI
        std::cerr << "STOPPING PROGRAM from rank " << rank << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }
    catch (std::exception &e)
    {
        std::cout << "Exception, STOPPING PROGRAM from rank " << rank << "\n";
        std::cout << e.what() << "\n";
#if ADIOS2_USE_MPI
        std::cerr << "STOPPING PROGRAM from rank " << rank << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}

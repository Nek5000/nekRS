/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * hdf5SubFile.cpp: Simple self-descriptive example of how to write a
 * variable to a parallel HDF5 File using MPI processes.
 *
 *  Created on: March 6, 2023
 *      Author: Junmin
 */

#include <adios2.h>
#include <ios>      //std::ios_base::failure
#include <iostream> //std::cout
#include <mpi.h>
#include <stdexcept> //std::invalid_argument std::exception
#include <vector>

void writeMe(adios2::IO &hdf5IO, int rank, int size, const char *testFileName)
{
    /** Application variable */
    int scale = 1;

    const char *temp = std::getenv("TEST_SCALE");
    if (NULL != temp)
    {
        int itemp = -1;
        sscanf(temp, "%d", &itemp);
        if (itemp > 1)
            scale = itemp;
    }

    const std::size_t Nx = 1024;
    const std::size_t Ny = 1024 * scale;

    std::vector<float> myFloats(Nx * Ny, 0.1f * rank);
    std::vector<int> myInts(Nx * Ny, (int)(1 + rank));

    hdf5IO.SetParameter("IdleH5Writer",
                        "true"); // set this if not all ranks are writting

    adios2::Variable<float> h5Floats = hdf5IO.DefineVariable<float>(
        "h5Floats", {size * Nx, Ny}, {rank * Nx, 0}, {Nx, Ny}, adios2::ConstantDims);

    adios2::Variable<int> h5Ints = hdf5IO.DefineVariable<int>(
        "h5Ints", {size * Nx, Ny}, {rank * Nx, 0}, {Nx, Ny}, adios2::ConstantDims);

    /** Engine derived class, spawned to start IO operations */
    adios2::Engine hdf5Writer = hdf5IO.Open(testFileName, adios2::Mode::Write);

    int nsteps = 5;

    if (size % 2 == 0)
    {
        // all Ranks must call Put
        /** Write variable for buffering */
        for (int i = 0; i < nsteps; i++)
        {
            hdf5Writer.BeginStep();
            hdf5Writer.Put<float>(h5Floats, myFloats.data());
            hdf5Writer.Put(h5Ints, myInts.data());
            hdf5Writer.EndStep();
        }
    }
    else
    {
        // using collective Begin/EndStep() to run the
        // collective HDF5 calls. Now Ranks can skip writting if no data
        // presented
        for (int i = 0; i < nsteps; i++)
        {
            hdf5Writer.BeginStep();
            if (rank == 0)
            {
                hdf5Writer.Put<float>(h5Floats, myFloats.data());
                hdf5Writer.Put(h5Ints, myInts.data());
            }
            hdf5Writer.EndStep();
        }
    }
    std::vector<int64_t> m_globalDims = {10, 20, 30, 40};
    hdf5IO.DefineAttribute<std::string>("adios2_schema/version_major",
                                        std::to_string(ADIOS2_VERSION_MAJOR));
    hdf5IO.DefineAttribute<std::string>("adios2_schema/version_minor",
                                        std::to_string(ADIOS2_VERSION_MINOR));
    hdf5IO.DefineAttribute<std::string>("/adios2_schema/mesh/type", "explicit");
    hdf5IO.DefineAttribute<std::int64_t>("adios2_schema/mesh/dimension0", m_globalDims[0]);
    hdf5IO.DefineAttribute<std::int64_t>("adios2_schema/mesh/dimension1", m_globalDims[1]);
    hdf5IO.DefineAttribute<std::int64_t>("adios2_schema/mesh/dimension2", m_globalDims[2]);
    hdf5IO.DefineAttribute<std::int64_t>("adios2_schema/mesh/dimension3", m_globalDims[3]);
    hdf5IO.DefineAttribute<std::int64_t>("adios2_schema/mesh/dimension-num", m_globalDims.size());

    hdf5Writer.Close();
}

template <class T>
void ReadVarData(adios2::IO h5IO, adios2::Engine &h5Reader, const std::string &name)
{
    adios2::Variable<T> var = h5IO.InquireVariable<T>(name);

    if (var)
    {
        int nDims = (int)var.Shape().size();
        size_t totalSize = 1;
        for (int i = 0; i < nDims; i++)
        {
            totalSize *= var.Shape()[i];
        }
        std::vector<T> myValues(totalSize);
        // myFloats.data is pre-allocated
        h5Reader.Get<T>(var, myValues.data(), adios2::Mode::Sync);

        // std::cout << "\tValues of "<<name<<": ";
        std::cout << "\tPeek Values: ";

        if (totalSize < 20)
        { // print all
            for (const auto number : myValues)
            {
                std::cout << number << " ";
            }
        }
        else
        {
            size_t counter = 0;
            for (const auto number : myValues)
            {
                if ((counter < 5) || (counter > totalSize - 5))
                {
                    std::cout << number << " ";
                }
                else if (counter == 5)
                {
                    std::cout << " ......  ";
                }
                counter++;
            }
        }
        std::cout << "\n";
    }
}

void readMe(adios2::IO &h5IO, int rank, int size, const char *fileName)
{
    /** Engine derived class, spawned to start IO operations */
    adios2::Engine h5Reader = h5IO.Open(fileName, adios2::Mode::Read);

    const std::map<std::string, adios2::Params> variables = h5IO.AvailableVariables();

    if (0 == rank)
        std::cout << " Num Vars: " << variables.size() << std::endl;

    for (const auto &variablePair : variables)
    {
        std::cout << "Name: " << variablePair.first;
        std::cout << std::endl;

        for (const auto &parameter : variablePair.second)
        {
            std::cout << "\t" << parameter.first << ": " << parameter.second << "\n";
            if (parameter.second == "double")
            {
                ReadVarData<double>(h5IO, h5Reader, variablePair.first);
            }
            else if (parameter.second == "float")
            {
                ReadVarData<float>(h5IO, h5Reader, variablePair.first);
            }
            else if (parameter.second == "unsigned int")
            {
                ReadVarData<unsigned int>(h5IO, h5Reader, variablePair.first);
            }
            else if (parameter.second == "int")
            {
                ReadVarData<int>(h5IO, h5Reader, variablePair.first);
            }
        }
    } // variables

    const std::map<std::string, adios2::Params> attributes = h5IO.AvailableAttributes();

    if (0 == rank)
        std::cout << "Num Attrs:" << attributes.size() << std::endl;

    for (const auto &attrPair : attributes)
    {
        std::cout << "AttrName: " << attrPair.first;
        std::cout << std::endl;

        for (const auto &parameter : attrPair.second)
        {
            std::cout << "\t" << parameter.first << ": " << parameter.second << "\n";

            if (parameter.second == "double")
            {
                // ReadVarData<double>(h5IO, h5Reader, variablePair.first);
            }
            else if (parameter.second == "float")
            {
                // ReadVarData<float>(h5IO, h5Reader, variablePair.first);
            }
            else if (parameter.second == "unsigned int")
            {
                // ReadVarData<unsigned int>(h5IO, h5Reader,
                // variablePair.first);
            }
            else if (parameter.second == "int")
            {
                // ReadVarData<int>(h5IO, h5Reader, variablePair.first);
            }
            //... add more types if needed
        }
    }

    h5Reader.Close();
}

int main(int argc, char *argv[])
{
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    if (provided < MPI_THREAD_MULTIPLE)
    {
        std::cout << "MPI_THREAD_MULTIPLE is not supported, not able to use "
                     "the subfile feature in HDF5. Aborting. \n"
                  << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    try
    {
        /** ADIOS class factory of IO class objects */
        adios2::ADIOS adios(MPI_COMM_WORLD);
        adios2::IO writerIO = adios.DeclareIO("HDFFileIOWriter");

        std::string testName = "test.h5";
        if (argc > 1)
            testName = argv[1];

        if (1 == testName.size())
        {
            writerIO.SetEngine("NullCore");
            writeMe(writerIO, rank, size, "null.bp");
        }
        else
        {
            writeMe(writerIO, rank, size, testName.c_str());
        }

        // read back  if required
        if (argc > 2)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            adios2::IO readerIO = adios.DeclareIO("HDFFileIOReader");
            readMe(readerIO, rank, size, testName.c_str());
        }
    }
    catch (std::invalid_argument &e)
    {
        std::cout << "Invalid argument exception, STOPPING PROGRAM from rank " << rank << "\n";
        std::cout << e.what() << "\n";
    }
    catch (std::ios_base::failure &e)
    {
        std::cout << "IO System base failure exception, STOPPING PROGRAM from rank " << rank
                  << "\n";
        std::cout << e.what() << "\n";
    }
    catch (std::exception &e)
    {
        std::cout << "Exception, STOPPING PROGRAM from rank " << rank << "\n";
        std::cout << e.what() << "\n";
    }

    MPI_Finalize();

    return 0;
}

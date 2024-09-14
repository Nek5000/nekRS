#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <ios>       //std::ios_base::failure
#include <iostream>  //std::cout
#include <stdexcept> //std::invalid_argument std::exception
#include <vector>

#include <adios2.h>

const std::string streamname = "skeleton_stream";

static void printDataStep(const float *data, const size_t start, const size_t count, const int rank,
                          const size_t step)
{
    std::ofstream myfile;
    std::string filename = "data." + std::to_string(rank);
    if (step == 0)
    {
        myfile.open(filename);
    }
    else
    {
        myfile.open(filename, std::ios::app);
    }

    myfile << "rank=" << rank << " size=" << count << " offsets=" << start << " step=" << step
           << std::endl;

    myfile << " step   row   columns " << start << "..." << start + count - 1 << std::endl;
    myfile << "        ";
    for (size_t j = 0; j < count; j++)
    {
        myfile << std::setw(9) << start + j;
    }
    myfile << std::endl;
    myfile << "------------------------------------------------------------"
              "--\n";
    for (size_t i = 0; i < count; i++)
    {
        myfile << std::setw(5) << step << std::setw(5) << start + i;
        for (size_t j = 0; j < count; j++)
        {
            myfile << std::setw(9) << std::setprecision(4) << data[i * count + j];
        }
        myfile << std::endl;
    }
    myfile.close();
}

int main(int argc, char *argv[])
{
    int rank = 0, nproc = 1;
    int retval = 0;

    std::ofstream out("TestSkeletonReaderOutput.txt");
    auto coutbuf = std::cout.rdbuf(out.rdbuf()); // save and redirect

#if ADIOS2_USE_MPI
    int wrank = 0, wnproc = 1;
    MPI_Comm mpiReaderComm;
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wnproc);

    const unsigned int color = 2;
    MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &mpiReaderComm);

    MPI_Comm_rank(mpiReaderComm, &rank);
    MPI_Comm_size(mpiReaderComm, &nproc);
#endif

    try
    {

/** ADIOS class factory of IO class objects */
#if ADIOS2_USE_MPI
        adios2::ADIOS adios(mpiReaderComm);
#else
        adios2::ADIOS adios;
#endif
        adios2::IO io = adios.DeclareIO("reader");
        io.SetEngine("Skeleton");
        io.SetParameter("VerBose", "5");
        adios2::Engine reader = io.Open(streamname, adios2::Mode::Read);

        while (true)
        {
            adios2::StepStatus status = reader.BeginStep(adios2::StepMode::Read, 60.0f);
            if (status != adios2::StepStatus::OK)
            {
                break;
            }

            // discover in the metadata that the variable exists
            adios2::Variable<float> vMyArray;
            size_t gndx;
            vMyArray = io.InquireVariable<float>("myArray");
            if (!vMyArray)
            {
                std::cout << "Missing 'myArray' variable. The Skeleton reader "
                             "engine must retrieve variables from the writer and "
                             "create Variable objects before they can be "
                             "inquired\n";
                // Let's fake the read from now on
                // so that we can test the rest of the read API
                gndx = (size_t)nproc;
                adios2::Variable<float> varArray = io.DefineVariable<float>(
                    "myArray", {gndx}, {gndx / (size_t)nproc}, {gndx / (size_t)nproc});
                (void)varArray;

                adios2::Variable<std::string> varSyncString =
                    io.DefineVariable<std::string>("mySyncString");
                (void)varSyncString;
            }
            else
            {
                // Get the variable dimension
                gndx = vMyArray.Shape()[0];
            }
            size_t ndx = gndx / (size_t)nproc;
            size_t offsx = ndx * (size_t)rank;
            if (rank == nproc - 1)
            {
                // right-most processes need to read all the rest
                ndx = gndx - ndx * (size_t)(nproc - 1);
            }
            size_t step = reader.CurrentStep();
            adios2::Dims count, start;
            count.push_back(ndx);
            start.push_back(offsx);

            if (vMyArray)
            {
                vMyArray.SetSelection({start, count});
            }

            std::vector<float> myArray(ndx);
            std::string s;
            reader.Get("mySyncString", s, adios2::Mode::Sync);
            reader.Get("myArray", myArray.data());
            reader.PerformGets();
            printDataStep(myArray.data(), offsx, ndx, rank, step);

            reader.EndStep();
        }
        reader.Close();
    }
    catch (std::invalid_argument &e)
    {
        std::cout << "Invalid argument exception, STOPPING PROGRAM from rank " << rank << "\n";
        std::cout << e.what() << "\n";
        retval = 1;
    }
    catch (std::ios_base::failure &e)
    {
        std::cout << "IO System base failure exception, STOPPING PROGRAM "
                     "from rank "
                  << rank << "\n";
        std::cout << e.what() << "\n";
        retval = 2;
    }
    catch (std::exception &e)
    {
        std::cout << "Exception, STOPPING PROGRAM from rank " << rank << "\n";
        std::cout << e.what() << "\n";
        retval = 3;
    }

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    std::cout.rdbuf(coutbuf); // reset to standard output again

    return retval;
}

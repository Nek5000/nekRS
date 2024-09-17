#include "adios2.h"
#include <mpi.h>

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// touched block ids are printed.
void queryIDs(adios2::IO &queryIO, std::string &dataFileName, std::string &queryFile)
{
    adios2::Engine reader = queryIO.Open(dataFileName, adios2::Mode::Read, MPI_COMM_WORLD);
    // adios2::QueryWorker* worker = NULL;
    queryIO.SetParameter("StreamReader", "true");
    std::vector<size_t> touched_blockIDs;

    while (reader.BeginStep() == adios2::StepStatus::OK)
    {
        adios2::QueryWorker w = adios2::QueryWorker(queryFile, reader);
        w.GetResultCoverage(touched_blockIDs);

        std::cout << " Num touched blocks =" << touched_blockIDs.size() << std::endl;
        for (auto n : touched_blockIDs)
        {
            std::cout << "\t[" << n << "] " << std::endl;
        }

        reader.EndStep();
    }
    reader.Close();
}

void queryWithStreaming(adios2::IO &queryIO, std::string &dataFileName, std::string &queryFile)
{
    adios2::Engine reader = queryIO.Open(dataFileName, adios2::Mode::Read, MPI_COMM_WORLD);
    // adios2::QueryWorker* worker = NULL;
    queryIO.SetParameter("StreamReader", "true");
    std::vector<adios2::Box<adios2::Dims>> touched_blocks;

    while (reader.BeginStep() == adios2::StepStatus::OK)
    {
        adios2::QueryWorker w = adios2::QueryWorker(queryFile, reader);
        w.GetResultCoverage(touched_blocks);

        std::cout << " Num touched regions ="
                  << touched_blocks.size()
                  // std::cout << " ... now can read out touched blocks ... size=" <<
                  // touched_blocks.size()
                  << std::endl;
        for (auto n : touched_blocks)
        {
            std::ostringstream startStr;
            std::ostringstream countStr;
            for (size_t k = 0; k < n.first.size(); k++)
            {
                startStr << n.first[k] << " ";
                countStr << n.second[k] << " ";
            }
            std::cout << "\t[" << startStr.str() << "]  [" << countStr.str() << "]" << std::endl;
        }
        reader.EndStep();
    }
    reader.Close();
}

int main(int argc, char *argv[])
{
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    try
    {
        std::string configFileName = "query.xml";
        std::string dataFileName = "/tmp/heatbp4.bp";
        if (argc <= 2)
        {
            std::cout << "Usage: " << argv[0] << " configFileName  dataFilePath (queryFile)"
                      << std::endl;
            std::cout << "    e.g.  " << argv[0] << " bp4io.xml heat_bp4.bp/ " << std::endl;
            std::cout << "    or    " << argv[0] << " bp4io.xml heat_bp4.bp/ q1.json" << std::endl;
            return 0;
        }

        configFileName = argv[1];
        dataFileName = argv[2];

        adios2::ADIOS ad = adios2::ADIOS(configFileName, MPI_COMM_WORLD);

        adios2::IO queryIO = ad.DeclareIO("query");

        std::string queryFile = configFileName;
        if (argc > 3)
        {
            queryFile = argv[3];
        }
        if (rank == 0)
        {
            std::cout << " using config file = " << configFileName << std::endl;
            std::cout << "         data file = " << dataFileName << std::endl;
            std::cout << "         queryfile = " << queryFile << std::endl;
        }

        queryIDs(queryIO, dataFileName, queryFile);

        std::cout << "\n" << std::endl;
        queryWithStreaming(queryIO, dataFileName, queryFile);

        return 0;
    }
    catch (std::exception &e)
    {
        std::cout << "Exception, STOPPING PROGRAM from rank " << rank << "\n";
        std::cout << e.what() << "\n";
    }

    MPI_Finalize();

    return 0;
}

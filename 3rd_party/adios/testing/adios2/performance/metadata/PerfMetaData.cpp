/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cctype>
#include <cstdint>
#include <cstring>
#include <ctime>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <stdexcept>
#include <thread>

#include <adios2.h>

#ifndef _WIN32
#include "strings.h"
#else
#define strcasecmp _stricmp
#endif

std::string engine = "sst";
adios2::Params engineParams = {}; // parsed from command line

int NSteps = 10;
int NumVars = 100;
int NumArrays = 100;
int NumAttrs = 100;
int NumBlocks = 1;
int ReaderDelay = 0;
int WriterSize;
int FileMode = 0;
int FileWriter = 0;
int FileReader = 0;
int UseMPI = 1;
int E3sm = 0;
int Warpx = 0;
char *ErrorStr = NULL;

int LastVarSize;
int LastArraySize;
int LastAttrsSize;

typedef enum
{
    WriterConsolidation,
    ReaderInstallation,
    ReaderTraversal
} TestModeEnum;

TestModeEnum TestMode = WriterConsolidation;

static std::string Trim(std::string &str)
{
    size_t first = str.find_first_not_of(' ');
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

/*
 * Engine parameters spec is a poor-man's JSON.  name:value pairs are separated
 * by equal.  White space is trimmed off front and back.  No quotes or anything
 * fancy allowed.
 */
static adios2::Params ParseEngineParams(std::string Input)
{
    std::istringstream ss(Input);
    std::string Param;
    adios2::Params Ret = {};

    while (std::getline(ss, Param, ','))
    {
        std::istringstream ss2(Param);
        std::string ParamName;
        std::string ParamValue;
        std::getline(ss2, ParamName, '=');
        if (!std::getline(ss2, ParamValue, '='))
        {
            throw std::invalid_argument("Engine parameter \"" + Param + "\" missing value");
        }
        Ret[Trim(ParamName)] = Trim(ParamValue);
    }
    return Ret;
}

static void Usage()
{
    std::cout << "PerfMetaData <opt args> " << std::endl;
    std::cout << "  --num_steps <steps>" << std::endl;
    std::cout << "  --num_vars <vars>" << std::endl;
    std::cout << "  --num_arrays <arrays>" << std::endl;
    std::cout << "  --num_attrs <attributes>" << std::endl;
    std::cout << "  --e3sm   (Approx E3SM characteristics)" << std::endl;
    std::cout << "  --warpx   (Approx WarpX characteristics)" << std::endl;
    std::cout << "  --engine <enginename>" << std::endl;
    std::cout << "  --engine_params <param=value,param=value>" << std::endl;
    std::cout << "  --file  (run in file mode)" << std::endl;
}

static void ParseArgs(int argc, char **argv)
{
    int bare_arg = 0;
    while (argc > 1)
    {
        if (std::string(argv[1]) == "--file_reader")
        {
            FileMode = 1;
            FileReader = 1;
        }
        else if (std::string(argv[1]) == "--file_writer")
        {
            FileMode = 1;
            FileWriter = 1;
        }
        else if (std::string(argv[1]) == "--no_mpi")
        {
            UseMPI = 0;
        }
        else if (std::string(argv[1]) == "--num_steps")
        {
            std::istringstream ss(argv[2]);
            if (!(ss >> NSteps))
                std::cerr << "Invalid number for num_steps " << argv[1] << '\n';
            argv++;
            argc--;
        }
        else if (std::string(argv[1]) == "--num_vars")
        {
            std::istringstream ss(argv[2]);
            if (!(ss >> NumVars))
                std::cerr << "Invalid number for number of variables " << argv[1] << '\n';
            argv++;
            argc--;
        }
        else if (std::string(argv[1]) == "--num_arrays")
        {
            std::istringstream ss(argv[2]);
            if (!(ss >> NumArrays))
                std::cerr << "Invalid number for number of arrays" << argv[1] << '\n';
            argv++;
            argc--;
        }
        else if (std::string(argv[1]) == "--num_attrs")
        {
            std::istringstream ss(argv[2]);
            if (!(ss >> NumAttrs))
                std::cerr << "Invalid number for number of attrs" << argv[1] << '\n';
            argv++;
            argc--;
        }
        else if (std::string(argv[1]) == "--num_blocks")
        {
            std::istringstream ss(argv[2]);
            if (!(ss >> NumBlocks))
                std::cerr << "Invalid number for number of blocks" << argv[1] << '\n';
            argv++;
            argc--;
        }
        else if (std::string(argv[1]) == "--e3sm")
        {
            NumArrays = 535;
            NumVars = 1000;
            NumAttrs = 3800;
            NumBlocks = 1;
            E3sm = 1;
        }
        else if (std::string(argv[1]) == "--warpx")
        {
            NumArrays = 24;
            NumVars = 1000;
            NumAttrs = 0;
            NumBlocks = 86;
            Warpx = 1;
        }
        else if (std::string(argv[1]) == "--test_mode")
        {
            std::string mode = argv[2];
            std::transform(mode.begin(), mode.end(), mode.begin(),
                           [](unsigned char c) { return std::tolower(c); });
            if (mode == "writerconsolidation")
            {
                TestMode = WriterConsolidation;
            }
            else if (mode == "readerinstallation")
            {
                TestMode = ReaderInstallation;
            }
            else if (mode == "readertraversal")
            {
                TestMode = ReaderTraversal;
            }
            else
            {
                std::cerr << "Unknown test mode \"" << argv[2]
                          << "\", must be one of WriterConsolidation, "
                             "ReaderInstallation, or ReaderTraversal"
                          << std::endl;
            }
            argv++;
            argc--;
        }
        else if (std::string(argv[1]) == "--engine_params")
        {
            engineParams = ParseEngineParams(argv[2]);
            argv++;
            argc--;
        }
        else if (std::string(argv[1]) == "--engine")
        {
            engine = std::string(argv[2]);
            argv++;
            argc--;
        }
        else if (std::string(argv[1]) == "--reader_delay")
        {
            std::istringstream ss(argv[2]);
            if (!(ss >> ReaderDelay))
                std::cerr << "Invalid number for ms_delay " << argv[1] << '\n';
            argv++;
            argc--;
        }
        else
        {
            if (std::string(argv[1], 2) == "--")
            {
                ErrorStr = argv[1];
            }
            if (bare_arg == 0)
            {
                /* first arg without -- is engine */
                engine = std::string(argv[1]);
                bare_arg++;
            }
            else if (bare_arg == 1)
            {
                engineParams = ParseEngineParams(argv[1]);
                bare_arg++;
            }
            else
            {

                throw std::invalid_argument("Unknown argument \"" + std::string(argv[1]) + "\"");
            }
        }
        argv++;
        argc--;
    }
}

MPI_Comm testComm;
std::chrono::duration<double> elapsed;
std::chrono::duration<double> InstallTime, TraversalTime;

void DoWriter(adios2::Params writerParams)
{
    // form a mpiSize * Nx 1D array
    int mpiRank = 0, mpiSize = 1;

#if ADIOS2_USE_MPI
    if (UseMPI)
    {
        MPI_Comm_rank(testComm, &mpiRank);
        MPI_Comm_size(testComm, &mpiSize);
    }
#endif
    if ((mpiRank == 0) && E3sm)
        std::cout << "E3SM mode.  Full run: 960 Timesteps, at 1344 ranks" << std::endl;
    if ((mpiRank == 0) && Warpx)
        std::cout << "Warpx mode.  Full run: 50 Timesteps, at 3000 ranks" << std::endl;
    if ((mpiRank == 0) && ErrorStr)
    {
        std::cout << "Unknown arg " << ErrorStr << std::endl;
        Usage();
        exit(1);
    }

    // Write test data using ADIOS2

    adios2::ADIOS *adios;
    if (UseMPI)
        adios = new adios2::ADIOS(testComm);
    else
        adios = new adios2::ADIOS();

    adios2::IO io = adios->DeclareIO("TestIO");

    io.SetEngine(engine);
    io.SetParameters(writerParams);
    adios2::Engine writer = io.Open("MetaDataTest", adios2::Mode::Write);
    std::chrono::time_point<std::chrono::high_resolution_clock> start, finish;
    std::vector<float> myFloats = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    adios2::Variable<float> *Floats = new adios2::Variable<float>[NumVars];
    adios2::Variable<float> *FloatArrays = new adios2::Variable<float>[NumArrays];
    if (mpiRank == 0)
    {
        // attributes and globals on rank 0
        for (int i = 0; i < NumVars; i++)
        {
            std::string varname = "Variable" + std::to_string(i);
            Floats[i] = io.DefineVariable<float>(varname);
        }
        for (int i = 0; i < NumAttrs; i++)
        {
            std::string varname = "Attribute" + std::to_string(i);
            io.DefineAttribute<float>(varname, 0.0);
        }
    }
    for (int i = 0; i < NumArrays; i++)
    {
        std::string varname = "Array" + std::to_string(i);
        if (NumBlocks == 1)
        {
            FloatArrays[i] =
                io.DefineVariable<float>(varname, {(unsigned long)mpiSize},
                                         {(unsigned long)mpiRank}, {1}, adios2::ConstantDims);
        }
        else
        {
            FloatArrays[i] = io.DefineVariable<float>(varname, {}, {}, {1}, adios2::ConstantDims);
        }
    }
    start = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < NSteps; j++)
    {
        writer.BeginStep();
        if (mpiRank == 0)
        {
            // globals on rank 0
            for (int i = 0; i < NumVars; i++)
            {
                writer.Put(Floats[i], myFloats.data()[0]);
            }
        }
        for (int i = 0; i < NumArrays; i++)
        {
            for (int j = 0; j < NumBlocks; j++)
            {
                writer.Put(FloatArrays[i], myFloats.data());
            }
        }
        writer.EndStep();
    }
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    writer.Close();
}

void DoReader()
{
    // form a mpiSize * Nx 1D array
    int mpiRank = 0, mpiSize = 1;

#if ADIOS2_USE_MPI
    if (UseMPI)
    {
        MPI_Comm_rank(testComm, &mpiRank);
        MPI_Comm_size(testComm, &mpiSize);
    }
#endif

    adios2::ADIOS *adios;
    if (UseMPI)
        adios = new adios2::ADIOS(testComm);
    else
        adios = new adios2::ADIOS();

    adios2::IO io = adios->DeclareIO("TestIO");

    io.SetEngine(engine);
    engineParams["SpeculativePreloadMode"] = "Off";
    engineParams["ReaderShortCircuitReads"] = "On";
    io.SetParameters(engineParams);
    adios2::Engine reader = io.Open("MetaDataTest", adios2::Mode::Read);
    std::chrono::time_point<std::chrono::high_resolution_clock> startTS, endBeginStep, finishTS;
    std::vector<float> myFloats = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::this_thread::sleep_for(std::chrono::seconds(ReaderDelay));
    std::vector<float> in(1);

    in.resize(WriterSize);
    while (1)
    {
        startTS = std::chrono::high_resolution_clock::now();
        adios2::StepStatus status = reader.BeginStep();
        endBeginStep = std::chrono::high_resolution_clock::now();
        if (status != adios2::StepStatus::OK)
        {
            break;
        }
        //	if (DoGets) {
        std::vector<adios2::Variable<float>> Floats;
        std::vector<adios2::Variable<float>> FloatArrays;
        std::vector<adios2::Attribute<float>> Attributes;
        int Cont = 1;
        int i = 0;
        while (Cont)
        {
            std::string varname = "Variable" + std::to_string(i++);
            adios2::Variable<float> tmp = io.InquireVariable<float>(varname);
            if (tmp)
            {
                Floats.push_back(tmp);
            }
            else
            {
                Cont = 0;
            }
        }
        Cont = 1;
        i = 0;
        while (Cont)
        {
            std::string varname = "Array" + std::to_string(i++);
            adios2::Variable<float> tmp = io.InquireVariable<float>(varname);
            if (tmp)
            {
                FloatArrays.push_back(tmp);
            }
            else
            {
                Cont = 0;
            }
        }
        i = 0;
        Cont = 1;
        while (Cont)
        {
            std::string varname = "Attribute" + std::to_string(i++);
            adios2::Attribute<float> tmp = io.InquireAttribute<float>(varname);
            if (tmp)
            {
                Attributes.push_back(tmp);
            }
            else
            {
                Cont = 0;
            }
        }
        LastVarSize = (int)Floats.size();
        for (auto Var : Floats)
        {
            reader.Get(Var, in.data());
        }
        LastArraySize = (int)FloatArrays.size();
        for (auto Var : FloatArrays)
        {
            if (Var.ShapeID() == adios2::ShapeID::GlobalArray)
            {
                reader.Get(Var, in.data());
            }
            else
            {
                // local, go through blocks
                for (int rank = 0; rank < WriterSize; rank++)
                {
                    for (auto blk : reader.BlocksInfo(Var, reader.CurrentStep()))
                    {
                        Var.SetBlockSelection(blk.BlockID);
                        reader.Get(Var, in.data());
                    }
                }
            }
        }
        LastAttrsSize = (int)Attributes.size();
        for (auto Attr : Attributes)
        {
            if (Attr.Data().front() != 0.0)
            {
                std::cerr << "Bad attr data" << std::endl;
            }
        }
        //	}
        reader.EndStep();
        finishTS = std::chrono::high_resolution_clock::now();
        InstallTime += (endBeginStep - startTS);
        TraversalTime += (finishTS - endBeginStep);
    }
    reader.Close();
}

void DoReaderOutput()
{
    std::cout << "Metadata Installation Time " << InstallTime.count() << " seconds." << std::endl;

    std::cout << "Metadata Traversal Time " << TraversalTime.count() << " seconds." << std::endl;

    std::cout << "Parameters Nsteps=" << NSteps << ", NumArrays=" << NumArrays
              << ", NumVArs=" << NumVars << ", NumAttrs=" << NumAttrs << ", NumBlocks=" << NumBlocks
              << std::endl;
    if ((NumArrays != LastArraySize) || (NumVars != LastVarSize) || (NumAttrs != LastAttrsSize))
    {
        std::cout << "Arrays=" << LastArraySize << ", Vars=" << LastVarSize
                  << ", Attrs=" << LastAttrsSize << ", NumBlocks=" << NumBlocks << std::endl;
        std::cout << "Inconsistency" << std::endl;
    }
}

int main(int argc, char **argv)
{
    int key;
    std::string MeasurementString;
    unsigned int color = 0;

    ParseArgs(argc, argv);

    if (UseMPI)
    {
        int provided;

        // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
        MPI_Init_thread(nullptr, nullptr, MPI_THREAD_MULTIPLE, &provided);

        MPI_Comm_rank(MPI_COMM_WORLD, &key);
        MPI_Comm_size(MPI_COMM_WORLD, &WriterSize);
        if (!FileMode)
        {
            if (key > 0)
            {
                color = 1;
                WriterSize--;
            }
        }
        MPI_Comm_split(MPI_COMM_WORLD, color, key, &testComm);
    }

    if (FileMode)
    {
        if (FileWriter)
        {
            DoWriter(engineParams);
            if (key == 0)
                std::cout << "File Writer Time " << elapsed.count() << " seconds." << std::endl;
        }
        else
        {
            DoReader();
            if (key == 0)
                DoReaderOutput();
        }
        if (UseMPI)
        {
            MPI_Finalize();
        }
        exit(0);
    }
    /*
     * Below here for staging solutions
     * Always run with MPI, 1 or more writers running with 1 reader
     */
    if (WriterSize < 2)
    {
        std::cerr << "PerfMetaData cannot run with MPI size < 2" << std::endl;
        exit(1);
    }
    if ((NumBlocks > 1) && (key == 0))
    {
        std::cerr << "Warning, metadata info for FFS not valid for num_blocks > 1" << std::endl;
    }
    // first all writer ranks do Writer calcs with no reader
    if (key > 0)
    {
        adios2::Params ConsolidationParams = engineParams; // parsed from command line
        ConsolidationParams["RendezvousReaderCount"] = "0";
        DoWriter(ConsolidationParams);
        if (key == 1)
        {
            std::cout << "Metadata Consolidation Time " << elapsed.count() << " seconds."
                      << std::endl;
        }
        ReaderDelay = (int)elapsed.count() + 5;
    }

    MPI_Bcast(&ReaderDelay, 1, MPI_INT, 1, MPI_COMM_WORLD);
    // next do writers and single reader
    if (key == 0)
    {
        DoReader();
    }
    else
    {
        DoWriter(engineParams);
    }

    if (key == 0)
    {
        DoReaderOutput();
    }
    MPI_Finalize();

    return 0;
}

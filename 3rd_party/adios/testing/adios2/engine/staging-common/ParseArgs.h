#ifndef _WIN32
#include "strings.h"
#else
#define strcasecmp _stricmp
#endif
std::string fname = "ADIOS2Common";
std::string engine = "sst";
adios2::Params engineParams = {}; // parsed from command line

#ifndef TESTING_ADIOS2_ENGINE_COMMON_TESTDATA_H_
// Usually we get this from TestData.h, but not needed everywhere.
std::size_t Nx = 10;
#endif
bool SharedIO = false;
bool SharedVar = false;
int NSteps = 10;
int DurationSeconds = 60 * 60 * 24 * 365; // one year default
int DelayMS = 1000;                       // one step per sec default
int Latest = 0;
int Discard = 0;
int IncreasingDelay = 0;
int NonBlockingBeginStep = 0;
int CompressSz = 0;
int CompressZfp = 0;
int TimeGapExpected = 0;
int IgnoreTimeGap = 1;
int ExpectWriterFailure = 0;
int ExpectOpenTimeout = 0;
int ZeroDataVar = 0;
int ZeroDataRank = 0;
int DelayWhileHoldingStep = 0;
int LongFirstDelay = 0;
int FirstTimestepMustBeZero = 0;
int LockGeometry = 0;
bool VaryingDataSize = false;
bool TestVarDestruction = false;
bool AdvancingAttrs = false;
int NoData = 0;
int NoDataNode = -1;
int Flush = 0;
int EarlyExit = 0;
int LocalCount = 1;
int DataSize = 5 * 1024 * 1024 / 8; /* DefaultMinDeferredSize is 4*1024*1024
                                       This should be more than that. */
bool ModifiableAttributes = false;
bool RoundRobin = false;
bool OnDemand = false;
bool DontClose = false;

std::string shutdown_name = "DieTest";
adios2::Mode GlobalWriteMode = adios2::Mode::Deferred;
adios2::Mode GlobalReadMode = adios2::Mode::Deferred;

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

void ParseArgs(int argc, char **argv)
{
    int bare_arg = 0;
    while (argc > 1)
    {
        if (std::string(argv[1]) == "--num_steps")
        {
            std::istringstream ss(argv[2]);
            if (!(ss >> NSteps))
                std::cerr << "Invalid number for num_steps " << argv[1] << '\n';
            argv++;
            argc--;
        }
        else if (std::string(argv[1]) == "--nx")
        {
            std::istringstream ss(argv[2]);
            if (!(ss >> Nx))
                std::cerr << "Invalid number for nx (base element count) " << argv[1] << '\n';
            argv++;
            argc--;
        }
        else if (std::string(argv[1]) == "--expect_time_gap")
        {

            TimeGapExpected++;
            IgnoreTimeGap = 0;
        }
        else if (std::string(argv[1]) == "--expect_writer_failure")
        {
            IncreasingDelay = 1;
            ExpectWriterFailure++;
        }
        else if (std::string(argv[1]) == "--non_blocking")
        {
            IncreasingDelay = 1;
            NonBlockingBeginStep = 1;
        }
        else if (std::string(argv[1]) == "--modifiable_attributes")
        {
            ModifiableAttributes = true;
        }
        else if (std::string(argv[1]) == "--expect_contiguous_time")
        {
            TimeGapExpected = 0;
            IgnoreTimeGap = 0;
        }
        else if (std::string(argv[1]) == "--discard")
        {
            IncreasingDelay = 1;
            Discard = 1;
        }
        else if (std::string(argv[1]) == "--delay_while_holding")
        {
            DelayWhileHoldingStep = 1;
        }
        else if (std::string(argv[1]) == "--precious_first")
        {
            FirstTimestepMustBeZero = 1;
        }
        else if (std::string(argv[1]) == "--lock_geometry")
        {
            LockGeometry = 1;
        }
        else if (std::string(argv[1]) == "--ignore_time_gap")
        {
            IgnoreTimeGap++;
        }
        else if (std::string(argv[1]) == "--compress_sz")
        {
            CompressSz++;
        }
        else if (std::string(argv[1]) == "--shared_io")
        {
            SharedIO = true;
        }
        else if (std::string(argv[1]) == "--var_destruction")
        {
            TestVarDestruction = true;
        }
        else if (std::string(argv[1]) == "--shared_var")
        {
            SharedVar = true;
            SharedIO = true;
        }
        else if (std::string(argv[1]) == "--compress_zfp")
        {
            CompressZfp++;
        }
        else if (std::string(argv[1]) == "--filename")
        {
            fname = std::string(argv[2]);
            argv++;
            argc--;
        }
        else if (std::string(argv[1]) == "--write_mode")
        {
            if (strcasecmp(argv[2], "sync") == 0)
            {
                GlobalWriteMode = adios2::Mode::Sync;
            }
            else if (strcasecmp(argv[2], "deferred") == 0)
            {
                GlobalWriteMode = adios2::Mode::Deferred;
            }
            else
            {
                std::cerr << "Invalid mode for --write_mode " << argv[2] << std::endl;
            }
            argv++;
            argc--;
        }
        else if (std::string(argv[1]) == "--read_mode")
        {
            if (strcasecmp(argv[2], "sync") == 0)
            {
                GlobalReadMode = adios2::Mode::Sync;
            }
            else if (strcasecmp(argv[2], "deferred") == 0)
            {
                GlobalReadMode = adios2::Mode::Deferred;
            }
            else
            {
                std::cerr << "Invalid mode for --write_mode " << argv[2] << std::endl;
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
        else if (std::string(argv[1]) == "--expect_timeout")
        {
            ExpectOpenTimeout++;
        }
        else if (std::string(argv[1]) == "--duration")
        {
            std::istringstream ss(argv[2]);
            if (!(ss >> DurationSeconds))
                std::cerr << "Invalid number for duration " << argv[1] << '\n';
            argv++;
            argc--;
        }
        else if (std::string(argv[1]) == "--shutdown_filename")
        {
            shutdown_name = std::string(argv[2]);
            argv++;
            argc--;
        }
        else if (std::string(argv[1]) == "--latest")
        {
            IncreasingDelay = 1;
            Latest = 1;
        }
        else if (std::string(argv[1]) == "--varying_data_size")
        {
            VaryingDataSize = true;
        }
        else if (std::string(argv[1]) == "--advancing_attributes")
        {
            AdvancingAttrs = true;
        }
        else if (std::string(argv[1]) == "--long_first_delay")
        {
            LongFirstDelay = 1;
        }
        else if (std::string(argv[1]) == "--ms_delay")
        {
            std::istringstream ss(argv[2]);
            if (!(ss >> DelayMS))
                std::cerr << "Invalid number for ms_delay " << argv[1] << '\n';
            argv++;
            argc--;
        }
        else if (std::string(argv[1]) == "--zero_data_var")
        {
            ZeroDataVar++;
        }
        else if (std::string(argv[1]) == "--round_robin")
        {
            if (OnDemand)
                std::cerr << "OnDemand already specified, round robin ignored" << std::endl;
            RoundRobin = true;
        }
        else if (std::string(argv[1]) == "--on_demand")
        {
            if (RoundRobin)
                std::cerr << "RoundRobin already specified, on_demand ignored" << std::endl;
            OnDemand = true;
        }
        else if (std::string(argv[1]) == "--no_data")
        {
            NoData++;
        }
        else if (std::string(argv[1]) == "--no_data_node")
        {
            std::istringstream ss(argv[2]);
            if (!(ss >> NoDataNode))
                std::cerr << "Invalid number for --no_data_node argument" << argv[1] << '\n';
            argv++;
            argc--;
        }
        else if (std::string(argv[1]) == "--local_count")
        {
            std::istringstream ss(argv[2]);
            if (!(ss >> LocalCount))
                std::cerr << "Invalid number for --local_count argument" << argv[1] << '\n';
            argv++;
            argc--;
        }
        else if (std::string(argv[1]) == "--data_size")
        {
            std::istringstream ss(argv[2]);
            if (!(ss >> DataSize))
                std::cerr << "Invalid number for --data_size argument" << argv[1] << '\n';
            argv++;
            argc--;
        }
        else if (std::string(argv[1]) == "--early_exit")
        {
            EarlyExit++;
        }
        else if (std::string(argv[1]) == "--flush")
        {
            Flush++;
        }
        else if (std::string(argv[1]) == "--dont_close")
        {
            DontClose = true;
        }
        else if (std::string(argv[1]) == "--disable_mpmd")
        {
            // someone else should have eaten this arg, but if it gets here,
            // ignore it
        }
        else
        {
            if (bare_arg == 0)
            {
                /* first arg without -- is engine */
                engine = std::string(argv[1]);
                bare_arg++;
            }
            else if (bare_arg == 1)
            {
                /* second arg without -- is filename */
                fname = std::string(argv[1]);
                bare_arg++;
            }
            else if (bare_arg == 2)
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

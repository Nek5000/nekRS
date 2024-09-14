/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*
 * List the content of a BP file. Rewritten from adios 1.x bpls C code with
 *adios 2.x C++ API
 *
 * Author: Norbert Podhorszki, pnorbert@ornl.gov
 *
 **/

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "bpls.h"
#include "verinfo.h"

#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <fstream>
#include <string>
#include <thread>
#include <vector>

#include <errno.h>

#include "adios2/helper/adiosLog.h"

#if defined(__GNUC__) && !(defined(__ICC) || defined(__INTEL_COMPILER))
#if (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__) < 40900
/* pre GCC 4.9 cannot handle the C++ regex implementation. Will use C-lib
 * regex"*/
#define USE_C_REGEX
#endif
#endif

#ifdef USE_C_REGEX
#include <regex.h> // regular expression matching
#else
#include <regex>
#endif

#ifdef _WIN32
#include "shlwapi.h"
#include "windows.h"
#pragma warning(disable : 4101) // unreferenced local variable
#else
#include <fnmatch.h>
#endif

#include "adios2/helper/adiosString.h" // EndsWith
#include "adios2/helper/adiosSystem.h" //isHDF5File
#include <adios2sys/CommandLineArguments.hxx>
#include <adios2sys/SystemTools.hxx>
#include <pugixml.hpp>

namespace adios2
{
namespace utils
{

using EntryMap = std::map<std::string, Entry>;

// global variables
// Values from the arguments or defaults

// output files' starting path (can be extended with subdirs,
// names, indexes)
std::string outpath;
char *varmask[MAX_MASKS];     // can have many -var masks (either shell patterns or
                              // extended regular expressions)
int nmasks;                   // number of masks specified
char *vfile;                  // file name to bpls
std::string start;            // dimension spec starting points
std::string count;            // dimension spec counts
std::string format;           // format string for one data element (e.g. %6.2f)
std::string transport_params; // Transport parameters (e.g. "Library=stdio,verbose=3")
std::string engine_name;      // Engine name (e.g. "BP5")
std::string engine_params;    // Engine parameters (e.g. "SelectSteps=0:5:2")
std::string accuracy_def;     // Accuracy definition (e.g. "accuracy="0.0,0.0,rel")

// Flags from arguments or defaults
bool dump; // dump data not just list info(flag == 1)
bool output_xml;
bool use_regexp;         // use varmasks as regular expressions
bool sortnames;          // sort names before listing
bool listattrs;          // do list attributes too
bool listmeshes;         // do list meshes too
bool attrsonly;          // do list attributes only
bool longopt;            // -l is turned on
bool timestep;           // read step by step
bool ignore_flatten;     // dont flatten steps to one
bool filestream = false; // are we using an engine through FileStream?
bool noindex;            // do no print array indices with data
bool printByteAsChar;    // print 8 bit integer arrays as string
bool plot;               // dump histogram related information
bool hidden_attrs;       // show hidden attrs in BP file
int hidden_attrs_flag;   // to be passed on in option struct
bool show_decomp;        // show decomposition of arrays
bool show_version;       // print binary version info of file before work
adios2::Accuracy accuracy;
bool accuracyWasSet = false;

// other global variables
char *prgname; /* argv[0] */
// long timefrom, timeto;
int64_t istart[MAX_DIMS], icount[MAX_DIMS]; // negative values are allowed
int ndimsspecified = 0;
#ifdef USE_C_REGEX
regex_t varregex[MAX_MASKS]; // compiled regular expressions of varmask
#else
std::vector<std::regex> varregex;
#endif
int ncols = 6; // how many values to print in one row (only for -p)
int verbose = 0;
FILE *outf; // file to print to or stdout
char commentchar;

// help function
void display_help()
{
    // printf( "Usage: %s  \n", prgname);
    printf("usage: bpls [OPTIONS] file [mask1 mask2 ...]\n"
           "\nList/dump content of a BP/HDF5 file. \n"
           "A mask can be a shell pattern like with 'ls' e.g. \"*/x?\".\n"
           "Variables with multiple timesteps are reported with an extra "
           "dimensions.\n"
           "The time dimension is the first dimension then.\n"
           "\n"
           "  --long      | -l           Print values of all scalars and "
           "attributes and\n"
           "                               min/max values of arrays (no overhead "
           "to get them!)\n"
           "  --attrs     | -a           List/match attributes too\n"
           "  --attrsonly | -A           List attributes only\n"
           "  --meshes    | -m           List meshes\n"
           /*
           "  --sort      | -r           Sort names before listing\n"
           */
           "  --timestep  | -t           Read content step by step (stream "
           "reading)\n"
           "  --ignore_flatten           Display steps as written (don't flatten, even if writer "
           "said to)\n"
           "  --dump      | -d           Dump matched variables/attributes\n"
           "                               To match attributes too, add option "
           "-a\n"
           "  --regexp    | -e           Treat masks as extended regular "
           "expressions\n"
           "  --output    | -o <path>    Print to a file instead of stdout\n"
           /*
              "  --xml    | -x            # print as xml instead of ascii text\n"
            */
           "  --start     | -s \"spec\"    Offset indices in each dimension \n"
           "                               (default is 0 for all dimensions) \n"
           "                               <0 is handled as in python (-1 is "
           "last)\n"
           "  --count     | -c \"spec\"    Number of elements in each dimension\n"
           "                               -1 denotes 'until end' of dimension\n"
           "                               (default is -1 for all dimensions)\n"
           "  --noindex   | -y           Print data without array indices\n"
           "  --string    | -S           Print 8bit integer arrays as strings\n"
           "  --columns   | -n \"cols\"    Number of data elements per row to "
           "print\n"
           "  --format    | -f \"str\"     Format string to use for one data item "
           "in print\n"
           "                               instead of the default. E.g. "
           "\"%%6.3f\"\n"
           "  --hidden_attrs             Show hidden ADIOS attributes in the "
           "file\n"
           "  --decomp    | -D           Show decomposition of variables as layed "
           "out in file\n"
           "  --error     | -X string    Specify read accuracy (error,norm,rel|abs)\n"
           "                             e.g. error=\"0.0,0.0,abs\"\n"
           "                             L2 norm = 0.0, Linf = inf\n"

           "  --transport-parameters | -T         Specify File transport "
           "parameters\n"
           "                                      e.g. \"Library=stdio\"\n"
           "  --engine               | -E <name>  Specify ADIOS Engine\n"
           "  --engine-params        | -P string  Specify ADIOS Engine "
           "Parameters\n"
           "                                      e.g. \"SelectSteps=0:n:2\""
           "\n"
           "  Examples for slicing:\n"
           "  -s \"0,0,0\"   -c \"1,99,1\":  Print 100 elements (of the 2nd "
           "dimension).\n"
           "  -s \"0,0\"     -c \"1,-1\":    Print the whole 2nd dimension "
           "however large it is.\n"
           "  -s \"-1,-1\"   -c \"1,1\":     Print the very last element (of a 2D "
           "array)\n"
           "\n"
           "Help options\n"
           "  --help      | -h           Print this help.\n"
           "  --verbose   | -v           Print log about what this program is "
           "doing.\n"
           "                               Use multiple -v to increase logging "
           "level.\n"
           "  --version   | -V           Print version information; compatible "
           " with\n"
           "                               --verbose for additional information, "
           "i.e.\n"
           "                               -v --version.\n"
           "\nTypical use: bpls -lav <file>\n");
}

bool option_help_was_called = false;
int optioncb_help(const char *argument, const char *value, void *call_data)
{
    // adios2sys::CommandLineArguments *arg =
    // static_cast<adios2sys::CommandLineArguments *>(call_data);
    // printf("%s\n", arg->GetHelp());
    display_help();
    option_help_was_called = true;
    return 1;
}

int optioncb_verbose(const char *argument, const char *value, void *call_data)
{
    verbose++;
    return 1;
}

void print_bpls_version()
{
    if (verbose == 0)
    {
        printf(ADIOS2_VERSION_STR "\n");
        option_help_was_called = true;
    }
    else
    {
        printf("blps: ADIOS file introspection utility\n");
        printf("\nBuild configuration:\n");
        if (strlen(ADIOS_INFO_VER_GIT) > 0)
        {
            printf("ADIOS version: %s (%s)\n", ADIOS2_VERSION_STR, ADIOS_INFO_VER_GIT);
        }
        else
        {
            printf("ADIOS version: %s\n", ADIOS2_VERSION_STR);
        }
        if (strlen(ADIOS_INFO_COMPILER_WRAP) > 0)
        {
            printf("C++ Compiler:  %s %s (%s)\n", ADIOS_INFO_COMPILER_ID, ADIOS_INFO_COMPILER_VER,
                   ADIOS_INFO_COMPILER_WRAP);
        }
        else
        {
            printf("C++ Compiler:  %s %s\n", ADIOS_INFO_COMPILER_ID, ADIOS_INFO_COMPILER_VER);
        }
        printf("Target OS:     %s\n", ADIOS_INFO_SYSTEM);
        printf("Target Arch:   %s\n", ADIOS_INFO_ARCH);

        size_t nengines;
        const char *const *list_engines;
        adios2_available_engines(&nengines, &list_engines);
        printf("Available engines = %zu:", nengines);
        for (size_t i = 0; i < nengines; ++i)
        {
            printf(" %s", list_engines[i]);
            if (i < nengines - 1)
            {
                printf(",");
            }
        }
        printf("\n");

        size_t noperators;
        const char *const *list_operators;
        adios2_available_operators(&noperators, &list_operators);
        printf("Available operators = %zu:", noperators);
        for (size_t i = 0; i < noperators; ++i)
        {
            printf(" %s", list_operators[i]);
            if (i < noperators - 1)
            {
                printf(",");
            }
        }
        printf("\n");

        size_t nfeatures;
        const char *const *list_features;
        adios2_available_features(&nfeatures, &list_features);
        printf("Available features = %zu:", nfeatures);
        for (size_t i = 0; i < nfeatures; ++i)
        {
            printf(" %s", list_features[i]);
            if (i < nfeatures - 1)
            {
                printf(",");
            }
        }
        printf("\n");
    }
}

bool introspectAsHDF5File(std::ifstream &f, const std::string &name) noexcept
{
    const unsigned char HDF5Header[8] = {137, 72, 68, 70, 13, 10, 26, 10};
    bool isHDF5 = false;
    char header[8] = "       ";
    f.read(header, 8);
    if (f && !std::memcmp(header, HDF5Header, 8))
    {
        printf("Hierarchical Data Format (version 5) data\n");
        isHDF5 = true;
    }
    return isHDF5;
}

bool introspectAsBPFile(std::ifstream &f, const std::string &name) noexcept
{
    const int MFOOTERSIZE = 56;
    std::vector<char> buffer(MFOOTERSIZE, 0);
    f.seekg(0, f.end);
    auto flength = f.tellg();
    if (flength < MFOOTERSIZE)
    {
        return false;
    }
    f.seekg(-MFOOTERSIZE, f.end);
    f.read(buffer.data(), MFOOTERSIZE);
    if (f)
    {
        const uint8_t endianness = static_cast<uint8_t>(buffer[52]);
        if (endianness > 1)
        {
            return false;
        }
        bool IsBigEndian = (endianness == 1) ? true : false;

        const int8_t fileType = static_cast<int8_t>(buffer[54]);
        if (fileType != 0 && fileType != 2 && fileType != 3)
        {
            return false;
        }
        const int8_t BPVersion = static_cast<uint8_t>(buffer[55]);
        if (BPVersion < 1 || BPVersion > 3)
        {
            return false;
        }

        size_t position = 0;
        std::string VersionTag(buffer.data(), 28);
        position = 28;

        if (!IsBigEndian)
        {
            uint64_t PGIndexStart = helper::ReadValue<uint64_t>(buffer, position, !IsBigEndian);
            uint64_t VarsIndexStart = helper::ReadValue<uint64_t>(buffer, position, !IsBigEndian);
            uint64_t AttributesIndexStart =
                helper::ReadValue<uint64_t>(buffer, position, !IsBigEndian);
            if (PGIndexStart >= VarsIndexStart || VarsIndexStart >= AttributesIndexStart ||
                AttributesIndexStart >= static_cast<uint64_t>(flength))
            {
                return false;
            }
        }

        if (BPVersion < 3)
        {
            printf("ADIOS-BP Version %d\n", BPVersion);
        }
        else
        {
            uint8_t major = static_cast<uint8_t>(buffer[24]);
            uint8_t minor = static_cast<uint8_t>(buffer[25]);
            uint8_t patch = static_cast<uint8_t>(buffer[26]);
            if (major > '0')
            {
                // ADIOS2 writes these as ASCII characters
                major -= '0';
                minor -= '0';
                patch -= '0';
            }

            /* Cleanup ADIOS2 bug here: VersionTag is not filled with 0s */
            int pos = 10;
            while (VersionTag[pos] == '.' || (VersionTag[pos] >= '0' && VersionTag[pos] <= '9'))
            {
                ++pos;
            }
            VersionTag[pos] = '\0';
            printf("ADIOS-BP Version %d %s - ADIOS v%d.%d.%d\n", BPVersion,
                   (IsBigEndian ? "Big Endian" : "Little Endian"), major, minor, patch);
        }
    }
    return true;
}

bool introspectAsBPDir(const std::string &name) noexcept
{
    /* Must have name/md.0 */
    const std::string MDName = name + PathSeparator + "md.0";
    std::ifstream md(MDName, std::ifstream::in | std::ifstream::binary);
    if (!md)
    {
        return false;
    }
    md.close();

    /* Must have name/md.idx */
    const std::string IdxName = name + PathSeparator + "md.idx";
    std::ifstream f(IdxName, std::ifstream::in | std::ifstream::binary);
    if (!f)
    {
        return false;
    }

    const int HEADERSIZE = 64;
    std::vector<char> buffer(HEADERSIZE, 0);
    f.seekg(0, f.end);
    auto flength = f.tellg();
    if (flength >= HEADERSIZE)
    {
        f.seekg(0, f.beg);
        f.read(buffer.data(), HEADERSIZE);
    }
    f.close();

    if (flength == 0)
    {
        printf("This could be an active ADIOS BP output just opened but not "
               "written "
               "to yet\n");
        return true;
    }
    else if (flength < HEADERSIZE)
    {
        return false;
    }

    std::string tag(buffer.data(), 9);
    if (tag != "ADIOS-BP ")
    {
        return false;
    }

    char major = buffer[32];
    char minor = buffer[33];
    char patch = buffer[34];
    bool isBigEndian = static_cast<bool>(buffer[36]);
    uint8_t BPVersion = static_cast<uint8_t>(buffer[37]);
    uint8_t flatten = static_cast<uint8_t>(buffer[41]);
    bool isActive = false;
    if (BPVersion == 4)
    {
        isActive = static_cast<bool>(buffer[38]);
        printf("ADIOS-BP Version %d %s - ADIOS v%c.%c.%c %s\n", BPVersion,
               (isBigEndian ? "Big Endian" : "Little Endian"), major, minor, patch,
               (isActive ? "- active" : ""));
    }
    else if (BPVersion == 5)
    {
        uint8_t minversion = static_cast<uint8_t>(buffer[38]);
        isActive = static_cast<bool>(buffer[39]);
        printf("ADIOS-BP Version %d.%d %s - ADIOS v%c.%c.%c %s%s\n", BPVersion, minversion,
               (isBigEndian ? "Big Endian" : "Little Endian"), major, minor, patch,
               (isActive ? "- active" : ""), (flatten ? "- flatten_steps " : ""));
    }
    else
    {
        printf("ADIOS-BP Version %d %s - ADIOS v%c.%c.%c %s\n", BPVersion,
               (isBigEndian ? "Big Endian" : "Little Endian"), major, minor, patch,
               (isActive ? "- active" : ""));
    }

    return true;
}

void introspect_file(const char *filename) noexcept
{
    if (adios2sys::SystemTools::FileIsDirectory(filename))
    {
        if (!introspectAsBPDir(filename))
        {
            printf("bpls does not recognize this directory as an ADIOS "
                   "dataset\n");
        }
    }
    else
    {
        std::ifstream f(filename, std::ifstream::in | std::ifstream::binary);
        if (f)
        {
            if (!introspectAsHDF5File(f, filename))
            {
                if (!introspectAsBPFile(f, filename))
                {
                    printf("bpls does not recognize this file\n");
                }
            }
            f.close();
        }
        else
        {
            printf("File cannot be opened: %s\n", filename);
        }
    }
}

int process_unused_args(adios2sys::CommandLineArguments &arg)
{
    int nuargs;
    char **uargs;
    arg.GetUnusedArguments(&nuargs, &uargs);

    std::vector<char *> retry_args;
    retry_args.push_back(new char[4]());

    // first arg is argv[0], so skip that
    for (int i = 1; i < nuargs; i++)
    {
        if (uargs[i] != NULL && uargs[i][0] == '-')
        {
            if (uargs[i][1] == '-')
            {
                fprintf(stderr, "Unknown long option: %s\n", uargs[i]);
                arg.DeleteRemainingArguments(nuargs, &uargs);
                return 1;
            }
            else
            {
                // Maybe -abc is -a -b -c?
                size_t len = strlen(uargs[i]);
                for (size_t j = 1; j < len; ++j)
                {
                    char *opt = new char[3];
                    opt[0] = '-';
                    opt[1] = uargs[i][j];
                    opt[2] = '\0';
                    retry_args.push_back(opt);
                }
            }
        }
        else if (vfile == NULL)
        {
            vfile = mystrndup(uargs[i], 4096);
            // fprintf(stderr, "Set file argument: %s\n", vfile);
        }
        else
        {
            varmask[nmasks] = mystrndup(uargs[i], 256);
            // fprintf(stderr, "Set mask %d argument: %s\n", nmasks,
            //        varmask[nmasks]);
            nmasks++;
        }
    }
    arg.DeleteRemainingArguments(nuargs, &uargs);

    if (retry_args.size() > 1)
    {
        // Run a new parse on the -a single letter arguments
        // fprintf(stderr, "Rerun parse on %zu options\n",
        // retry_args.size());
        arg.Initialize(static_cast<int>(retry_args.size()), retry_args.data());
        arg.StoreUnusedArguments(false);
        if (!arg.Parse())
        {
            fprintf(stderr, "Parsing arguments failed\n");
            return 1;
        }
        for (size_t j = 0; j < retry_args.size(); ++j)
        {
            delete[] retry_args[j];
        }
    }
    else
    {
        delete[] retry_args[0];
    }

    return 0;
}

/** Main */
int bplsMain(int argc, char *argv[])
{
    int retval = 0;

    init_globals();

    adios2sys::CommandLineArguments arg;
    arg.Initialize(argc, argv);
    typedef adios2sys::CommandLineArguments argT;
    arg.StoreUnusedArguments(true);
    arg.AddCallback("-v", argT::NO_ARGUMENT, optioncb_verbose, nullptr, "");
    arg.AddCallback("--verbose", argT::NO_ARGUMENT, optioncb_verbose, nullptr,
                    "Print information about what bpls is doing");
    arg.AddCallback("--help", argT::NO_ARGUMENT, optioncb_help, &arg, "Help");
    arg.AddCallback("-h", argT::NO_ARGUMENT, optioncb_help, &arg, "");
    arg.AddBooleanArgument("--dump", &dump, "Dump matched variables/attributes");
    arg.AddBooleanArgument("-d", &dump, "");
    arg.AddBooleanArgument("--long", &longopt,
                           "Print values of all scalars and attributes and min/max "
                           "values of arrays");
    arg.AddBooleanArgument("-l", &longopt, "");
    arg.AddBooleanArgument("--regexp", &use_regexp,
                           "| -e Treat masks as extended regular expressions");
    arg.AddBooleanArgument("-e", &use_regexp, "");
    arg.AddArgument("--output", argT::SPACE_ARGUMENT, &outpath,
                    "| -o opt    Print to a file instead of stdout");
    arg.AddArgument("-o", argT::SPACE_ARGUMENT, &outpath, "");
    arg.AddArgument("--start", argT::SPACE_ARGUMENT, &start,
                    "| -s opt    Offset indices in each dimension (default is "
                    "0 for all dimensions).  opt<0 is handled as in python (-1 "
                    "is last)");
    arg.AddArgument("-s", argT::SPACE_ARGUMENT, &start, "");
    arg.AddArgument("--count", argT::SPACE_ARGUMENT, &count,
                    "| -c opt    Number of elements in each dimension. -1 "
                    "denotes 'until end' of dimension. default is -1 for all "
                    "dimensions");
    arg.AddArgument("-c", argT::SPACE_ARGUMENT, &count, "");
    arg.AddBooleanArgument("--noindex", &noindex, " | -y Print data without array indices");
    arg.AddBooleanArgument("-y", &noindex, "");
    arg.AddBooleanArgument("--timestep", &timestep, " | -t Print values of timestep elements");
    arg.AddBooleanArgument("--ignore_flatten", &ignore_flatten, " Don't flatten steps to one");
    arg.AddBooleanArgument("-t", &timestep, "");
    arg.AddBooleanArgument("--attrs", &listattrs, " | -a List/match attributes too");
    arg.AddBooleanArgument("-a", &listattrs, "");
    arg.AddBooleanArgument("--attrsonly", &attrsonly,
                           " | -A List/match attributes only (no variables)");
    arg.AddBooleanArgument("-A", &attrsonly, "");
    arg.AddBooleanArgument("--meshes", &listmeshes, " | -m List meshes");
    arg.AddBooleanArgument("-m", &listmeshes, "");
    arg.AddBooleanArgument("--string", &printByteAsChar,
                           " | -S Print 8bit integer arrays as strings");
    arg.AddBooleanArgument("-S", &printByteAsChar, "");
    arg.AddArgument("--columns", argT::SPACE_ARGUMENT, &ncols,
                    "| -n opt    Number of data elements per row to print");
    arg.AddArgument("-n", argT::SPACE_ARGUMENT, &ncols, "");
    arg.AddArgument("--format", argT::SPACE_ARGUMENT, &format,
                    "| -f opt    Format string to use for one data item ");
    arg.AddArgument("-f", argT::SPACE_ARGUMENT, &format, "");
    arg.AddBooleanArgument("--hidden_attrs", &hidden_attrs,
                           "  Show hidden ADIOS attributes in the file");
    arg.AddBooleanArgument("--decompose", &show_decomp,
                           "| -D Show decomposition of variables as layed out in file");
    arg.AddBooleanArgument("-D", &show_decomp, "");
    arg.AddBooleanArgument("--version", &show_version,
                           "Print version information (add -verbose for additional"
                           " information)");
    arg.AddBooleanArgument("-V", &show_version, "");
    arg.AddArgument("--error", argT::SPACE_ARGUMENT, &accuracy_def,
                    "| -X string    Specify read accuracy (error,norm,rel|abs)");
    arg.AddArgument("-X", argT::SPACE_ARGUMENT, &accuracy_def, "");
    arg.AddArgument("--transport-parameters", argT::SPACE_ARGUMENT, &transport_params,
                    "| -T string    Specify File transport parameters manually");
    arg.AddArgument("-T", argT::SPACE_ARGUMENT, &transport_params, "");
    arg.AddArgument("--engine", argT::SPACE_ARGUMENT, &engine_name,
                    "| -E string    Specify ADIOS Engine manually");
    arg.AddArgument("-E", argT::SPACE_ARGUMENT, &engine_name, "");
    arg.AddArgument("--engine-params", argT::SPACE_ARGUMENT, &engine_params,
                    "| -P string    Specify ADIOS Engine Parameters manually");
    arg.AddArgument("-P", argT::SPACE_ARGUMENT, &engine_params, "");

    if (!arg.Parse())
    {
        fprintf(stderr, "Parsing arguments failed\n");
        return 1;
    }
    if (option_help_was_called)
        return 0;

    retval = process_unused_args(arg);
    if (retval)
    {
        return retval;
    }
    if (option_help_was_called)
        return 0;

    if (show_version)
    {
        if (vfile == NULL)
        {
            print_bpls_version();
        }
        else
        {
            introspect_file(vfile);
        }
        return 0;
    }

    /* Check if we have a file defined */
    if (vfile == NULL)
    {
        fprintf(stderr, "Missing file name\n");
        return 1;
    }

    /* Process dimension specifications */
    parseDimSpec(start, istart);
    parseDimSpec(count, icount);

    // process the regular expressions
    if (use_regexp)
    {
        retval = compile_regexp_masks();
        if (retval)
            return retval;
    }

    if (noindex)
        commentchar = ';';
    else
        commentchar = ' ';

    if (hidden_attrs_flag)
        hidden_attrs = true;

    if (attrsonly)
        listattrs = true;

    retval = parseAccuracy();
    if (retval)
        return retval;

    if (verbose > 1)
        printSettings();

    retval = print_start(outpath);
    if (retval)
        return retval;

    /* Start working */
    size_t len = strlen(vfile);
    if (vfile[len - 1] == '/')
    {
        vfile[len - 1] = '\0';
    }
    retval = doList(vfile);

    print_stop();

    /* Free allocated memories */
    for (int i = 0; i < nmasks; i++)
    {
        myfree(varmask[i]);
#ifdef USE_C_REGEX
        regfree(&(varregex[i]));
#else
        varregex.clear();
#endif
    }
    myfree(vfile);

    return retval;
}

void init_globals()
{
    int i;
    // variables for arguments
    for (i = 0; i < MAX_MASKS; i++)
        varmask[i] = NULL;
    nmasks = 0;
    vfile = NULL;
    verbose = 0;
    ncols = 6; // by default when printing ascii, print "X Y", not X: Y1 Y2...
    dump = false;
    output_xml = false;
    noindex = false;
    timestep = false;
    ignore_flatten = false;
    sortnames = false;
    listattrs = false;
    listmeshes = false;
    attrsonly = false;
    longopt = false;
    // timefrom             = 1;
    // timeto               = -1;
    use_regexp = false;
    plot = false;
    hidden_attrs = false;
    hidden_attrs_flag = 0;
    printByteAsChar = false;
    show_decomp = false;
    show_version = false;
    for (i = 0; i < MAX_DIMS; i++)
    {
        istart[i] = 0LL;
        icount[i] = -1LL; // read full var by default
    }
    ndimsspecified = 0;
}

#define PRINT_DIMS_INT(str, v, n, loopvar)                                                         \
    printf("%s = { ", str);                                                                        \
    for (loopvar = 0; loopvar < n; loopvar++)                                                      \
        printf("%d ", v[loopvar]);                                                                 \
    printf("}")

#define PRINT_DIMS_UINT64(str, v, n, loopvar)                                                      \
    printf("%s = { ", str);                                                                        \
    for (loopvar = 0; loopvar < n; loopvar++)                                                      \
        printf("%" PRIu64 " ", v[loopvar]);                                                        \
    printf("}")

#define PRINT_DIMS_INT64(str, v, n, loopvar)                                                       \
    printf("%s = { ", str);                                                                        \
    for (loopvar = 0; loopvar < n; loopvar++)                                                      \
        printf("%" PRId64 " ", v[loopvar]);                                                        \
    printf("}")

#define PRINT_DIMS_SIZET(str, v, n, loopvar)                                                       \
    printf("%s = { ", str);                                                                        \
    for (loopvar = 0; loopvar < n; loopvar++)                                                      \
        printf("%zu ", v[loopvar]);                                                                \
    printf("}")

void printSettings(void)
{
    int i;
    printf("Settings :\n");
    printf("  masks  : %d ", nmasks);
    for (i = 0; i < nmasks; i++)
        printf("%s ", varmask[i]);
    printf("\n");
    printf("  file   : %s\n", vfile);
    printf("  output : %s\n", (outpath.empty() ? "stdout" : outpath.c_str()));

    if (start.size())
    {
        PRINT_DIMS_INT64("  start", istart, ndimsspecified, i);
        printf("\n");
    }
    if (count.size())
    {
        PRINT_DIMS_INT64("  count", icount, ndimsspecified, i);
        printf("\n");
    }

    if (longopt)
        printf("      -l : show scalar values and min/max/avg of arrays\n");
    if (sortnames)
        printf("      -r : sort names before listing\n");
    if (attrsonly)
        printf("      -A : list attributes only\n");
    else if (listattrs)
        printf("      -a : list attributes too\n");
    if (listmeshes)
        printf("      -m : list meshes too\n");
    if (dump)
        printf("      -d : dump matching variables and attributes\n");
    if (use_regexp)
        printf("      -e : handle masks as regular expressions\n");
    if (format.size())
        printf("      -f : dump using printf format \"%s\"\n", format.c_str());
    if (output_xml)
        printf("      -x : output data in XML format\n");
    if (show_decomp)
        printf("      -D : show decomposition of variables in the file\n");
    if (show_version)
        printf("      -V : show binary version info of file\n");
    if (timestep)
        printf("      -t : read step-by-step\n");
    if (ignore_flatten)
        printf("      --ignore_flatten : ignore FlattenSteps writer specification\n");

    if (hidden_attrs)
    {
        printf("         : show hidden attributes in the file\n");
    }
}

void bpexit(int code, core::Engine *fp)
{
    if (fp != nullptr)
        fp->Close();
    exit(code);
}

void print_file_size(uint64_t size)
{
    static const char *sm[] = {"bytes", "KB", "MB", "GB", "TB", "PB", "EB"};
    uint64_t s = size, r = 0;
    int idx = 0;
    while (s / 1024 > 0)
    {
        r = s % 1024;
        s = s / 1024;
        idx++;
    }
    if (r > 511)
        s++;
    printf("  file size:     %" PRIu64 " %s\n", s, sm[idx]);
}

static inline int ndigits(size_t n)
{
    static char digitstr[32];
    return snprintf(digitstr, 32, "%zu", n);
}

template <class T>
int printAttributeValue(core::Engine *fp, core::IO *io, core::Attribute<T> *attribute)
{
    DataType adiosvartype = attribute->m_Type;
    if (attribute->m_IsSingleValue)
    {
        print_data((void *)&attribute->m_DataSingleValue, 0, adiosvartype, true);
    }
    else
    {
        fprintf(outf, "{");
        size_t nelems = attribute->m_DataArray.size();
        for (size_t j = 0; j < nelems; j++)
        {
            print_data((void *)&attribute->m_DataArray[j], 0, adiosvartype, true);
            if (j < nelems - 1)
            {
                fprintf(outf, ", ");
            }
        }
        fprintf(outf, "}");
    }
    return 0;
}

template <>
int printAttributeValue(core::Engine *fp, core::IO *io, core::Attribute<std::string> *attribute)
{
    DataType adiosvartype = attribute->m_Type;
    bool xmlprint = helper::EndsWith(attribute->m_Name, "xml", false);
    bool printDataAnyway = true;

    if (attribute->m_IsSingleValue)
    {
        if (xmlprint)
        {
            printDataAnyway = print_data_xml(attribute->m_DataSingleValue.data(),
                                             attribute->m_DataSingleValue.length());
        }
        if (printDataAnyway)
        {
            print_data((void *)&attribute->m_DataSingleValue, 0, adiosvartype, true);
        }
    }
    else
    {
        fprintf(outf, "{");
        size_t nelems = attribute->m_DataArray.size();
        for (size_t j = 0; j < nelems; j++)
        {
            if (xmlprint)
            {
                printDataAnyway = print_data_xml(attribute->m_DataArray[j].data(),
                                                 attribute->m_DataArray[j].length());
            }
            if (printDataAnyway)
            {
                print_data((void *)&attribute->m_DataArray[j], 0, adiosvartype, true);
            }
            if (j < nelems - 1)
            {
                fprintf(outf, ", ");
            }
        }
        fprintf(outf, "}");
    }
    return 0;
}

int nEntriesMatched = 0;

int doList_vars(core::Engine *fp, core::IO *io)
{

    const core::VarMap &variables = io->GetVariables();
    const core::AttrMap &attributes = io->GetAttributes();

    // make a sorted list of all variables and attributes
    EntryMap entries;
    if (!attrsonly)
    {
        for (const auto &vpair : variables)
        {
            Entry e(vpair.second->m_Type, vpair.second.get());
            bool valid = true;
            if (timestep && !filestream)
            {
                valid = e.var->IsValidStep(fp->CurrentStep() + 1);
                // fprintf(stdout, "Entry: ptr = %p valid = %d\n", e.var,
                // valid);
            }
            if (valid)
            {
                entries.emplace(vpair.first, e);
            }
        }
    }
    if (listattrs)
    {
        for (const auto &apair : attributes)
        {
            Entry e(apair.second->m_Type, apair.second.get());
            entries.emplace(apair.first, e);
        }
    }

    // size_t nNames = entries.size();

    // calculate max length of variable names and type names in the first
    // round
    int maxlen = 4; // need int for printf formatting
    int maxtypelen = 7;
    for (const auto &entrypair : entries)
    {
        int len = static_cast<int>(entrypair.first.size());
        if (len > maxlen)
            maxlen = len;
        len = static_cast<int>(ToString(entrypair.second.typeName).size());
        if (len > maxtypelen)
            maxtypelen = len;
    }

    /* VARIABLES */
    for (const auto &entrypair : entries)
    {
        int retval = 0;
        const std::string &name = entrypair.first;
        const Entry &entry = entrypair.second;
        bool matches = matchesAMask(name.c_str());
        if (matches)
        {
            nEntriesMatched++;

            // print definition of variable
            fprintf(outf, "%c %-*s  %-*s", commentchar, maxtypelen,
                    ToString(entry.typeName).c_str(), maxlen, name.c_str());
            if (!entry.isVar)
            {
                // list (and print) attribute
                if (longopt || dump)
                {
                    fprintf(outf, "  attr   = ");
                    if (entry.typeName == DataType::Struct)
                    {
                        // not supported
                    }
#define declare_template_instantiation(T)                                                          \
    else if (entry.typeName == helper::GetDataType<T>())                                           \
    {                                                                                              \
        core::Attribute<T> *a = static_cast<core::Attribute<T> *>(entry.attr);                     \
        retval = printAttributeValue(fp, io, a);                                                   \
    }
                    ADIOS2_FOREACH_ATTRIBUTE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation
                    fprintf(outf, "\n");
                    matches = false; // already printed
                }
                else
                {
                    fprintf(outf, "  attr\n");
                }
            }
            else
            {
                if (entry.typeName == DataType::Struct)
                {
                    // not supported
                }
#define declare_template_instantiation(T)                                                          \
    else if (entry.typeName == helper::GetDataType<T>())                                           \
    {                                                                                              \
        core::Variable<T> *v = static_cast<core::Variable<T> *>(entry.var);                        \
        retval = printVariableInfo(fp, io, v);                                                     \
    }
                ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation
            }
        }

        if (retval && retval != 10) // do not return after unsupported type
            return retval;
    }

    entries.clear();
    return 0;
}

template <class T>
int printVariableInfo(core::Engine *fp, core::IO *io, core::Variable<T> *variable)
{
    size_t nsteps = variable->GetAvailableStepsCount();
    if (timestep)
    {
        nsteps = 1;
    }
    DataType adiosvartype = variable->m_Type;
    int retval = 0;

    bool isGlobalValue = (nsteps == 1);
    isGlobalValue &= variable->m_SingleValue;
    isGlobalValue &= (variable->m_ShapeID != ShapeID::GlobalArray);

    if (!isGlobalValue)
    {
        fprintf(outf, "  ");

        if (nsteps > 1)
            fprintf(outf, "%zu*", nsteps);

        if (variable->m_ShapeID == ShapeID::GlobalArray)
        {
            Dims d = get_global_array_signature(fp, io, variable);
            fprintf(outf, "{%s", d[0] > 0 ? std::to_string(d[0]).c_str() : "__");
            for (size_t j = 1; j < variable->m_Shape.size(); j++)
            {
                fprintf(outf, ", %s", d[j] > 0 ? std::to_string(d[j]).c_str() : "__");
            }
            fprintf(outf, "}");
        }
        else if (variable->m_ShapeID == ShapeID::LocalArray)
        {
            std::pair<size_t, Dims> signo = get_local_array_signature(fp, io, variable);
            Dims &d = signo.second;
            fprintf(outf, "[%s]*", signo.first > 0 ? std::to_string(signo.first).c_str() : "__");
            fprintf(outf, "{%s", d[0] > 0 ? std::to_string(d[0]).c_str() : "__");
            for (size_t j = 1; j < variable->m_Count.size(); j++)
            {
                fprintf(outf, ", %s", d[j] > 0 ? std::to_string(d[j]).c_str() : "__");
            }
            fprintf(outf, "}");
        }
        else
        {
            fprintf(outf, "scalar");
        }
#if 0
        if (longopt || plot)
        {
            adios_inq_var_stat(fp, vi, timestep && timed, show_decomp);
        }

        if (plot && vi->statistics && vi->statistics->histogram)
        {
            print_data_hist(vi, &names[n][1]);
        }
#endif
        if (longopt /* TODO: && variable->has_statistics */)
        {
            if (timestep == false)
            {
                MinMaxStruct MinMax;
                try
                {
                    if (fp->VariableMinMax(*variable, DefaultSizeT, MinMax))
                    {
                        fprintf(outf, " = ");
                        print_data(&MinMax.MinUnion, 0, adiosvartype, false);
                        fprintf(outf, " / ");
                        print_data(&MinMax.MaxUnion, 0, adiosvartype, false);
                    }
                    else
                    {
                        fprintf(outf, " = ");
                        print_data(&variable->m_Min, 0, adiosvartype, false);
                        fprintf(outf, " / ");
                        print_data(&variable->m_Max, 0, adiosvartype, false);
                    }
                    // fprintf(outf," {MIN / MAX} ");
                }
                catch (std::logic_error &)
                {
                }
            }
#if 0
            else
            {
                int time_start = 0, time_end = vi->nsteps;

                if (start != NULL)
                {
                    if (istart[0] >= 0)
                        time_start = istart[0];
                    else
                        time_start = vi->nsteps - 1 + (int)istart[0];
                }

                if (count != NULL)
                {
                    if (icount[0] > 0)
                        time_end = time_start + (int)icount[0];
                    else
                        time_end = vi->nsteps + (int)icount[0] + 1;
                }

                if (time_start < 0 || time_start >= vi->nsteps)
                {
                    fprintf(stderr, "Error when reading variable %s. "
                                    "errno=%d : Variable (id=%d) has "
                                    "no data at %d time step\n",
                            names[n], 15, vi->varid, time_start);
                    bpexit(15, fp);
                }

                if (time_end < 0 || time_end > vi->nsteps)
                {
                    fprintf(stderr, "Error when reading variable %s. "
                                    "errno=%d : Variable (id=%d) has "
                                    "no data at %d time step\n",
                            names[n], 15, vi->varid, time_end);
                    bpexit(16, fp);
                }

                static char *indent_char = " ";
                int indent_len = 11;

                /* Start - Print the headers of statistics first */
                fprintf(outf, "\n%-*s", indent_len + 7, indent_char);
                fprintf(outf, "%10s  ", "MIN");
                fprintf(outf, "%10s  ", "MAX");
                fprintf(outf, "%10s  ", "AVG");
                fprintf(outf, "%10s  ", "STD DEV");

                /* End - Print the headers of statistics first */

                void *min, *max, *avg, *std_dev;
                DataType vt = vartype;
                struct ADIOS_STAT_STEP *s = vi->statistics->steps;
                if (vi->type == DataType::FloatComplex ||
                    vi->type == DataType::DoubleComplex)
                    vt = DataType::Double;
                fprintf(outf, "\n%-*sglobal:", indent_len, indent_char);
                print_data_characteristics(
                    vi->statistics->min, vi->statistics->max,
                    vi->statistics->avg, vi->statistics->std_dev, vt, false);

                for (i = time_start; i < time_end; i++)
                {
                    min = max = avg = std_dev = 0;
                    if (s->maxs && s->maxs[i])
                        max = s->maxs[i];
                    if (s->mins && s->mins[i])
                        min = s->mins[i];
                    if (s->avgs && s->avgs[i])
                        avg = s->avgs[i];
                    if (s->std_devs && s->std_devs[i])
                        std_dev = s->std_devs[i];

                    // Align the output, previous lines has atleast
                    // (maxlen + strlen(names[n])) characters
                    // Better way to printf N spaces?
                    fprintf(outf, "\n%-*st%-5d:", indent_len, indent_char, i);
                    print_data_characteristics(min, max, avg, std_dev, vt,
                                               false);
                }
                fprintf(outf, "\n");
            }
#endif
        } // longopt && vi->statistics
        fprintf(outf, "\n");

        if (show_decomp)
        {
            if (timestep)
            {
                print_decomp_singlestep(fp, io, variable);
            }
            else
            {
                print_decomp(fp, io, variable);
            }
        }
    }
    else
    {
        // single GlobalValue without timesteps
        fprintf(outf, "  scalar");
        if (longopt && !timestep)
        {
            fprintf(outf, " = ");
            T value;
            fp->Get(*variable, value, adios2::Mode::Sync);
            print_data(&value, 0, adiosvartype, false);
        }
        fprintf(outf, "\n");

        if (show_decomp)
        {
            if (timestep)
            {
                print_decomp_singlestep(fp, io, variable);
            }
            else
            {
                print_decomp(fp, io, variable);
            }
        }
    }

    if (dump && !show_decomp)
    {
        variable->SetAccuracy(accuracy);
        // print variable content
        if (variable->m_ShapeID == ShapeID::LocalArray)
        {
            if (timestep)
            {
                print_decomp_singlestep(fp, io, variable);
            }
            else
            {
                print_decomp(fp, io, variable);
            }
        }
        else
        {
            retval = readVar(fp, io, variable);
        }
        fprintf(outf, "\n");
    }
    return retval;
}

#define PRINT_ARRAY(str, ndim, dims, loopvar, format)                                              \
    fprintf(outf, "%s", str);                                                                      \
    if (ndim > 0)                                                                                  \
    {                                                                                              \
        fprintf(outf, "{%" #format, dims[0]);                                                      \
        for (loopvar = 1; loopvar < ndim; loopvar++)                                               \
        {                                                                                          \
            fprintf(outf, ", %" #format, dims[loopvar]);                                           \
        }                                                                                          \
        fprintf(outf, "}\n");                                                                      \
    }                                                                                              \
    else                                                                                           \
    {                                                                                              \
        fprintf(outf, "empty\n");                                                                  \
    }

#define PRINT_ARRAY64(str, ndim, dims, loopvar)                                                    \
    fprintf(outf, "%s", str);                                                                      \
    if (ndim > 0)                                                                                  \
    {                                                                                              \
        fprintf(outf, "{%" PRIu64, dims[0]);                                                       \
        for (loopvar = 1; loopvar < ndim; loopvar++)                                               \
        {                                                                                          \
            fprintf(outf, ", %" PRIu64, dims[loopvar]);                                            \
        }                                                                                          \
        fprintf(outf, "}\n");                                                                      \
    }                                                                                              \
    else                                                                                           \
    {                                                                                              \
        fprintf(outf, "empty\n");                                                                  \
    }

void printMeshes(core::Engine *fp)
{
    fprintf(outf, "Mesh info: is not implemented in adios 2.x at the moment\n");
    return;
    /*
    int meshid, i, j; // loop vars
    int mpi_comm_dummy = 0;
    if (fp->nmeshes == 0)
    {
        fprintf(outf, "Mesh info: There are no meshes defined in this
    file\n");
        return;
    }
    fprintf(outf, "Mesh info: \n");
    for (meshid = 0; meshid < fp->nmeshes; meshid++)
    {
        fprintf(outf, "  %s\n", fp->mesh_namelist[meshid]);
        ADIOS_MESH *mi = adios_inq_mesh_byid(fp, meshid);
        if (mi)
        {
            if (meshid != mi->id)
                fprintf(
                    outf,
                    "  bpls warning: meshid (=%d) != inquired mesh id
    (%d)\n",
                    meshid, mi->id);
            if (strcmp(fp->mesh_namelist[meshid], mi->name))
                fprintf(outf, "  bpls warning: mesh name in list (=\"%s\")
    != "
                              "inquired mesh name (\"%s\")\n",
                        fp->mesh_namelist[meshid], mi->name);
            if (mi->file_name)
            {
                IO &io = adios.DeclareIO("mesh");
                Engine &meshfp = io->Open(mi->file_name, Mode::Read);
                adios_complete_meshinfo(fp, meshfp, mi);
                meshfp.Close();
            }
            fprintf(outf, "    type:         ");

            switch (mi->type)
            {
            case ADIOS_MESH_UNIFORM:
                fprintf(outf, "uniform\n");
                PRINT_ARRAY64("    dimensions:   ",
    mi->uniform->num_dimensions,
                              mi->uniform->dimensions, j)
                if (mi->uniform->origins)
                {
                    PRINT_ARRAY("    origins:      ",
                                mi->uniform->num_dimensions,
                                mi->uniform->origins, j, g)
                }
                if (mi->uniform->spacings)
                {
                    PRINT_ARRAY("    spacings:     ",
                                mi->uniform->num_dimensions,
                                mi->uniform->spacings, j, g)
                }
                if (mi->uniform->maximums)
                {
                    PRINT_ARRAY("    maximums:     ",
                                mi->uniform->num_dimensions,
                                mi->uniform->maximums, j, g)
                }
                break;

            case ADIOS_MESH_RECTILINEAR:
                fprintf(outf, "rectilinear\n");
                PRINT_ARRAY64("    dimensions:   ",
                              mi->rectilinear->num_dimensions,
                              mi->rectilinear->dimensions, j)
                if (mi->rectilinear->use_single_var)
                {
                    fprintf(outf, "    coordinates:  single-var: \"%s\"\n",
                            mi->rectilinear->coordinates[0]);
                }
                else
                {
                    fprintf(outf, "    coordinates:  multi-var: \"%s\"",
                            mi->rectilinear->coordinates[0]);
                    for (i = 1; i < mi->rectilinear->num_dimensions; i++)
                    {
                        fprintf(outf, ", \"%s\"",
                                mi->rectilinear->coordinates[i]);
                    }
                    fprintf(outf, "\n");
                }
                break;

            case ADIOS_MESH_STRUCTURED:
                fprintf(outf, "structured\n");
                PRINT_ARRAY64("    dimensions:   ",
                              mi->structured->num_dimensions,
                              mi->structured->dimensions, j);
                if (mi->structured->use_single_var)
                {
                    fprintf(outf, "    points:       single-var: \"%s\"\n",
                            mi->structured->points[0]);
                }
                else
                {
                    fprintf(outf, "    points:       multi-var: \"%s\"",
                            mi->structured->points[0]);
                    for (i = 1; i < mi->structured->num_dimensions; i++)
                    {
                        fprintf(outf, ", \"%s\"",
    mi->structured->points[i]);
                    }
                    fprintf(outf, "\n");
                }
                fprintf(outf, "    nspaces:      %d\n",
                        mi->structured->nspaces);
                break;

            case ADIOS_MESH_UNSTRUCTURED:
                fprintf(outf, "unstructured\n");
                if (mi->unstructured->nvar_points <= 1)
                {
                    fprintf(outf, "    npoints:      %" PRIu64 "\n",
                            mi->unstructured->npoints);
                    fprintf(outf, "    points:       single-var: \"%s\"\n",
                            mi->unstructured->points[0]);
                }
                else
                {
                    fprintf(outf, "    points:       multi-var: \"%s\"",
                            mi->unstructured->points[0]);
                    for (i = 1; i < mi->unstructured->nvar_points; i++)
                    {
                        fprintf(outf, ", \"%s\"",
    mi->unstructured->points[i]);
                    }
                    fprintf(outf, "\n");
                }
                fprintf(outf, "    ncsets:       %d\n",
                        mi->unstructured->ncsets);
                for (i = 0; i < mi->unstructured->ncsets; i++)
                {
                    fprintf(outf, "    cell set %d:\n", i);
                    fprintf(outf, "      cell type:  %d\n",
                            mi->unstructured->ctypes[i]);
                    fprintf(outf, "      ncells:     %" PRIu64 "\n",
                            mi->unstructured->ccounts[i]);
                    fprintf(outf, "      cells var:  \"%s\"\n",
                            mi->unstructured->cdata[i]);
                }
                fprintf(outf, "    nspaces:      %d\n",
                        mi->unstructured->nspaces);
                break;

            default:
                fprintf(outf, "undefined\n");
            }
            fprintf(outf, "    time varying: %s\n",
                    (mi->time_varying ? "yes" : "no"));
            adios_free_meshinfo(mi);
        }
    }
    fprintf(outf, "\n");
    */
}

std::vector<std::string> getEnginesList(const std::string path)
{
    std::vector<std::string> list;
#ifdef ADIOS2_HAVE_HDF5
    size_t slen = path.length();
    if (slen >= 3 && path.compare(slen - 3, 3, ".h5") == 0)
    {
        list.push_back("HDF5");
        if (timestep)
        {
            list.push_back("FileStream");
            list.push_back("BP3");
        }
        else
        {
            list.push_back("BPFile");
        }
    }
    else
    {
        if (timestep)
        {
            list.push_back("FileStream");
            list.push_back("BP3");
        }
        else
        {
            list.push_back("BPFile");
        }
        list.push_back("HDF5");
    }
#else
    if (timestep)
    {
        list.push_back("FileStream");
        list.push_back("BP3");
    }
    else
    {
        list.push_back("BPFile");
    }
#endif
    return list;
}

int doList(std::string path)
{
    char init_params[128];
    int adios_verbose = 2;

    if (verbose > 1)
        printf("\nADIOS Open: read header info from %s\n", path.c_str());

    // initialize BP reader
    if (verbose > 1)
        adios_verbose = 3; // print info lines
    if (verbose > 2)
        adios_verbose = 4; // print debug lines
    snprintf(init_params, sizeof(init_params), "verbose=%d", adios_verbose);
    if (hidden_attrs)
        strcat(init_params, ";show_hidden_attrs");

    core::ADIOS adios("C++");
    const adios2::UserOptions userOptions = adios.GetUserOptions();

    std::string tpl = helper::LowerCase(transport_params);
    bool remoteFile =
        (tpl.find("awssdk") != std::string::npos) || (tpl.find("daos") != std::string::npos);
    if (remoteFile)
    {
        if (engine_name.empty())
        {
            fprintf(stderr, "\nError: For remote file access you must specify the engine "
                            "explicitly with -E parameter, e.g. -E bp5 or -E bp4 or -E "
                            "bp3.\nVirtual engines like BPFile or FileStream do not know "
                            "which engine to use.\n");
            return 8;
        }
    }
    else
    {
        bool exists = adios2sys::SystemTools::FileExists(path);
        if (!exists && !userOptions.campaign.campaignstorepath.empty() && path[0] != PathSeparator)
        {
            std::string path2 = userOptions.campaign.campaignstorepath + PathSeparator + path;
            exists = adios2sys::SystemTools::FileExists(path2);
            if (exists)
            {
                ; // path = path2.c_str();
            }
            else
            {
                std::string path3 =
                    userOptions.campaign.campaignstorepath + PathSeparator + path + ".aca";
                exists = adios2sys::SystemTools::FileExists(path3);
                if (exists)
                {
                    path += ".aca"; // path = path3.c_str();
                }
            }
        }

        if (!exists)
        {
            fprintf(stderr, "\nError: input path %s does not exist\n", path.c_str());
            return 4;
        }
    }

    core::IO &io = adios.DeclareIO("bpls");
    if (timestep)
    {
        // BP4 can process metadata in chuncks to conserve memory
        io.SetParameter("StreamReader", "true");
    }
    core::Engine *fp = nullptr;

    if (!transport_params.empty())
    {
        auto p = helper::BuildParametersMap(transport_params, '=', ',');
        io.AddTransport("File", p);
    }

    std::vector<std::string> engineList;
    if (engine_name.empty())
    {
        engineList = getEnginesList(path);
    }
    else
    {
        engineList.push_back(engine_name);
    }

    if (!engine_params.empty())
    {
        auto p = helper::BuildParametersMap(engine_params, '=', ',');
        io.SetParameters(p);
    }

    if (ignore_flatten)
    {
        io.SetParameters("IgnoreFlattenSteps=on");
    }

    for (auto &engineName : engineList)
    {
        if (verbose > 2)
            printf("Try %s engine to open the file...\n", engineName.c_str());
        io.SetEngine(engineName);
        try
        {
            if (timestep)
            {
                fp = &io.Open(path, Mode::Read);
            }
            else
            {
                fp = &io.Open(path, Mode::ReadRandomAccess);
            }
            if (engineName == "FileStream")
            {
                filestream = true;
            }
        }
        catch (std::exception &e)
        {
            printf("Failed to open with %s engine: %s\n", engineName.c_str(), e.what());
        }
        if (fp != nullptr)
            break;
    }

    if (fp == nullptr)
    {
        fprintf(stderr, "\nError: Could not open this file with any ADIOS2 "
                        "file reading engines\n");
        return 4;
    }

    {
        //, variables, timesteps, and attributes
        // all parameters are integers,
        // besides the last parameter, which is an array of strings for
        // holding
        // the
        // list of group names
        // ntsteps = fp->tidx_stop - fp->tidx_start + 1;
        if (verbose)
        {
            printf("File info:\n");
            if (!timestep)
            {
                printf("  of variables:  %zu\n", io.GetVariables().size());
                printf("  of attributes: %zu\n", io.GetAttributes().size());
            }
            // printf("  of meshes:     %d\n", fp->nmeshes);
            // print_file_size(fp->file_size);
            // printf("  bp version:    %d\n", fp->version);
            // printf("  endianness:    %s\n",
            //       (fp->endianness ? "Big Endian" : "Little Endian"));
            if (longopt)
                printf("  statistics:    Min / Max \n");
            printf("\n");
        }

        // Print out the meshes in the file
        if (listmeshes)
        {
            printMeshes(fp);
        }

        if (timestep)
        {
            while (true)
            {
                adios2::StepStatus status = fp->BeginStep(adios2::StepMode::Read);
                if (status == adios2::StepStatus::EndOfStream)
                {
                    break;
                }
                else if (status == adios2::StepStatus::NotReady)
                {
                    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
                    continue;
                }
                else if (status == adios2::StepStatus::OtherError)
                {
                    fprintf(stderr, "\nError: Cannot read more steps due to errors\n");
                    break;
                }
                fprintf(stdout, "Step %zu:\n", fp->CurrentStep());
                doList_vars(fp, &io);
                fp->EndStep();
            }
        }
        else
        {
            doList_vars(fp, &io);
        }

        if (nmasks > 0 && nEntriesMatched == 0)
        {
            fprintf(stderr, "\nError: None of the variables/attributes matched any "
                            "name/regexp you provided\n");
            return 4;
        }
        fp->Close();
    }
    return 0;
}

/*
int print_data_hist(ADIOS_VARINFO *vi, char *varname)
{
    char hist_file[256], gnuplot_file[256];
    int i;
    char xtics[512], str[512];
    FILE *out_hist, *out_plot;
    struct ADIOS_HIST *h = vi->statistics->histogram;

    memcpy(hist_file, varname, strlen(varname) + 1);
    strcat(hist_file, ".hist");

    if ((out_hist = fopen(hist_file, "w")) == NULL)
    {
        fprintf(stderr, "Error at opening for writing file %s: %s\n",
hist_file,
                strerror(errno));
        return 30;
    }

    memcpy(gnuplot_file, varname, strlen(varname) + 1);
    strcat(gnuplot_file, ".gpl");

    if ((out_plot = fopen(gnuplot_file, "w")) == NULL)
    {
        fprintf(stderr, "Error at opening for writing file %s: %s\n",
                gnuplot_file, strerror(errno));
        return 30;
    }

    xtics[0] = '\0';
    strcat(xtics, "set xtics offset start axis (");
    for (i = 0; i <= h->num_breaks; i++)
    {
        if (i == 0)
        {
            fprintf(out_hist, "-Inf %.2lf %u\n", h->breaks[i],
                    h->gfrequencies[i]);
            sprintf(str, "\"-Inf\" pos(%d)", i);
        }
        else if (i < h->num_breaks)
        {
            fprintf(out_hist, "%.2lf %.2lf %u\n", h->breaks[i - 1],
                    h->breaks[i], h->gfrequencies[i]);
            sprintf(str, ", \"%.2lf\" pos(%d)", h->breaks[i - 1], i);
        }
        else
        {
            fprintf(out_hist, "%.2lf Inf %u\n", h->breaks[i],
                    h->gfrequencies[i]);
            sprintf(str, ", \"Inf\" pos(%d)", i);
        }
        strcat(xtics, str);
    }
    strcat(xtics, ")\n");

    fprintf(out_plot, "start = -0.5\npos(x) = start + x * 1\nset boxwidth "
                      "1\nset style fill solid border 5#5lt6#6\n");
    fputs(xtics, out_plot);
    fprintf(out_plot, "plot '%s' using 3 smooth frequency w boxes\n",
            hist_file);
    fprintf(out_plot, "pause -1 'Press Enter to quit'\n");
    return 0;
}
*/

int cmpstringp(const void *p1, const void *p2)
{
    /* The actual arguments to this function are "pointers to
       pointers to char", but strcmp() arguments are "pointers
       to char", hence the following cast plus dereference */
    return strcmp(*(char *const *)p1, *(char *const *)p2);
}
/** Merge listV with listA if listattrs=true, return listA if
  attrsonly=true,
  otherwise just return listV.
  If sortnames=true, quicksort the result list.
 */
void mergeLists(int nV, char **listV, int nA, char **listA, char **mlist, bool *isVar)
{
    int v, a, idx;
    if (sortnames && listattrs && !attrsonly)
    {
        // merge sort the two lists
        v = 0;
        a = 0;
        while (v < nV || a < nA)
        {
            if (a < nA && (v >= nV || strcmp(listV[v], listA[a]) > 0))
            {
                // fully consumed var list or
                // next item in attr list is less than next item in var list
                mlist[v + a] = listA[a];
                isVar[v + a] = false;
                a++;
            }
            else
            {
                mlist[v + a] = listV[v];
                isVar[v + a] = true;
                v++;
            }
        }
    }
    else
    {
        // first add vars then attrs (if ask ed)
        idx = 0;
        if (!attrsonly)
        {
            for (v = 0; v < nV; v++)
            {
                mlist[idx] = listV[v];
                isVar[idx] = true;
                idx++;
            }
        }
        if (listattrs)
        {
            for (a = 0; a < nA; a++)
            {
                mlist[idx] = listA[a];
                isVar[idx] = false;
                idx++;
            }
        }
    }
}

int getTypeInfo(DataType adiosvartype, int *elemsize)
{
    switch (adiosvartype)
    {
    case DataType::UInt8:
        *elemsize = 1;
        break;
    case DataType::Int8:
        *elemsize = 1;
        break;
    case DataType::String:
        *elemsize = 1;
        break;

    case DataType::UInt16:
        *elemsize = 2;
        break;
    case DataType::Int16:
        *elemsize = 2;
        break;

    case DataType::UInt32:
        *elemsize = 4;
        break;
    case DataType::Int32:
        *elemsize = 4;
        break;

    case DataType::UInt64:
        *elemsize = 8;
        break;
    case DataType::Int64:
        *elemsize = 8;
        break;

    case DataType::Float:
        *elemsize = 4;
        break;

    case DataType::Double:
        *elemsize = 8;
        break;

    case DataType::FloatComplex:
        *elemsize = 8;
        break;

    case DataType::DoubleComplex:
        *elemsize = 16;
        break;

    case DataType::LongDouble: // do not know how to print
    //*elemsize = 16;
    default:
        return 1;
    }
    return 0;
}

/** Read data of a variable and print
 * Return: 0: ok, != 0 on error
 */
template <class T>
int readVar(core::Engine *fp, core::IO *io, core::Variable<T> *variable)
{
    int i, j;
    uint64_t start_t[MAX_DIMS],
        count_t[MAX_DIMS]; // processed <0 values in start/count
    uint64_t s[MAX_DIMS],
        c[MAX_DIMS]; // for block reading of smaller chunks
    int tdims;       // number of dimensions including time
    int tidx;        // 0 or 1 to account for time dimension
    uint64_t nelems; // number of elements to read
    // size_t elemsize;                   // size in bytes of one element
    uint64_t stepStart, stepCount;
    std::vector<T> dataV;
    uint64_t sum;             // working var to sum up things
    uint64_t maxreadn;        // max number of elements to read once up to a limit
                              // (10MB of data)
    uint64_t actualreadn;     // our decision how much to read at once
    uint64_t readn[MAX_DIMS]; // how big chunk to read in in each dimension?
    int ndigits_dims[32];     // # of digits (to print) of each dimension

    const size_t elemsize = variable->m_ElementSize;
    const int nsteps = (timestep ? 1 : static_cast<int>(variable->GetAvailableStepsCount()));
    const int ndim = static_cast<int>(variable->m_Shape.size());
    // create the counter arrays with the appropriate lengths
    // transfer start and count arrays to format dependent arrays

    nelems = 1;
    tidx = 0;
    stepStart = 0;
    stepCount = 1;

    if (!timestep && nsteps > 1)
    {
        // user selection must start with step selection
        // calculate the starting step and number of steps requested
        if (istart[0] < 0) // negative index means last-|index|
            stepStart = nsteps + istart[0];
        else
            stepStart = istart[0];
        if (icount[0] < 0) // negative index means last-|index|+1-start
            stepCount = nsteps + icount[0] + 1 - stepStart;
        else
            stepCount = icount[0];

        if (verbose > 2)
            printf("    j=0, stepStart=%" PRIu64 " stepCount=%" PRIu64 "\n", stepStart, stepCount);

        if (stepStart + stepCount > static_cast<uint64_t>(nsteps))
        {
            printf("ERROR: The sum of start step (%" PRIu64 ") and step count (%" PRIu64
                   ") is larger "
                   "than the number of steps available (%d)\n",
                   stepStart, stepCount, nsteps);
            return -1;
        }

        start_t[0] = stepStart;
        count_t[0] = stepCount;
        nelems *= stepCount;
        if (verbose > 1)
            printf("    s[0]=%" PRIu64 ", c[0]=%" PRIu64 ", n=%" PRIu64 "\n", start_t[0],
                   count_t[0], nelems);
        tidx = 1;
    }

    tdims = ndim + tidx;

    // Absolute step is needed to access Shape of variable in a step
    // and also for StepSelection
    size_t absstep = relative_to_absolute_step(fp, variable, stepStart);

    // Get the shape of the variable for the starting step
    Dims shape;
    if (timestep)
    {
        shape = variable->Shape();
    }
    else
    {
        shape = variable->Shape(absstep);
    }
    if (verbose > 2)
    {
        printf("    starting step=%" PRIu64 " absolute step=%zu"
               ", dims={",
               stepStart, absstep);
        for (auto dim : shape)
        {
            printf("%zu", dim);
        }
        printf("}\n");
    }

    if (tidx)
    {
        // If shape is changing we still can support printing a single step
        auto dimsign = get_global_array_signature(fp, io, variable);
        bool changingShape = false;
        for (auto d : dimsign)
        {
            if (d == 0)
            {
                changingShape = true;
                break;
            }
        }
        if (changingShape && stepCount > 1)
        {
            printf("ERROR: This variable has a changing shape over time, "
                   "so bpls cannot dump it as a global array. \n");
            printf("You can dump a single step with\n"
                   "    bpls -d %s -s \"T",
                   variable->m_Name.c_str());
            for (j = 0; j < ndim; j++)
            {
                printf(",0");
            }
            printf("\" -c \"1");
            for (j = 0; j < ndim; j++)
            {
                printf(",-1");
            }
            printf("\"\nwhere T is a number between 0 and %d,\n", nsteps - 1);
            printf("or dump each block separately with \n"
                   "    bpls -dD %s\n",
                   variable->m_Name.c_str());
            return -1;
        }
    }

    for (j = 0; j < ndim; j++)
    {
        uint64_t st, ct;
        if (istart[j + tidx] < 0) // negative index means last-|index|
            st = shape[j] + istart[j + tidx];
        else
            st = istart[j + tidx];
        if (icount[j + tidx] < 0) // negative index means last-|index|+1-start
            ct = shape[j] + icount[j + tidx] + 1 - st;
        else
            ct = icount[j + tidx];

        if (verbose > 2)
            printf("    j=%d, st=%" PRIu64 " ct=%" PRIu64 "\n", j + tidx, st, ct);

        start_t[j + tidx] = st;
        count_t[j + tidx] = ct;
        nelems *= ct;
        if (verbose > 1)
            printf("    s[%d]=%" PRIu64 ", c[%d]=%" PRIu64 ", n=%" PRIu64 "\n", j + tidx,
                   start_t[j + tidx], j + tidx, count_t[j + tidx], nelems);
    }

    if (verbose > 1)
    {
        printf(" total size of data to read = %" PRIu64 "\n", nelems * elemsize);
    }

    print_slice_info(variable, (tidx == 1), start_t, count_t, shape);

    maxreadn = (uint64_t)MAX_BUFFERSIZE / elemsize;
    if (nelems < maxreadn)
        maxreadn = nelems;

    bool xmlprint = helper::EndsWith(variable->m_Name, ".xml", false);
    if (xmlprint && nelems > maxreadn)
        maxreadn = nelems;

    // special case: string. Need to use different elemsize
    /*if (vi->type == DataType::String)
    {
        if (vi->value)
            elemsize = strlen(vi->value) + 1;
        maxreadn = elemsize;
    }*/

    // allocate data array
    // data = (T *)malloc(maxreadn * elemsize);

    // determine strategy how to read in:
    //  - at once
    //  - loop over 1st dimension
    //  - loop over 1st & 2nd dimension
    //  - etc
    if (verbose > 1)
        printf("Read size strategy:\n");
    sum = (uint64_t)1;
    actualreadn = (uint64_t)1;
    for (i = tdims - 1; i >= 0; i--)
    {
        if (sum >= (uint64_t)maxreadn)
        {
            readn[i] = 1;
        }
        else
        {
            readn[i] = maxreadn / (int)sum; // sum is small for 4 bytes here
            // this may be over the max count for this dimension
            if (readn[i] > count_t[i])
                readn[i] = count_t[i];
        }
        if (verbose > 1)
            printf("    dim %d: read %" PRIu64 " elements\n", i, readn[i]);
        sum = sum * (uint64_t)count_t[i];
        actualreadn = actualreadn * readn[i];
    }
    if (verbose > 1)
        printf("    read %" PRIu64 " elements at once, %" PRIu64 " in total (nelems=%" PRIu64 ")\n",
               actualreadn, sum, nelems);

    // init s and c
    // and calculate ndigits_dims
    for (j = 0; j < tdims; j++)
    {
        s[j] = start_t[j];
        c[j] = readn[j];

        ndigits_dims[j] =
            ndigits(start_t[j] + count_t[j] - 1); // -1: dim=100 results in 2 digits (0..99)
    }

    // read until read all 'nelems' elements
    sum = 0;
    while (sum < nelems)
    {

        // how many elements do we read in next?
        actualreadn = 1;
        for (j = 0; j < tdims; j++)
            actualreadn *= c[j];

        if (verbose > 2)
        {
            printf("adios_read_var name=%s ", variable->m_Name.c_str());
            PRINT_DIMS_UINT64("  start", s, tdims, j);
            PRINT_DIMS_UINT64("  count", c, tdims, j);
            printf("  read %" PRIu64 " elems\n", actualreadn);
        }

        // read a slice finally
        const Dims startv = variable->m_ShapeID == ShapeID::GlobalArray
                                ? helper::Uint64ArrayToSizetVector(tdims - tidx, s + tidx)
                                : Dims();
        const Dims countv = variable->m_ShapeID == ShapeID::GlobalArray
                                ? helper::Uint64ArrayToSizetVector(tdims - tidx, c + tidx)
                                : Dims();

        if (verbose > 2)
        {
            printf("set selection: ");
            PRINT_DIMS_SIZET("  start", startv.data(), tdims - tidx, j);
            PRINT_DIMS_SIZET("  count", countv.data(), tdims - tidx, j);
            printf("\n");
        }

        if (variable->m_ShapeID == ShapeID::GlobalArray)
        {
            variable->SetSelection({startv, countv});
        }

        if (tidx)
        {
            if (verbose > 2)
            {
                printf("set Step selection: from relative step %" PRIu64 " read %" PRIu64
                       " steps\n",
                       s[0], c[0]);
            }
            variable->SetStepSelection({s[0], c[0]});
        }

        dataV.resize(variable->SelectionSize());
        fp->Get(*variable, dataV, adios2::Mode::Sync);

        // print slice
        print_dataset(dataV.data(), variable->m_Type, s, c, tdims, ndigits_dims);

        // prepare for next read
        sum += actualreadn;
        bool incdim = true; // last dim should be increased
        for (j = tdims - 1; j >= 0; j--)
        {
            if (incdim)
            {
                if (s[j] + c[j] == start_t[j] + count_t[j])
                {
                    // reached the end of this dimension
                    s[j] = start_t[j];
                    c[j] = readn[j];
                    incdim = true; // previous dim can increase too
                }
                else
                {
                    // move up in this dimension up to total count
                    s[j] += readn[j];
                    if (s[j] + c[j] > start_t[j] + count_t[j])
                    {
                        // do not reach over the limit
                        c[j] = start_t[j] + count_t[j] - s[j];
                    }
                    incdim = false;
                }
            }
        }
    } // end while sum < nelems
    print_endline();

    if (accuracyWasSet)
    {
        adios2::Accuracy a = variable->GetAccuracy();
        std::cout << "Read accuracy was (error = " << a.error << ", norm = " << a.norm << ", "
                  << (a.relative ? "rel" : "abs") << ")\n";
    }

    return 0;
}

/** Read one writeblock of a variable and print
 * Return: 0: ok, != 0 on error
 */
template <class T>
int readVarBlock(core::Engine *fp, core::IO *io, core::Variable<T> *variable, size_t step,
                 size_t blockid, Dims Count, Dims Start)
{
    int i, j;
    uint64_t start_t[MAX_DIMS],
        count_t[MAX_DIMS]; // processed <0 values in start/count
    uint64_t s[MAX_DIMS],
        c[MAX_DIMS]; // for block reading of smaller chunks
    uint64_t nelems; // number of elements to read
    int tidx;
    uint64_t st, ct;
    std::vector<T> dataV;
    uint64_t sum;             // working var to sum up things
    uint64_t maxreadn;        // max number of elements to read once up to a limit
                              // (10MB of data)
    uint64_t actualreadn;     // our decision how much to read at once
    uint64_t readn[MAX_DIMS]; // how big chunk to read in in each dimension?
    bool incdim;              // used in incremental reading in
    int ndigits_dims[32];     // # of digits (to print) of each dimension

    const size_t elemsize = variable->m_ElementSize;
    const int nsteps = (timestep ? 1 : static_cast<int>(variable->GetAvailableStepsCount()));
    const int ndim = static_cast<int>(variable->m_Count.size());
    // create the counter arrays with the appropriate lengths
    // transfer start and count arrays to format dependent arrays

    nelems = 1;
    tidx = 0;

    if (nsteps > 1)
    {
        if (istart[0] < 0) // negative index means last-|index|
            st = nsteps + istart[0];
        else
            st = istart[0];
        if (icount[0] < 0) // negative index means last-|index|+1-start
            ct = nsteps + icount[0] + 1 - st;
        else
            ct = icount[0];

        if (verbose > 2)
            printf("    time: st=%" PRIu64 " ct=%" PRIu64 "\n", st, ct);

        // check if this block falls into the requested time
        if (step < st || step >= st + ct)
            return 0;
        tidx = 1;
    }

    int out_of_bound = 0;
    for (j = 0; j < ndim; j++)
    {
        if (istart[j + tidx] < 0) // negative index means last-|index|
            st = Count[j] + istart[j + tidx];
        else
            st = istart[j + tidx];
        if (icount[j + tidx] < 0) // negative index means last-|index|+1-start
            ct = Count[j] + icount[j + tidx] + 1 - st;
        else
            ct = icount[j + tidx];

        if (st > Count[j])
        {
            out_of_bound = 1;
        }
        else if (ct > Count[j] - st)
        {
            ct = Count[j] - st;
        }
        if (verbose > 2)
            printf("    j=%d, st=%" PRIu64 " ct=%" PRIu64 "\n", j, st, ct);

        start_t[j] = st;
        count_t[j] = ct;
        nelems *= ct;
        if (verbose > 1)
            printf("    s[%d]=%" PRIu64 ", c[%d]=%" PRIu64 ", n=%" PRIu64 "\n", j, start_t[j], j,
                   count_t[j], nelems);
    }

    if (verbose > 1)
    {
        printf(" total size of data to read = %" PRIu64 "\n", nelems * elemsize);
    }

    print_slice_info(variable, false, start_t, count_t, Count);

    if (out_of_bound)
        return 0;

    maxreadn = MAX_BUFFERSIZE / elemsize;
    if (nelems < maxreadn)
        maxreadn = nelems;

    // allocate data array
    // data = (void *)malloc(maxreadn * elemsize + 8); // +8 for just to be
    // sure

    // determine strategy how to read in:
    //  - at once
    //  - loop over 1st dimension
    //  - loop over 1st & 2nd dimension
    //  - etc
    if (verbose > 1)
        printf("Read size strategy:\n");
    sum = (uint64_t)1;
    actualreadn = (uint64_t)1;
    for (i = ndim - 1; i >= 0; i--)
    {
        if (sum >= (uint64_t)maxreadn)
        {
            readn[i] = 1;
        }
        else
        {
            readn[i] = maxreadn / (int)sum; // sum is small for 4 bytes here
            // this may be over the max count for this dimension
            if (readn[i] > count_t[i])
                readn[i] = count_t[i];
        }
        if (verbose > 1)
            printf("    dim %d: read %" PRIu64 " elements\n", i, readn[i]);
        sum = sum * (uint64_t)count_t[i];
        actualreadn = actualreadn * readn[i];
    }
    if (verbose > 1)
        printf("    read %" PRIu64 " elements at once, %" PRIu64 " in total (nelems=%" PRIu64 ")\n",
               actualreadn, sum, nelems);

    // init s and c
    // and calculate ndigits_dims
    for (j = 0; j < ndim; j++)
    {
        s[j] = start_t[j];
        c[j] = readn[j];

        ndigits_dims[j] =
            ndigits(start_t[j] + count_t[j] - 1); // -1: dim=100 results in 2 digits (0..99)
    }

    // read until read all 'nelems' elements
    sum = 0;
    while (sum < nelems)
    {

        // how many elements do we read in next?
        actualreadn = 1;
        for (j = 0; j < ndim; j++)
            actualreadn *= c[j];

        uint64_t startoffset = s[ndim - 1];
        uint64_t tmpprod = c[ndim - 1];
        for (i = ndim - 2; i >= 0; i--)
        {
            startoffset += s[i] * tmpprod;
            tmpprod *= c[i];
        }

        if (verbose > 2)
        {
            printf("adios_read_var name=%s ", variable->m_Name.c_str());
            PRINT_DIMS_UINT64("  start", s, ndim, j);
            PRINT_DIMS_UINT64("  count", c, ndim, j);
            printf("  read %" PRIu64 " elems\n", actualreadn);
        }
        if (verbose > 1)
            printf("    read block %zu from offset %" PRIu64 " nelems %" PRIu64 ")\n", blockid,
                   startoffset, actualreadn);

        // read a slice finally
        Dims startv = helper::Uint64ArrayToSizetVector(ndim, s);
        Dims countv = helper::Uint64ArrayToSizetVector(ndim, c);

        /* In current implementation we read with global selection, so
         * we need to adjust start_t for global offsets here.
         * TODO: this will change in the future to block reading with
         * relative
         * selection
         */
        if (variable->m_ShapeID == ShapeID::GlobalArray)
        {
            for (j = 0; j < ndim; j++)
            {
                startv[j] += Start[j];
            }
        }

        if (verbose > 2)
        {
            printf("set selection: ");
            PRINT_DIMS_SIZET("  start", startv.data(), ndim, j);
            PRINT_DIMS_SIZET("  count", countv.data(), ndim, j);
            printf("\n");
        }

        if (!variable->m_SingleValue)
        {
            if (variable->m_ShapeID == ShapeID::LocalArray)
            {
                variable->SetBlockSelection(blockid);
            }
            variable->SetSelection({startv, countv});
        }

        if (nsteps > 1)
        {
            if (verbose > 2)
            {
                printf("set Step selection: from %" PRIu64 " read %" PRIu64 " steps\n", s[0], c[0]);
            }
            variable->SetStepSelection({step, 1});
        }

        dataV.resize(variable->SelectionSize());
        fp->Get(*variable, dataV, adios2::Mode::Sync);
        // print slice
        print_dataset(dataV.data(), variable->m_Type, s, c, ndim, ndigits_dims);

        // prepare for next read
        sum += actualreadn;
        incdim = true; // largest dim should be increased
        for (j = ndim - 1; j >= 0; j--)
        {
            if (incdim)
            {
                if (s[j] + c[j] == start_t[j] + count_t[j])
                {
                    // reached the end of this dimension
                    s[j] = start_t[j];
                    c[j] = readn[j];
                    incdim = true; // next smaller dim can increase too
                }
                else
                {
                    // move up in this dimension up to total count
                    s[j] += readn[j];
                    if (s[j] + c[j] > start_t[j] + count_t[j])
                    {
                        // do not reach over the limit
                        c[j] = start_t[j] + count_t[j] - s[j];
                    }
                    incdim = false;
                }
            }
        }
    } // end while sum < nelems
    print_endline();
    return 0;
}

bool matchesAMask(const char *name)
{
    int startpos = 0; // to match with starting / or without
#ifdef USE_C_REGEX
    regmatch_t pmatch[1] = {{(regoff_t)-1, (regoff_t)-1}};
#else
#endif

    if (nmasks == 0)
        return true;

    for (int i = 0; i < nmasks; i++)
    {
        if (use_regexp)
        {
#ifdef USE_C_REGEX
            int excode = regexec(&(varregex[i]), name, 1, pmatch, 0);
            if (name[0] == '/') // have to check if it matches from the
                                // second character too
                startpos = 1;
            if (excode == 0 &&                                           // matches
                (pmatch[0].rm_so == 0 || pmatch[0].rm_so == startpos) && // from the beginning
                static_cast<size_t>(pmatch[0].rm_eo) == strlen(name) // to the very end of the name
            )
#else
            bool matches = std::regex_match(name, varregex[i]);
            if (!matches && name[0] == '/')
                matches = std::regex_match(name + 1, varregex[i]);
            if (matches)
#endif
            {
                if (verbose > 1)
                    printf("Name %s matches regexp %i %s\n", name, i, varmask[i]);
                return true;
            }
        }
        else
        {
            // use shell pattern matching
            if (varmask[i][0] != '/' && name[0] == '/')
                startpos = 1;
#ifdef _WIN32
            if (PathMatchSpec(name + startpos, varmask[i]))
#else
            if (fnmatch(varmask[i], name + startpos, FNM_FILE_NAME) == 0)
#endif
            {
                if (verbose > 1)
                    printf("Name %s matches varmask %i %s\n", name, i, varmask[i]);
                return true;
            }
        }
    }
    return false;
}

int print_start(const std::string &fname)
{
    if (fname.empty())
    {
        outf = stdout;
    }
    else
    {
        if ((outf = fopen(fname.c_str(), "w")) == NULL)
        {
            fprintf(stderr, "Error at opening for writing file %s: %s\n", fname.c_str(),
                    strerror(errno));
            return 30;
        }
    }
    return 0;
}

void print_stop() { fclose(outf); }

static int nextcol = 0; // column index to start with (can have lines split in two calls)

void print_slice_info(core::VariableBase *variable, bool timed, uint64_t *s, uint64_t *c,
                      Dims count)
{
    // print the slice info in indexing is on and
    // not the complete variable is read
    size_t ndim = variable->m_Shape.size();
    size_t nsteps = variable->m_AvailableStepsCount;
    bool isaslice = false;
    int tidx = (timed ? 1 : 0);
    size_t tdim = ndim + tidx;
    if (timed)
    {
        if (c[0] < nsteps)
            isaslice = true;
    }
    for (size_t i = 0; i < ndim; i++)
    {
        if (c[i + tidx] < count[i])
            isaslice = true;
    }
    if (isaslice)
    {
        fprintf(outf, "%c   slice (%" PRIu64 ":%" PRIu64, commentchar, s[0], s[0] + c[0] - 1);
        for (size_t i = 1; i < tdim; i++)
        {
            fprintf(outf, ", %" PRIu64 ":%" PRIu64, s[i], s[i] + c[i] - 1);
        }
        fprintf(outf, ")\n");
    }
}

int print_data_as_string(const void *data, int maxlen, DataType adiosvartype)
{
    const char *str = (const char *)data;
    int len = maxlen;
    switch (adiosvartype)
    {
    case DataType::UInt8:
    case DataType::Int8:
    case DataType::String:
        while (str[len - 1] == 0)
        {
            len--;
        } // go backwards on ascii 0s
        if (len < maxlen)
        {
            // it's a C string with terminating \0
            fprintf(outf, "\"%s\"", str);
        }
        else
        {
            // fortran VARCHAR, lets trim from right padded zeros
            while (str[len - 1] == ' ')
            {
                len--;
            }
            fprintf(outf, "\"%*.*s\"", len, len, (char *)data);
            if (len < maxlen)
                fprintf(outf, " + %d spaces", maxlen - len);
        }
        break;
    default:
        fprintf(stderr,
                "Error in bpls code: cannot use print_data_as_string() "
                "for type \"%d\"\n",
                static_cast<typename std::underlying_type<DataType>::type>(adiosvartype));
        return -1;
    }
    return 0;
}

int print_data_characteristics(void *min, void *max, double *avg, double *std_dev,
                               DataType adiosvartype, bool allowformat)
{
    bool f = format.size() && allowformat;
    const char *fmt = format.c_str();

    switch (adiosvartype)
    {
    case DataType::UInt8:
        if (min)
            fprintf(outf, (f ? fmt : "%10hhu  "), *((unsigned char *)min));
        else
            fprintf(outf, "      null  ");
        if (max)
            fprintf(outf, (f ? fmt : "%10hhu  "), *((unsigned char *)max));
        else
            fprintf(outf, "      null  ");
        if (avg)
            fprintf(outf, "%10.2f  ", *avg);
        else
            fprintf(outf, "      null  ");
        if (std_dev)
            fprintf(outf, "%10.2f  ", *std_dev);
        else
            fprintf(outf, "      null  ");
        break;
    case DataType::Int8:
        if (min)
            fprintf(outf, (f ? fmt : "%10hhd  "), *((char *)min));
        else
            fprintf(outf, "      null  ");
        if (max)
            fprintf(outf, (f ? fmt : "%10hhd  "), *((char *)max));
        else
            fprintf(outf, "      null  ");
        if (avg)
            fprintf(outf, "%10.2f  ", *avg);
        else
            fprintf(outf, "      null  ");
        if (std_dev)
            fprintf(outf, "%10.2f  ", *std_dev);
        else
            fprintf(outf, "      null  ");
        break;
    case DataType::String:
        break;

    case DataType::UInt16:
        if (min)
            fprintf(outf, (f ? fmt : "%10hu  "), (*(unsigned short *)min));
        else
            fprintf(outf, "      null  ");
        if (max)
            fprintf(outf, (f ? fmt : "%10hu  "), (*(unsigned short *)max));
        else
            fprintf(outf, "      null  ");
        if (avg)
            fprintf(outf, "%10.2f  ", *avg);
        else
            fprintf(outf, "      null  ");
        if (std_dev)
            fprintf(outf, "%10.2f  ", *std_dev);
        else
            fprintf(outf, "      null  ");
        break;
    case DataType::Int16:
        if (min)
            fprintf(outf, (f ? fmt : "%10hd  "), (*(short *)min));
        else
            fprintf(outf, "      null  ");
        if (max)
            fprintf(outf, (f ? fmt : "%10hd  "), (*(short *)max));
        else
            fprintf(outf, "      null  ");
        if (avg)
            fprintf(outf, "%10.2f  ", *avg);
        else
            fprintf(outf, "      null  ");
        if (std_dev)
            fprintf(outf, "%10.2f  ", *std_dev);
        else
            fprintf(outf, "      null  ");
        break;

    case DataType::UInt32:
        if (min)
            fprintf(outf, (f ? fmt : "%10u  "), (*(unsigned int *)min));
        else
            fprintf(outf, "      null  ");
        if (max)
            fprintf(outf, (f ? fmt : "%10u  "), (*(unsigned int *)max));
        else
            fprintf(outf, "      null  ");
        if (avg)
            fprintf(outf, "%10.2f  ", *avg);
        else
            fprintf(outf, "      null  ");
        if (std_dev)
            fprintf(outf, "%10.2f  ", *std_dev);
        else
            fprintf(outf, "      null  ");
        break;
    case DataType::Int32:
        if (min)
            fprintf(outf, (f ? fmt : "%10d  "), (*(int *)min));
        else
            fprintf(outf, "      null  ");
        if (max)
            fprintf(outf, (f ? fmt : "%10d  "), (*(int *)max));
        else
            fprintf(outf, "      null  ");
        if (avg)
            fprintf(outf, "%10.2f  ", *avg);
        else
            fprintf(outf, "      null  ");
        if (std_dev)
            fprintf(outf, "%10.2f  ", *std_dev);
        else
            fprintf(outf, "      null  ");
        break;

    case DataType::UInt64:
        if (min)
            fprintf(outf, (f ? fmt : "%10llu  "), (*(unsigned long long *)min));
        else
            fprintf(outf, "      null  ");
        if (max)
            fprintf(outf, (f ? fmt : "%10llu  "), (*(unsigned long long *)max));
        else
            fprintf(outf, "      null  ");
        if (avg)
            fprintf(outf, "%10.2f  ", *avg);
        else
            fprintf(outf, "      null  ");
        if (std_dev)
            fprintf(outf, "%10.2f  ", *std_dev);
        else
            fprintf(outf, "      null  ");
        break;
    case DataType::Int64:
        if (min)
            fprintf(outf, (f ? fmt : "%10lld  "), (*(long long *)min));
        else
            fprintf(outf, "      null  ");
        if (max)
            fprintf(outf, (f ? fmt : "%10lld  "), (*(long long *)max));
        else
            fprintf(outf, "      null  ");
        if (avg)
            fprintf(outf, "%10.2f  ", *avg);
        else
            fprintf(outf, "      null  ");
        if (std_dev)
            fprintf(outf, "%10.2f  ", *std_dev);
        else
            fprintf(outf, "      null  ");
        break;

    case DataType::Float:
        if (min)
            fprintf(outf, (f ? fmt : "%10.2g  "), (*(float *)min));
        else
            fprintf(outf, "      null  ");
        if (max)
            fprintf(outf, (f ? fmt : "%10.2g  "), (*(float *)max));
        else
            fprintf(outf, "      null  ");
        if (avg)
            fprintf(outf, "%10.2g  ", *avg);
        else
            fprintf(outf, "      null  ");
        if (std_dev)
            fprintf(outf, "%10.2g  ", *std_dev);
        else
            fprintf(outf, "      null  ");
        break;
    case DataType::Double:
        if (min)
            fprintf(outf, (f ? fmt : "%10.2g  "), (*(double *)min));
        else
            fprintf(outf, "      null  ");
        if (max)
            fprintf(outf, (f ? fmt : "%10.2g  "), (*(double *)max));
        else
            fprintf(outf, "      null  ");
        if (avg)
            fprintf(outf, "%10.2g  ", *avg);
        else
            fprintf(outf, "      null  ");
        if (std_dev)
            fprintf(outf, "%10.2g  ", *std_dev);
        else
            fprintf(outf, "      null  ");
        break;

    case DataType::LongDouble:
        fprintf(outf, "????????");
        break;

    // TO DO
    /*
       case DataType::FloatComplex:
       fprintf(outf,(f ? format : "(%g,i%g) "), ((float *) data)[2*item],
       ((float *) data)[2*item+1]);
       break;

       case DataType::DoubleComplex:
       fprintf(outf,(f ? format : "(%g,i%g)" ), ((double *) data)[2*item],
       ((double *) data)[2*item+1]);
       break;
     */
    default:
        break;
    } // end switch
    return 0;
}

/* s is a character array not necessarily null terminated.
 * return false on OK print, true if it not XML (not printed)*/
bool print_data_xml(const char *s, const size_t length)
{
    pugi::xml_document document;
    auto parse_result = document.load_buffer(s, length);
    if (!parse_result)
    {
        //   helper::Throw<std::invalid_argument>( "Utils", "bpls",
        //   "print_data_xml", "ERROR: XML: parse error in XML string,
        //   description: " + std::string(parse_result.description()) + ", check
        //   with any XML editor if format is ill-formed");
        return true;
    }
    std::cout << std::endl;
    document.save(std::cout, PUGIXML_TEXT("  "),
                  pugi::format_default | pugi::format_no_declaration);
    std::cout << std::flush;
    return false;
}

int print_data(const void *data, int item, DataType adiosvartype, bool allowformat,
               bool char_star_string)
{
    bool f = format.size() && allowformat;
    const char *fmt = format.c_str();
    if (data == NULL)
    {
        fprintf(outf, "null ");
        return 0;
    }
    // print next data item
    switch (adiosvartype)
    {
    case DataType::Char:
        fprintf(outf, (f ? fmt : "%c"), ((char *)data)[item]);
        break;
    case DataType::UInt8:
        fprintf(outf, (f ? fmt : "%hhu"), ((unsigned char *)data)[item]);
        break;
    case DataType::Int8:
        fprintf(outf, (f ? fmt : "%hhd"), ((signed char *)data)[item]);
        break;

    case DataType::String: {
        if (char_star_string)
        {
            fprintf(outf, (f ? fmt : "\"%s\""), *((char **)data));
        }
        else
        {
            const std::string *dataStr = reinterpret_cast<const std::string *>(data);
            fprintf(outf, (f ? fmt : "\"%s\""), dataStr[item].c_str());
        }
        break;
    }

    case DataType::UInt16:
        fprintf(outf, (f ? fmt : "%hu"), ((unsigned short *)data)[item]);
        break;
    case DataType::Int16:
        fprintf(outf, (f ? fmt : "%hd"), ((signed short *)data)[item]);
        break;

    case DataType::UInt32:
        fprintf(outf, (f ? fmt : "%u"), ((unsigned int *)data)[item]);
        break;
    case DataType::Int32:
        fprintf(outf, (f ? fmt : "%d"), ((signed int *)data)[item]);
        break;

    case DataType::UInt64:
        fprintf(outf, (f ? fmt : "%llu"), ((unsigned long long *)data)[item]);
        break;
    case DataType::Int64:
        fprintf(outf, (f ? fmt : "%lld"), ((signed long long *)data)[item]);
        break;

    case DataType::Float:
        fprintf(outf, (f ? fmt : "%g"), ((float *)data)[item]);
        break;

    case DataType::Double:
        fprintf(outf, (f ? fmt : "%g"), ((double *)data)[item]);
        break;

    case DataType::LongDouble:
        fprintf(outf, (f ? fmt : "%Lg"), ((long double *)data)[item]);
        // fprintf(outf,(f ? fmt : "????????"));
        break;

    case DataType::FloatComplex:
        fprintf(outf, (f ? fmt : "(%g,i%g)"), ((float *)data)[2 * item],
                ((float *)data)[2 * item + 1]);
        break;

    case DataType::DoubleComplex:
        fprintf(outf, (f ? fmt : "(%g,i%g)"), ((double *)data)[2 * item],
                ((double *)data)[2 * item + 1]);
        break;

    default:
        break;
    } // end switch
    return 0;
}

int print_dataset(const void *data, const DataType vartype, uint64_t *s, uint64_t *c, int tdims,
                  int *ndigits)
{
    int i, item, steps;
    char idxstr[128], buf[16];
    uint64_t ids[MAX_DIMS]; // current indices
    bool roll;
    DataType adiosvartype = vartype;

    // init current indices
    steps = 1;
    for (i = 0; i < tdims; i++)
    {
        ids[i] = s[i];
        steps *= c[i];
    }

    item = 0; // index to *data
    // loop through each data item and print value
    while (item < steps)
    {

        // print indices if needed into idxstr;
        idxstr[0] = '\0'; // empty idx string
        if (nextcol == 0)
        {
            if (!noindex && tdims > 0)
            {
                snprintf(idxstr, sizeof(idxstr), "    (%*" PRIu64, ndigits[0], ids[0]);
                for (i = 1; i < tdims; i++)
                {
                    snprintf(buf, sizeof(buf), ",%*" PRIu64, ndigits[i], ids[i]);
                    strcat(idxstr, buf);
                }
                strcat(idxstr, ")    ");
            }
        }

        // print item
        fprintf(outf, "%s", idxstr);
        if (printByteAsChar && (adiosvartype == DataType::Int8 || adiosvartype == DataType::UInt8))
        {
            /* special case: k-D byte array printed as (k-1)D array of
             * strings
             */
            if (tdims == 0)
            {
                print_data_as_string(data, steps, adiosvartype);
            }
            else
            {
                print_data_as_string((char *)data + item, c[tdims - 1],
                                     adiosvartype); // print data of last dim as string
                item += c[tdims - 1] - 1;           // will be ++-ed once below
                ids[tdims - 1] = s[tdims - 1] + c[tdims - 1] - 1; // will be rolled below
            }
            nextcol = ncols - 1; // force new line, will be ++-ed once below
        }
        else
        {
            print_data(data, item, adiosvartype, true);
        }

        // increment/reset column index
        nextcol++;
        if (nextcol == ncols)
        {
            fprintf(outf, "\n");
            nextcol = 0;
        }
        else
        {
            fprintf(outf, " ");
        }

        // increment indices
        item++;
        roll = true;
        for (i = tdims - 1; i >= 0; i--)
        {
            if (roll)
            {
                if (ids[i] == s[i] + c[i] - 1)
                {
                    // last index in this dimension, roll upward
                    ids[i] = s[i];
                }
                else
                {
                    ids[i]++;
                    roll = false;
                }
            }
        }
    }
    return 0;
}

void print_endline(void)
{
    if (nextcol != 0)
        fprintf(outf, "\n");
    nextcol = 0;
}

template <class T>
size_t relative_to_absolute_step(core::Engine *fp, core::Variable<T> *variable,
                                 const size_t relstep)
{
    auto minBlocks = fp->MinBlocksInfo(*variable, relstep);

    if (minBlocks)
    {
        size_t Step = minBlocks->Step;
        delete minBlocks;
        return Step;
    }

    const std::map<size_t, std::vector<size_t>> &indices =
        variable->m_AvailableStepBlockIndexOffsets;
    auto itStep = indices.begin();
    size_t absstep = itStep->first - 1;

    for (size_t step = 0; step < relstep; step++)
    {
        ++itStep;
        absstep = itStep->first - 1;
    }
    return absstep;
}

template <class T>
Dims get_global_array_signature(core::Engine *fp, core::IO *io, core::Variable<T> *variable)
{
    const size_t ndim = variable->m_Shape.size();
    Dims dims(ndim, 0);
    if (timestep)
    {
        dims = variable->Shape();
    }
    else
    {
        const size_t nsteps = variable->GetAvailableStepsCount();
        bool firstStep = true;

        // looping over the absolute step indexes
        // is not supported by a simple API function
        auto minBlocks = fp->MinBlocksInfo(*variable, fp->CurrentStep());

        if (minBlocks)
        {
            for (size_t step = 0; step < nsteps; step++)
            {
                if (step > 0)
                {
                    minBlocks = fp->MinBlocksInfo(*variable, step);
                }
                if (minBlocks->Shape)
                {
                    for (size_t k = 0; k < ndim; k++)
                    {
                        size_t n =
                            (minBlocks->WasLocalValue ? reinterpret_cast<size_t>(minBlocks->Shape)
                                                      : minBlocks->Shape[k]);
                        if (firstStep)
                        {
                            dims[k] = n;
                        }
                        else if (dims[k] != n)
                        {
                            dims[k] = 0;
                        }
                    }
                    firstStep = false;
                }
                delete minBlocks;
            }
        }
        else
        {
            const std::map<size_t, std::vector<size_t>> &indices =
                variable->m_AvailableStepBlockIndexOffsets;
            auto itStep = indices.begin();
            if (itStep != indices.end())
            {
                for (size_t step = 0; step < nsteps; step++)
                {
                    const size_t absstep = itStep->first;
                    Dims d = variable->Shape(absstep - 1);
                    if (d.empty())
                    {
                        continue;
                    }

                    for (size_t k = 0; k < ndim; k++)
                    {
                        if (firstStep)
                        {
                            dims[k] = d[k];
                        }
                        else if (dims[k] != d[k])
                        {
                            dims[k] = 0;
                        }
                    }
                    firstStep = false;
                    ++itStep;
                }
            }
        }
    }
    return dims;
}

template <class T>
std::pair<size_t, Dims> get_local_array_signature(core::Engine *fp, core::IO *io,
                                                  core::Variable<T> *variable)
{
    const size_t ndim = variable->m_Count.size();
    size_t nblocks = 0;
    Dims dims(ndim, 0);

    if (timestep)
    {
        const auto minBlocks = fp->MinBlocksInfo(*variable, fp->CurrentStep());

        if (minBlocks)
        {
            bool firstBlock = true;
            nblocks = minBlocks->BlocksInfo.size();
            for (size_t j = 0; j < nblocks; j++)
            {
                for (size_t k = 0; k < ndim; k++)
                {
                    if (firstBlock)
                    {
                        dims[k] = minBlocks->BlocksInfo[j].Count[k];
                    }
                    else if (dims[k] != minBlocks->BlocksInfo[j].Count[k])
                    {
                        dims[k] = 0;
                    }
                }
                firstBlock = false;
            }
        }
        else
        {
            std::vector<typename core::Variable<T>::BPInfo> blocks =
                fp->BlocksInfo(*variable, fp->CurrentStep());
            if (!blocks.empty())
            {
                nblocks = blocks.size();
                bool firstBlock = true;
                for (size_t j = 0; j < nblocks; j++)
                {
                    for (size_t k = 0; k < ndim; k++)
                    {
                        if (firstBlock)
                        {
                            dims[k] = blocks[j].Count[k];
                        }
                        else if (dims[k] != blocks[j].Count[k])
                        {
                            dims[k] = 0;
                        }
                    }
                    firstBlock = false;
                }
            }
        }
    }
    else
    {
        bool firstStep = true;
        bool firstBlock = true;

        MinVarInfo *minBlocksInfo = nullptr;
        minBlocksInfo = fp->MinBlocksInfo(*variable, 0 /* relative step 0 */);

        // first step
        if (minBlocksInfo)
        {
            dims.resize(minBlocksInfo->Dims);
            size_t RelStep = 0;
            for (RelStep = 0; RelStep < variable->m_AvailableStepsCount; RelStep++)
            {
                if (RelStep > 0)
                {
                    minBlocksInfo = fp->MinBlocksInfo(*variable, RelStep);
                }

                auto &coreBlocksInfo = minBlocksInfo->BlocksInfo;
                if (firstStep)
                {
                    nblocks = coreBlocksInfo.size();
                }
                else if (nblocks != coreBlocksInfo.size())
                {
                    nblocks = 0;
                }
                for (auto &coreBlockInfo : coreBlocksInfo)
                {
                    if (firstBlock)
                    {
                        for (size_t k = 0; k < dims.size(); k++)
                        {
                            dims[k] = coreBlockInfo.Count[k];
                        }
                    }
                    else
                    {
                        for (size_t k = 0; k < dims.size(); k++)
                        {
                            if (dims[k] != coreBlockInfo.Count[k])
                            {
                                dims[k] = 0;
                            }
                        }
                    }
                    firstBlock = false;
                }
                firstStep = false;
            }
        }
        else
        {
            std::map<size_t, std::vector<typename core::Variable<T>::BPInfo>> allblocks =
                fp->AllStepsBlocksInfo(*variable);

            for (auto &blockpair : allblocks)
            {
                std::vector<typename adios2::core::Variable<T>::BPInfo> &blocks = blockpair.second;
                const size_t blocksSize = blocks.size();
                if (firstStep)
                {
                    nblocks = blocksSize;
                }
                else if (nblocks != blocksSize)
                {
                    nblocks = 0;
                }

                for (size_t j = 0; j < blocksSize; j++)
                {
                    for (size_t k = 0; k < ndim; k++)
                    {
                        if (firstBlock)
                        {
                            dims[k] = blocks[j].Count[k];
                        }
                        else if (dims[k] != blocks[j].Count[k])
                        {
                            dims[k] = 0;
                        }
                    }
                    firstBlock = false;
                }
                firstStep = false;
            }
        }
    }
    return std::make_pair(nblocks, dims);
}

static int ndigits_dims[32] = {
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
};
template <class T>
void print_decomp(core::Engine *fp, core::IO *io, core::Variable<T> *variable)
{
    /* Print block info */
    DataType adiosvartype = variable->m_Type;

    MinVarInfo *mBI = nullptr;
    mBI = fp->MinBlocksInfo(*variable, variable->m_AvailableStepsCount - 1 /* relative step 0 */);

    // first step
    if (mBI)
    {
        size_t laststep = mBI->Step; // used relative last step above
        delete mBI;
        int ndigits_nsteps = ndigits(laststep);
        if (variable->m_ShapeID == ShapeID::GlobalValue ||
            variable->m_ShapeID == ShapeID::LocalValue)
        {
            for (size_t RelStep = 0; RelStep < variable->m_AvailableStepsCount; RelStep++)
            {
                auto minBlocksInfo = fp->MinBlocksInfo(*variable, RelStep);
                auto blocks = minBlocksInfo->BlocksInfo;
                fprintf(outf, "%c       step %*zu: ", commentchar, ndigits_nsteps,
                        minBlocksInfo->Step);
                if (blocks.size() == 1)
                {
                    fprintf(outf, " = ");
                    print_data(blocks[0].BufferP, 0, adiosvartype, true, /* MBI */ true);
                    fprintf(outf, "\n");
                }
                else
                {
                    fprintf(outf, "%zu instances available\n", blocks.size());
                }
                if (dump)
                {
                    fprintf(outf, "               ");
                    int col = 0;
                    for (size_t j = 0; j < blocks.size(); j++)
                    {
                        print_data(blocks[j].BufferP, 0, adiosvartype, true, /* MBI */ true);
                        ++col;
                        if (j < blocks.size() - 1)
                        {
                            if (col < ncols)
                            {
                                fprintf(outf, " ");
                            }
                            else
                            {
                                fprintf(outf, "\n               ");
                                col = 0;
                            }
                        }
                    }
                    fprintf(outf, "\n");
                }
                delete minBlocksInfo;
            }
        }
        else
        {
            // arrays
            for (size_t RelStep = 0; RelStep < variable->m_AvailableStepsCount; RelStep++)
            {
                auto minBlocksInfo = fp->MinBlocksInfo(*variable, RelStep);
                auto blocks = minBlocksInfo->BlocksInfo;
                size_t ndim = variable->m_Count.size();
                int ndigits_nblocks;
                for (size_t k = 0; k < ndim; k++)
                {
                    // get digit lengths for each dimension
                    if (variable->m_ShapeID == ShapeID::GlobalArray)
                    {
                        ndigits_dims[k] = ndigits(variable->m_Shape[k] - 1);
                    }
                    else
                    {
                        ndigits_dims[k] = ndigits(variable->m_Count[k] - 1);
                    }
                }

                size_t stepAbsolute = minBlocksInfo->Step;

                fprintf(outf, "%c       step %*zu: ", commentchar, ndigits_nsteps, stepAbsolute);
                fprintf(outf, "\n");
                const size_t blocksSize = blocks.size();
                ndigits_nblocks = ndigits(blocksSize - 1);

                for (size_t j = 0; j < blocksSize; j++)
                {
                    fprintf(outf, "%c         block %*zu: [", commentchar, ndigits_nblocks, j);

                    // just in case ndim for a block changes in LocalArrays:
                    ndim = variable->m_Count.size();

                    for (size_t k = 0; k < ndim; k++)
                    {
                        size_t c = (minBlocksInfo->WasLocalValue
                                        ? reinterpret_cast<size_t>(blocks[j].Count)
                                        : blocks[j].Count[k]);

                        if (c)
                        {
                            if (variable->m_ShapeID == ShapeID::GlobalArray)
                            {
                                size_t s = (minBlocksInfo->WasLocalValue
                                                ? reinterpret_cast<size_t>(blocks[j].Start)
                                                : blocks[j].Start[k]);
                                fprintf(outf, "%*zu:%*zu", ndigits_dims[k], s, ndigits_dims[k],
                                        s + c - 1);
                            }
                            else
                            {
                                // blockStart is empty vector for LocalArrays
                                fprintf(outf, "0:%*zu", ndigits_dims[k], c - 1);
                            }
                        }
                        else
                        {
                            fprintf(outf, "%-*s", 2 * ndigits_dims[k] + 1, "null");
                        }
                        if (k < ndim - 1)
                            fprintf(outf, ", ");
                    }
                    fprintf(outf, "]");

                    /* Print per-block statistics if available */
                    if (longopt)
                    {
                        if (true /* TODO: variable->has_minmax */)
                        {
                            fprintf(outf, " = ");
                            print_data(&blocks[j].MinMax.MinUnion, 0, adiosvartype, false);

                            fprintf(outf, " / ");
                            print_data(&blocks[j].MinMax.MaxUnion, 0, adiosvartype, false);
                        }
                        else
                        {
                            fprintf(outf, "N/A / N/A");
                        }
                    }
                    fprintf(outf, "\n");
                    if (dump)
                    {
                        Dims s, c;
                        if (variable->m_ShapeID == ShapeID::GlobalArray)
                        {
                            if (minBlocksInfo->WasLocalValue)
                            {
                                c = {reinterpret_cast<size_t>(blocks[j].Count)};
                                s = {reinterpret_cast<size_t>(blocks[j].Start)};
                            }
                            else
                            {
                                c = Dims(blocks[j].Count, blocks[j].Count + ndim);
                                s = Dims(blocks[j].Start, blocks[j].Start + ndim);
                            }
                        }
                        else
                        {
                            c = Dims(blocks[j].Count, blocks[j].Count + ndim);
                            s = Dims(0, ndim);
                        }
                        readVarBlock(fp, io, variable, RelStep, j, c, s);
                    }
                }
                delete minBlocksInfo;
            }
        }
        return;
    }

    std::map<size_t, std::vector<typename core::Variable<T>::BPInfo>> allblocks =
        fp->AllStepsBlocksInfo(*variable);
    if (allblocks.empty())
    {
        return;
    }
    size_t laststep = allblocks.rbegin()->first;
    int ndigits_nsteps = ndigits(laststep);
    if (variable->m_ShapeID == ShapeID::GlobalValue || variable->m_ShapeID == ShapeID::LocalValue)
    {
        // scalars
        for (auto &blockpair : allblocks)
        {
            size_t step = blockpair.first;
            std::vector<typename adios2::core::Variable<T>::BPInfo> &blocks = blockpair.second;
            fprintf(outf, "%c       step %*zu: ", commentchar, ndigits_nsteps, step);
            if (blocks.size() == 1)
            {
                fprintf(outf, " = ");
                print_data(&blocks[0].Value, 0, adiosvartype, true);
                fprintf(outf, "\n");
            }
            else
            {
                fprintf(outf, "%zu instances available\n", blocks.size());
            }
            if (dump)
            {
                fprintf(outf, "               ");
                int col = 0;
                for (size_t j = 0; j < blocks.size(); j++)
                {
                    print_data(&blocks[j].Value, 0, adiosvartype, true);
                    ++col;
                    if (j < blocks.size() - 1)
                    {
                        if (col < ncols)
                        {
                            fprintf(outf, " ");
                        }
                        else
                        {
                            fprintf(outf, "\n               ");
                            col = 0;
                        }
                    }
                }
                fprintf(outf, "\n");
            }
        }
        return;
    }
    else
    {
        // arrays
        size_t ndim = variable->m_Count.size();
        int ndigits_nblocks;
        for (size_t k = 0; k < ndim; k++)
        {
            // get digit lengths for each dimension
            if (variable->m_ShapeID == ShapeID::GlobalArray)
            {
                ndigits_dims[k] = ndigits(variable->m_Shape[k] - 1);
            }
            else
            {
                ndigits_dims[k] = ndigits(variable->m_Count[k] - 1);
            }
        }

        size_t stepRelative = 0;
        for (auto &blockpair : allblocks)
        {
            size_t stepAbsolute = blockpair.first;
            std::vector<typename adios2::core::Variable<T>::BPInfo> &blocks = blockpair.second;
            const size_t blocksSize = blocks.size();
            fprintf(outf, "%c       step %*zu: ", commentchar, ndigits_nsteps, stepAbsolute);
            fprintf(outf, "\n");
            ndigits_nblocks = ndigits(blocksSize - 1);

            for (size_t j = 0; j < blocksSize; j++)
            {
                fprintf(outf, "%c         block %*zu: [", commentchar, ndigits_nblocks, j);

                // just in case ndim for a block changes in LocalArrays:
                ndim = variable->m_Count.size();

                for (size_t k = 0; k < ndim; k++)
                {
                    if (blocks[j].Count[k])
                    {
                        if (variable->m_ShapeID == ShapeID::GlobalArray)
                        {
                            fprintf(outf, "%*zu:%*zu", ndigits_dims[k], blocks[j].Start[k],
                                    ndigits_dims[k], blocks[j].Start[k] + blocks[j].Count[k] - 1);
                        }
                        else
                        {
                            // blockStart is empty vector for LocalArrays
                            fprintf(outf, "0:%*zu", ndigits_dims[k], blocks[j].Count[k] - 1);
                        }
                    }
                    else
                    {
                        fprintf(outf, "%-*s", 2 * ndigits_dims[k] + 1, "null");
                    }
                    if (k < ndim - 1)
                        fprintf(outf, ", ");
                }
                fprintf(outf, "]");

                /* Print per-block statistics if available */
                if (longopt)
                {
                    if (true /* TODO: variable->has_minmax */)
                    {
                        fprintf(outf, " = ");
                        print_data(&blocks[j].Min, 0, adiosvartype, false);

                        fprintf(outf, " / ");
                        print_data(&blocks[j].Max, 0, adiosvartype, false);
                    }
                    else
                    {
                        fprintf(outf, "N/A / N/A");
                    }
                }
                fprintf(outf, "\n");
                if (dump)
                {
                    readVarBlock(fp, io, variable, stepRelative, j, blocks[j].Count,
                                 blocks[j].Start);
                }
            }
            ++stepRelative;
        }
    }
}

template <class T>
void print_decomp_singlestep(core::Engine *fp, core::IO *io, core::Variable<T> *variable)
{
    /* Print block info */
    DataType adiosvartype = variable->m_Type;
    const auto minBlocks = fp->MinBlocksInfo(*variable, fp->CurrentStep());

    std::vector<typename core::Variable<T>::BPInfo> coreBlocks;

    if (!minBlocks)
    {
        coreBlocks = fp->BlocksInfo(*variable, fp->CurrentStep());
    }
    if (!minBlocks && coreBlocks.empty())
    {
        return;
    }

    size_t blocksSize;
    if (minBlocks)
    {
        blocksSize = minBlocks->BlocksInfo.size();
    }
    else
    {
        blocksSize = coreBlocks.size();
    }
    const int ndigits_nblocks = ndigits(blocksSize - 1);

    if (variable->m_ShapeID == ShapeID::GlobalValue || variable->m_ShapeID == ShapeID::LocalValue)
    {
        // scalars
        if (dump)
        {
            int col = 0;
            int maxcols = ncols;
            if (adiosvartype == DataType::String)
            {
                maxcols = 1;
            }
            for (size_t j = 0; j < blocksSize; j++)
            {
                if (col == 0 && !noindex)
                {
                    fprintf(outf, "    (%*zu)    ", ndigits_nblocks, j);
                }
                if (!minBlocks)
                    print_data(&coreBlocks[j].Value, 0, adiosvartype, true);
                else
                    print_data(&minBlocks->BlocksInfo[j].BufferP, 0, adiosvartype, true);
                ++col;
                if (j < blocksSize - 1)
                {
                    if (col < maxcols)
                    {
                        fprintf(outf, " ");
                    }
                    else
                    {
                        col = 0;
                        fprintf(outf, "\n");
                    }
                }
            }
            fprintf(outf, "\n");
        }
        return;
    }
    else
    {
        // arrays
        size_t ndim = variable->m_Count.size();
        int ndigits_dims[32] = {
            8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
            8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        };
        for (size_t k = 0; k < ndim; k++)
        {
            // get digit lengths for each dimension
            if (variable->m_ShapeID == ShapeID::GlobalArray)
            {
                ndigits_dims[k] = ndigits(variable->m_Shape[k] - 1);
            }
            else
            {
                ndigits_dims[k] = ndigits(variable->m_Count[k] - 1);
            }
        }

        size_t stepRelative = 0;

        for (size_t j = 0; j < blocksSize; j++)
        {
            fprintf(outf, "%c         block %*zu: [", commentchar, ndigits_nblocks, j);

            // just in case ndim for a block changes in LocalArrays:
            ndim = variable->m_Count.size();

            for (size_t k = 0; k < ndim; k++)
            {
                if (!minBlocks)
                {
                    if (coreBlocks[j].Count[k])
                    {
                        if (variable->m_ShapeID == ShapeID::GlobalArray)
                        {
                            fprintf(outf, "%*zu:%*zu", ndigits_dims[k], coreBlocks[j].Start[k],
                                    ndigits_dims[k],
                                    coreBlocks[j].Start[k] + coreBlocks[j].Count[k] - 1);
                        }
                        else
                        {
                            // blockStart is empty vector for LocalArrays
                            fprintf(outf, "0:%*zu", ndigits_dims[k], coreBlocks[j].Count[k] - 1);
                        }
                    }
                    else
                    {
                        fprintf(outf, "%-*s", 2 * ndigits_dims[k] + 1, "null");
                    }
                }
                else
                {
                    if (minBlocks->BlocksInfo[j].Count[k])
                    {
                        if (variable->m_ShapeID == ShapeID::GlobalArray)
                        {
                            fprintf(outf, "%*zu:%*zu", ndigits_dims[k],
                                    minBlocks->BlocksInfo[j].Start[k], ndigits_dims[k],
                                    minBlocks->BlocksInfo[j].Start[k] +
                                        minBlocks->BlocksInfo[j].Count[k] - 1);
                        }
                        else
                        {
                            // blockStart is empty vector for LocalArrays
                            fprintf(outf, "0:%*zu", ndigits_dims[k],
                                    minBlocks->BlocksInfo[j].Count[k] - 1);
                        }
                    }
                    else
                    {
                        fprintf(outf, "%-*s", 2 * ndigits_dims[k] + 1, "null");
                    }
                }
                if (k < ndim - 1)
                    fprintf(outf, ", ");
            }
            fprintf(outf, "]");

            /* Print per-block statistics if available */
            if (longopt)
            {
                if (true /* TODO: variable->has_minmax */)
                {
                    if (!minBlocks)
                    {
                        fprintf(outf, " = ");
                        print_data(&coreBlocks[j].Min, 0, adiosvartype, false);

                        fprintf(outf, " / ");
                        print_data(&coreBlocks[j].Max, 0, adiosvartype, false);
                    }
                    else
                    {
                        fprintf(outf, " = ");
                        print_data(&minBlocks->BlocksInfo[j].MinMax.MinUnion, 0, adiosvartype,
                                   false);

                        fprintf(outf, " / ");
                        print_data(&minBlocks->BlocksInfo[j].MinMax.MaxUnion, 0, adiosvartype,
                                   false);
                    }
                }
                else
                {
                    fprintf(outf, "N/A / N/A");
                }
            }
            fprintf(outf, "\n");
            if (dump)
            {
                if (!minBlocks)
                {
                    readVarBlock(fp, io, variable, stepRelative, j, coreBlocks[j].Count,
                                 coreBlocks[j].Start);
                }
                else
                {
                    // GSE todo
                }
            }
        }
        ++stepRelative;
    }
    if (minBlocks)
        delete minBlocks;
}

int parseAccuracy()
{
    if (accuracy_def.empty())
        return 0;

    // Vector of string to save tokens
    std::vector<std::string> tokens;

    // stringstream class check1
    std::stringstream ss(accuracy_def);

    std::string intermediate;

    // Tokenizing w.r.t. space ','
    while (getline(ss, intermediate, ','))
    {
        tokens.push_back(intermediate);
    }

    if (tokens.size() != 3)
    {
        std::cout
            << "ERROR: --error definition needs 3 values: error, norm, and relative|absolute\n"
            << " error >= 0.0, norm >= 0.0 or inf, third arg is \"rel\" or \"abs\"\n";
        return 1;
    }

    accuracy.error = adios2::helper::StringTo<double>(tokens[0], "error value for accuracy");
    std::string normval = adios2::helper::LowerCase(tokens[2]);
    if (normval == "inf")
        accuracy.norm = std::numeric_limits<double>::infinity();
    else
        accuracy.norm = adios2::helper::StringTo<double>(tokens[1], "norm value for accuracy");
    std::string relval = adios2::helper::LowerCase(tokens[2]);
    if (relval == "rel")
        accuracy.relative = true;
    else if (relval == "abs")
        accuracy.relative = false;
    else
    {
        std::cout << "ERROR: --error third value must be \"rel\" or \"abs\" but it was \""
                  << tokens[2] << "\"\n";
        return 1;
    }

    accuracyWasSet = true;
    if (verbose > 0)
    {
        std::cout << "Read accuracy is set to (error = " << accuracy.error
                  << ", norm = " << accuracy.norm << ", " << (accuracy.relative ? "rel" : "abs")
                  << ")\n";
    }
    return 0;
}

// parse a string "0, 3; 027" into an integer array
// of [0,3,27]
// exits if parsing failes
void parseDimSpec(const std::string &str, int64_t *dims)
{
    if (str.empty())
        return;

    char *token;
    char *s; // copy of s; strtok modifies the string
    int i = 0;

    s = mystrndup(str.c_str(), 1024);
    // token = strtok_r(s, " ,;x\t\n", &saveptr);
    token = strtok(s, " ,;x\t\n");
    while (token != NULL && i < MAX_DIMS)
    {
        // printf("\t|%s|", token);
        errno = 0;
        dims[i] = (int64_t)strtoll(token, (char **)NULL, 0);
        if (errno)
        {
            fprintf(stderr,
                    "Error: could not convert field into a value: "
                    "%s from \"%s\"\n",
                    token, str.c_str());
            exit(200);
        }

        // get next item
        // token = strtok_r(NULL, " ,;x\t\n", &saveptr);
        token = strtok(NULL, " ,;x\t\n");
        i++;
    }
    // if (i>0) printf("\n");

    if (i > ndimsspecified)
        ndimsspecified = i;

    // check if number of dims specified is larger than we can handle
    if (token != NULL)
    {
        fprintf(stderr,
                "Error: More dimensions specified in \"%s\" than we "
                "can handle (%d)\n",
                str.c_str(), MAX_DIMS);
        exit(200);
    }
    free(s);
}

int compile_regexp_masks(void)
{
#ifdef USE_C_REGEX
    int errcode;
    char buf[256];
    for (int i = 0; i < nmasks; i++)
    {
        errcode = regcomp(&(varregex[i]), varmask[i], REG_EXTENDED);
        if (errcode)
        {
            regerror(errcode, &(varregex[i]), buf, sizeof(buf));
            fprintf(stderr,
                    "Error: \"%s\" is an invalid extended regular "
                    "expression: %s\n",
                    varmask[i], buf);
            return 2;
        }
    }
#else
    varregex.reserve(nmasks);
    for (int i = 0; i < nmasks; i++)
    {
        try
        {
            varregex.push_back(std::regex(varmask[i]));
        }
        catch (std::regex_error &e)
        {
            fprintf(stderr,
                    "Error: \"%s\" is an invalid extended regular "
                    "expression. C++ regex error code: %d\n",
                    varmask[i], e.code() == std::regex_constants::error_badrepeat);
            return 2;
        }
    }
#endif
    return 0;
}

char *mystrndup(const char *s, size_t n)
{
    char *t = nullptr;
    if (n > 0)
    {
        size_t slen = strlen(s);
        size_t len = (slen <= n ? slen : n);
        t = (char *)malloc(len + 1);
        if (t)
        {
            memcpy(t, s, len);
            t[len] = '\0';
        }
    }
    return t;
}

// end of namespace
}
}

int main(int argc, char *argv[])
{
    int retval = 1;
    try
    {
        retval = adios2::utils::bplsMain(argc, argv);
    }
    catch (std::exception &e)
    {
        std::cout << "\nbpls caught an exception\n";
        std::cout << e.what() << std::endl;
    }
    return retval;
}

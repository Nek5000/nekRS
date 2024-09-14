/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Reorganize.h
 *
 *  Created on: Mar 7, 2018
 *      Author: Norbert Podhorszki, pnorbert@ornl.gov
 */

#ifndef UTILS_REORGANIZE_REORGANIZE_H_
#define UTILS_REORGANIZE_REORGANIZE_H_

#include "adios2/core/IO.h"
#include "adios2/helper/adiosComm.h"
#include "utils/Utils.h"

namespace adios2
{
namespace utils
{

struct VarInfo
{
    core::VariableBase *v = nullptr;
    std::string type;
    Dims start;
    Dims count;
    size_t writesize = 0;    // size of subset this process writes, 0: do not write
    void *readbuf = nullptr; // read in buffer
};

class Reorganize : public Utils
{
public:
    Reorganize(int argc, char *argv[]);

    ~Reorganize() = default;

    void Run() final;

private:
    static const int m_CommSplitColor = 23731; // color in Comm::Split() call
    static const std::string m_HelpMessage;
    static const Params m_Options;

    std::string m_FileName;

    void ParseArguments() final;
    void ProcessParameters() final;
    void PrintUsage() const noexcept final;
    void PrintExamples() const noexcept final;
    void SetParameters(const std::string argument, const bool isLong) final;

    void CleanUpStep(core::IO &io);

    template <typename T>
    std::string VectorToString(const T &v);

    size_t Decompose(int numproc, int rank, VarInfo &vi,
                     const int *np // number of processes in each dimension
    );
    int ProcessMetadata(core::Engine &rStream, core::IO &io, const core::VarMap &variables,
                        const core::AttrMap &attributes, int step);
    int ReadWrite(core::Engine &rStream, core::Engine &wStream, core::IO &io,
                  const core::VarMap &variables, int step);
    Params parseParams(const std::string &param_str);

    // Input arguments
    std::string infilename;       // File/stream to read
    std::string outfilename;      // File to write
    std::string wmethodname;      // ADIOS write method
    std::string wmethodparam_str; // ADIOS write method parameter string
    std::string rmethodname;      // ADIOS read method
    std::string rmethodparam_str; // ADIOS read method parameter string

    static const size_t max_read_buffer_size = 16 * 1024 * 1024 * 1024ULL;
    static const size_t max_write_buffer_size = 16 * 1024 * 1024 * 1024ULL;

    // will stop if no data found for this time (-1: never stop)
    static const int timeout_sec = 300;

    // Global variables
    int m_Rank = 0;
    int m_Size = 1;
    helper::Comm m_Comm;

    // Read/write method parameters
    Params rmethodparams;
    Params wmethodparams;

    uint64_t write_total = 0;   // data size read/written by one processor
    uint64_t largest_block = 0; // the largest variable block one process reads

    // Timeout handling:
    // BPFile, BP3, File, "" (default), HDF5 is handled as file, stop at timeout
    // BP4, FileStream and everything else is handled as stream, wait forever
    // for stream termination
    bool handleAsStream = true;

    int decomp_values[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    template <typename Arg, typename... Args>
    void print0(Arg &&arg, Args &&...args);

    template <typename Arg, typename... Args>
    void osprint0(std::ostream &out, Arg &&arg, Args &&...args);
};

} // end namespace utils
} // end namespace adios2

#endif /* UTILS_REORGANIZE_REORGANIZE_H_ */

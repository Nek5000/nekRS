/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Reorganize.cpp
 *
 *  Created on: Mar 7, 2018
 *      Author: Norbert Podhorszki, pnorbert@ornl.gov
 *
 * Reorganize global arrays
   Assumptions:
     - one output step fits into the memory of the reorganizer.
       Actually, this means, even more memory is needed than the size of output.
       We need to read each variable while also buffering all of them for
 output.
     - output steps contain the same variable set (no changes in variables)
     - attributes are the same for all steps (will write only once here)
 */

#include "Reorganize.h"

#include <assert.h>
#include <iomanip>
#include <string>

#include "adios2/common/ADIOSMacros.h"
#include "adios2/core/ADIOS.h"
#include "adios2/core/Engine.h"
#include "adios2/core/IO.h"
#include "adios2/helper/adiosComm.h"
#include "adios2/helper/adiosFunctions.h"
#include "adios2/helper/adiosLog.h"
#include "adios2/helper/adiosString.h"

#if ADIOS2_USE_MPI
#include "adios2/helper/adiosCommMPI.h"
#else
#include "adios2/helper/adiosCommDummy.h"
#endif

// C headers
#include <cerrno>
#include <cstdlib>

namespace adios2
{
namespace utils
{

Reorganize::Reorganize(int argc, char *argv[]) : Utils("adios_reorganize", argc, argv)
{
#if ADIOS2_USE_MPI
    {
        auto commWorld = helper::CommWithMPI(MPI_COMM_WORLD);
        m_Comm = commWorld.Split(m_CommSplitColor, 0);
    }
#else
    m_Comm = helper::CommDummy();
#endif
    m_Rank = m_Comm.Rank();
    m_Size = m_Comm.Size();

    if (argc < 7)
    {
        PrintUsage();
        helper::Throw<std::invalid_argument>("Utils", "AdiosReorganize", "Reorganize",
                                             "Not enough arguments. At least 6 are required");
    }
    infilename = std::string(argv[1]);
    outfilename = std::string(argv[2]);
    rmethodname = std::string(argv[3]);
    rmethodparam_str = std::string(argv[4]);
    wmethodname = std::string(argv[5]);
    wmethodparam_str = std::string(argv[6]);

    int nd = 0;
    int j = 7;
    char *end;
    while (argc > j && j < 13)
    { // get max 6 dimensions
        errno = 0;
        decomp_values[nd] = std::strtol(argv[j], &end, 10);
        if (errno || (end != 0 && *end != '\0'))
        {
            std::string errmsg("ERROR: Invalid decomposition number in argument " +
                               std::to_string(j) + ": '" + std::string(argv[j]) + "'\n");
            PrintUsage();
            helper::Throw<std::invalid_argument>("Utils", "AdiosReorganize", "Reorganize", errmsg);
        }
        nd++;
        j++;
    }

    if (argc > j)
    {
        helper::Throw<std::invalid_argument>("Utils", "AdiosReorganize", "Reorganize",
                                             "Up to 6 decomposition arguments are supported");
    }

    int prod = 1;
    for (int i = 0; i < nd; i++)
    {
        prod *= decomp_values[i];
    }

    if (prod > m_Size)
    {
        print0("ERROR: Product of decomposition numbers %d > number of "
               "processes %d\n",
               prod, m_Size);
        std::string errmsg("ERROR: The product of decomposition numbers " + std::to_string(prod) +
                           " > number of processes " + std::to_string(m_Size) + "\n");
        PrintUsage();
        helper::Throw<std::invalid_argument>("Utils", "AdiosReorganize", "Reorganize", errmsg);
    }
}

void Reorganize::Run()
{
    ParseArguments();
    ProcessParameters();
    int retval = 0;

    print0("Input stream            = ", infilename);
    print0("Output stream           = ", outfilename);
    print0("Read method             = ", rmethodname);
    print0("Read method parameters  = ", rmethodparam_str);
    print0("Write method            = ", wmethodname);
    print0("Write method parameters = ", wmethodparam_str);

    core::ADIOS adios(m_Comm.Duplicate(), "C++");
    core::IO &io = adios.DeclareIO("group");

    print0("Waiting to open stream ", infilename, "...");

    io.SetEngine(rmethodname);
    io.SetParameters(rmethodparams);
    core::Engine &rStream = io.Open(infilename, adios2::Mode::Read);
    // rStream.FixedSchedule();

    io.ClearParameters();
    io.SetEngine(wmethodname);
    io.SetParameters(wmethodparams);
    core::Engine &wStream = io.Open(outfilename, adios2::Mode::Write);

    int steps = 0;
    int curr_step = -1;
    while (true)
    {
        adios2::StepStatus status = rStream.BeginStep(adios2::StepMode::Read, 10.0);
        if (status == adios2::StepStatus::NotReady)
        {
            if (handleAsStream)
            {
                if (!m_Rank)
                {
                    std::cout << " No new steps arrived in a while " << std::endl;
                }
                continue;
            }
            else
            {
                if (!m_Rank)
                {
                    std::cout << " Timeout waiting for next step. If this is "
                                 "a live stream through file, use a different "
                                 "reading engine, like FileStream or BP4. "
                                 "If it is an unclosed BP file, you may manually "
                                 "close it with using adios_deactive_bp.sh."
                              << std::endl;
                }
                break;
            }
        }
        else if (status != adios2::StepStatus::OK)
        {
            break;
        }

        steps++; // start counting from 1

        if (rStream.CurrentStep() != static_cast<size_t>(curr_step + 1))
        {
            // we missed some steps
            std::cout << "rank " << m_Rank << " WARNING: steps " << curr_step << ".."
                      << rStream.CurrentStep() - 1 << "were missed when advancing." << std::endl;
        }

        curr_step = static_cast<int>(rStream.CurrentStep());
        const core::VarMap &variables = io.GetVariables();
        const core::AttrMap &attributes = io.GetAttributes();

        print0("____________________\n\nFile info:");
        print0("  current step:   ", curr_step);
        print0("  # of variables: ", variables.size());
        print0("  # of attributes: ", attributes.size());

        retval = ProcessMetadata(rStream, io, variables, attributes, steps);
        if (retval)
            break;

        retval = ReadWrite(rStream, wStream, io, variables, steps);
        if (retval)
            break;

        CleanUpStep(io);
    }

    rStream.Close();
    wStream.Close();
    print0("Bye after processing ", steps, " steps");
}

// PRIVATE
template <typename Arg, typename... Args>
void Reorganize::osprint0(std::ostream &out, Arg &&arg, Args &&...args)
{
    if (!m_Rank)
    {
        out << std::forward<Arg>(arg);
        using expander = int[];
        (void)expander{0, (void(out << std::forward<Args>(args)), 0)...};
        std::cout << std::endl;
    }
}

template <typename Arg, typename... Args>
void Reorganize::print0(Arg &&arg, Args &&...args)
{
    if (!m_Rank)
    {
        std::cout << std::forward<Arg>(arg);
        using expander = int[];
        (void)expander{0, (void(std::cout << std::forward<Args>(args)), 0)...};
        std::cout << std::endl;
    }
}

Params Reorganize::parseParams(const std::string &param_str)
{
    std::istringstream ss(param_str);
    std::vector<std::string> kvs;
    std::string kv;

    while (std::getline(ss, kv, ','))
    {
        kvs.push_back(kv);
    }

    return helper::BuildParametersMap(kvs, '=');
}

void Reorganize::ParseArguments()
{
    rmethodparams = parseParams(rmethodparam_str);
    wmethodparams = parseParams(wmethodparam_str);
}

void Reorganize::ProcessParameters()
{
    if (rmethodname.empty())
    {
        handleAsStream = false;
        return;
    }
    std::string s(rmethodname);
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    if (s == "file" || s == "bpfile" || s == "bp3" || s == "hdf5")
    {
        handleAsStream = false;
    }
    else
    {
        handleAsStream = true;
    }
}

void Reorganize::PrintUsage() const noexcept
{
    std::cout << "Usage: adios_reorganize input output rmethod \"params\" wmethod "
                 "\"params\" "
                 "<decomposition>\n"
                 "    input   Input stream path\n"
                 "    output  Output file path\n"
                 "    rmethod ADIOS method to read with\n"
                 "            Supported read methods: BPFile, HDF5, SST, SSC, "
                 "DataMan\n"
                 "    params  Read method parameters (in quotes; comma-separated "
                 "list)\n"
                 "    wmethod ADIOS method to write with\n"
                 "    params  Write method parameters (in quotes; comma-separated "
                 "list)\n"
                 "    <decomposition>    list of numbers e.g. 32 8 4\n"
                 "            Decomposition values in each dimension of an array\n"
                 "            The product of these number must be less then the "
                 "number\n"
                 "            of processes. Processes whose rank is higher than the\n"
                 "            product, will not write anything.\n"
                 "               Arrays with less dimensions than the number of "
                 "values,\n"
                 "            will be decomposed with using the appropriate number "
                 "of\n"
                 "            values."
              << std::endl;
}

void Reorganize::PrintExamples() const noexcept {}

void Reorganize::SetParameters(const std::string argument, const bool isLong) {}

std::vector<VarInfo> varinfo;

// cleanup all info from previous step except
// do
//   remove all variable and attribute definitions from output group
//   free all varinfo (will be inquired again at next step)
//   free read buffer (required size may change at next step)
// do NOT
//   destroy group
//
void Reorganize::CleanUpStep(core::IO &io)
{
    for (auto &vi : varinfo)
    {
        if (vi.readbuf != nullptr)
        {
            free(vi.readbuf);
        }
    }
    varinfo.clear();
    // io.RemoveAllVariables();
    // io.RemoveAllAttributes();
}

template <typename T>
std::string Reorganize::VectorToString(const T &v)
{
    std::string s;
    for (const auto e : v)
    {
        s += std::to_string(e) + ", ";
    }
    s.pop_back();
    s.pop_back();
    return s;
}

size_t Reorganize::Decompose(int numproc, int rank, VarInfo &vi,
                             const int *np // number of processes in each dimension
)
{
    size_t writesize = 0;
    if (vi.v == nullptr)
    {
        return writesize;
    }

    /* Handle local array: only one block for now */
    if (vi.v->m_ShapeID == adios2::ShapeID::LocalArray)
    {
        if (rank == 0)
        {
            writesize = 1;
            for (size_t i = 0; i < vi.v->m_Count.size(); i++)
            {
                writesize *= vi.v->m_Count[i];
                // vi.start.push_back(0);
                vi.count.push_back(vi.v->m_Count[i]);
            }
        }
        else
        {
            writesize = 0;
        }
        return writesize;
    }

    size_t ndim = vi.v->Shape().size();

    /* Scalars */
    if (ndim == 0)
    {
        // scalars -> rank 0 writes them
        if (rank == 0)
            writesize = 1;
        else
            writesize = 0;
        return writesize;
    }

    /* Global Arrays */
    /* calculate this process' position in the n-dim space
    0 1 2
    3 4 5
    6 7 8

    for 1D:
    posx = rank/1             ! 1st dim: 0, 1, 2...,rank-1 are in the same X
    position

    for 2D:
    posx = mod(rank, npx)     ! 1st dim: 0, npx, 2npx... are in the same X
    position
    posy = rank/(npx)         ! 2nd dim: npx processes belong into one dim

    for 3D:
    posx = mod(rank, npx)     ! 1st dim: 0, npx, 2npx... are in the same X
    position
    posy = mod(rank/npx, npy) ! 2nd dim: (0, npx-1) have the same dim (so divide
    with npx first)
    posz = rank/(npx*npy)     ! 3rd dim: npx*npy processes belong into one dim
    */
    int nps = 1;
    std::vector<int> pos(ndim); // rank's position in each dimensions
    vi.start.reserve(ndim);
    vi.count.reserve(ndim);

    size_t i = 0;
    for (i = 0; i < ndim - 1; i++)
    {
        pos[i] = (rank / nps) % np[i];
        nps *= np[i];
    }
    pos[i] = rank / nps;

    std::string ints = VectorToString(pos);
    if (pos[ndim - 1] >= np[ndim - 1])
    {
        std::cout << "rank " << rank << ": position in " << ndim << "-D decomposition = " << ints
                  << " ---> Out of bound process" << std::endl;
    }
    else
    {
        std::cout << "rank " << rank << ": position in " << ndim << "-D decomposition = " << ints
                  << std::endl;
    }

    /* Decompose each dimension according to the position */
    writesize = 1;
    for (i = 0; i < ndim; i++)
    {
        size_t start, count;
        if (pos[ndim - 1] >= np[ndim - 1])
        {
            // this process gets nothing to read
            count = 0;
            start = 0;
        }
        else
        {
            count = vi.v->Shape()[i] / np[i];
            start = count * pos[i];
            if (pos[i] == np[i] - 1)
            {
                // last one in the dimension may need to read more than the rest
                count = vi.v->Shape()[i] - count * (np[i] - 1);
            }
        }
        vi.start.push_back(start);
        vi.count.push_back(count);
        writesize *= count;
    }
    ints = VectorToString(vi.count);
    std::cout << "rank " << rank << ": ldims in " << ndim << "-D space = {" << ints << "}"
              << std::endl;
    ints = VectorToString(vi.start);
    std::cout << "rank " << rank << ": offsets in " << ndim << "-D space = {" << ints << "}"
              << std::endl;
    return writesize;
}

int Reorganize::ProcessMetadata(core::Engine &rStream, core::IO &io, const core::VarMap &variables,
                                const core::AttrMap &attributes, int step)
{
    int retval = 0;

    varinfo.resize(variables.size());
    write_total = 0;
    largest_block = 0;

    // Decompose each variable and calculate output buffer size
    int varidx = 0;
    for (const auto &variablePair : variables)
    {
        const std::string &name(variablePair.first);
        const DataType type(variablePair.second->m_Type);
        core::VariableBase *variable = nullptr;
        print0("Get info on variable ", varidx, ": ", name);
        size_t nBlocks = 1;

        if (type == DataType::Struct)
        {
            // not supported
        }
#define declare_template_instantiation(T)                                                          \
    else if (type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        core::Variable<T> *v = io.InquireVariable<T>(variablePair.first);                          \
        if (v->m_ShapeID == adios2::ShapeID::LocalArray)                                           \
        {                                                                                          \
                                                                                                   \
            const auto minBlocks = rStream.MinBlocksInfo(*v, step);                                \
            if (minBlocks)                                                                         \
            {                                                                                      \
                nBlocks = minBlocks->BlocksInfo.size();                                            \
            }                                                                                      \
            else                                                                                   \
            {                                                                                      \
                auto blocks = rStream.BlocksInfo(*v, rStream.CurrentStep());                       \
                nBlocks = blocks.size();                                                           \
            }                                                                                      \
        }                                                                                          \
        variable = v;                                                                              \
    }
        ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

        varinfo[varidx].v = variable;

        if (variable != nullptr)
        {

            // print variable type and dimensions
            if (!m_Rank)
            {
                std::cout << "    " << ToString(type) << " " << name;
            }
            // if (variable->Shape().size() > 0)
            if (variable->m_ShapeID == adios2::ShapeID::GlobalArray)
            {
                if (!m_Rank)
                {
                    std::cout << "[" << variable->Shape()[0];
                    for (size_t j = 1; j < variable->Shape().size(); j++)
                    {
                        std::cout << ", " << variable->Shape()[j];
                    }
                    std::cout << "]" << std::endl;
                }
            }

            else if (variable->m_ShapeID == adios2::ShapeID::GlobalValue)
            {
                print0("\tscalar");
            }
            else if (variable->m_ShapeID == adios2::ShapeID::LocalArray)
            {
                print0("\t local array ");
                if (nBlocks > 1)
                {
                    print0("ERROR: adios_reorganize does not support Local Arrays "
                           "except when there is only 1 written block in each "
                           "step. This one has ",
                           nBlocks, " blocks in this step ");
                    return 1;
                }
            }
            else
            {
                print0("\n *** Unidentified object ", name, " ***\n");
                return 1;
            }

            // determine subset we will write
            size_t sum_count = Decompose(m_Size, m_Rank, varinfo[varidx], decomp_values);
            varinfo[varidx].writesize = sum_count * variable->m_ElementSize;

            if (varinfo[varidx].writesize != 0)
            {
                write_total += varinfo[varidx].writesize;
                if (largest_block < varinfo[varidx].writesize)
                    largest_block = varinfo[varidx].writesize;
            }
        }
        else
        {
            print0("    Not available in this step");
        }
        ++varidx;
    }

    // determine output buffer size
    size_t bufsize = write_total + variables.size() * 200 + attributes.size() * 32 + 1024;
    if (bufsize > max_write_buffer_size)
    {
        helper::Log("Util", "Reorganize", "ProcessMetadata",
                    "write buffer size needs to hold about " + std::to_string(bufsize) +
                        " bytes but max is set to " + std::to_string(max_write_buffer_size),
                    m_Rank, m_Rank, 0, 0, helper::FATALERROR);
        return 1;
    }

    if (largest_block > max_read_buffer_size)
    {
        helper::Log("Util", "Reorganize", "ProcessMetadata",
                    "read buffer size needs to hold at least " + std::to_string(largest_block) +
                        " bytes but max is set to " + std::to_string(max_read_buffer_size),
                    m_Rank, m_Rank, 0, 0, helper::FATALERROR);
        return 1;
    }
    return retval;
}

int Reorganize::ReadWrite(core::Engine &rStream, core::Engine &wStream, core::IO &io,
                          const core::VarMap &variables, int step)
{
    int retval = 0;

    size_t nvars = variables.size();
    if (nvars != varinfo.size())
    {
        helper::Log("Util", "Reorganize", "ReadWrite",
                    "Invalid program state, number of variables (" + std::to_string(nvars) +
                        ") to read does not match the number of processed variables (" +
                        std::to_string(varinfo.size()) + ")",
                    m_Rank, m_Rank, 0, 0, helper::FATALERROR);
    }

    /*
     * Read all variables into memory
     */
    for (size_t varidx = 0; varidx < nvars; ++varidx)
    {
        if (varinfo[varidx].v != nullptr)
        {
            const std::string &name = varinfo[varidx].v->m_Name;
            assert(varinfo[varidx].readbuf == nullptr);
            if (varinfo[varidx].writesize != 0)
            {
                // read variable subset
                std::cout << "rank " << m_Rank << ": Read variable " << name << std::endl;
                const DataType type = variables.at(name)->m_Type;
                if (type == DataType::Struct)
                {
                    // not supported
                }
#define declare_template_instantiation(T)                                                          \
    else if (type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        varinfo[varidx].readbuf = calloc(1, varinfo[varidx].writesize);                            \
        if (varinfo[varidx].count.size() == 0)                                                     \
        {                                                                                          \
            rStream.Get<T>(name, reinterpret_cast<T *>(varinfo[varidx].readbuf),                   \
                           adios2::Mode::Sync);                                                    \
        }                                                                                          \
        else                                                                                       \
        {                                                                                          \
            varinfo[varidx].v->SetSelection({varinfo[varidx].start, varinfo[varidx].count});       \
            rStream.Get<T>(name, reinterpret_cast<T *>(varinfo[varidx].readbuf));                  \
        }                                                                                          \
    }
                ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation
            }
        }
    }
    rStream.EndStep(); // read in data into allocated pointers

    /*
     * Write all variables
     */
    wStream.BeginStep();
    for (size_t varidx = 0; varidx < nvars; ++varidx)
    {
        if (varinfo[varidx].v != nullptr)
        {
            const std::string &name = varinfo[varidx].v->m_Name;
            if (varinfo[varidx].writesize != 0)
            {
                // Write variable subset
                std::cout << "rank " << m_Rank << ": Write variable " << name << std::endl;
                const DataType type = variables.at(name)->m_Type;
                if (type == DataType::Struct)
                {
                    // not supported
                }
#define declare_template_instantiation(T)                                                          \
    else if (type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        if (varinfo[varidx].count.size() == 0)                                                     \
        {                                                                                          \
            wStream.Put<T>(name, reinterpret_cast<T *>(varinfo[varidx].readbuf),                   \
                           adios2::Mode::Sync);                                                    \
        }                                                                                          \
        else if (varinfo[varidx].v->m_ShapeID == adios2::ShapeID::LocalArray)                      \
        {                                                                                          \
            wStream.Put<T>(name, reinterpret_cast<T *>(varinfo[varidx].readbuf),                   \
                           adios2::Mode::Sync);                                                    \
        }                                                                                          \
        else                                                                                       \
        {                                                                                          \
            varinfo[varidx].v->SetSelection({varinfo[varidx].start, varinfo[varidx].count});       \
            wStream.Put<T>(name, reinterpret_cast<T *>(varinfo[varidx].readbuf));                  \
        }                                                                                          \
    }
                ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation
            }
        }
    }
    wStream.EndStep(); // write output buffer to file
    return retval;
}

} // end namespace utils
} // end namespace adios2

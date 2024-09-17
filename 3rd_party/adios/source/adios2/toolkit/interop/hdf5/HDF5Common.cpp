/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * HDF5Common.cpp
 *
 *  Created on: April 20, 2017
 *      Author: Junmin
 */

#include "HDF5Common.h"
#include "HDF5Common.tcc"

#include <complex>
#include <ios>
#include <iostream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <vector>

#include "adios2/helper/adiosFunctions.h" // IsRowMajor
#include <cstring>                        // strlen

namespace adios2
{
namespace interop
{

std::mutex HDF5Common_MPI_API_Mutex;
HDF5Common::MPI_API const *HDF5Common_MPI_API;

namespace
{

HDF5Common::MPI_API const *GetHDF5Common_MPI_API()
{
    std::lock_guard<std::mutex> guard(HDF5Common_MPI_API_Mutex);
    (void)guard; // Workaround to silence compiler warning about unused variable
    return HDF5Common_MPI_API;
}

} // end anonymous namespace

const std::string HDF5Common::ATTRNAME_NUM_STEPS = "NumSteps";
const std::string HDF5Common::ATTRNAME_GIVEN_ADIOSNAME = "ADIOSName";
const std::string HDF5Common::PREFIX_BLOCKINFO = "ADIOS_BLOCKINFO_";
const std::string HDF5Common::PREFIX_STAT = "ADIOS_STAT_";
const std::string HDF5Common::PARAMETER_COLLECTIVE = "H5CollectiveMPIO";
const std::string HDF5Common::PARAMETER_CHUNK_FLAG = "H5ChunkDim";
const std::string HDF5Common::PARAMETER_CHUNK_VARS = "H5ChunkVars";
const std::string HDF5Common::PARAMETER_HAS_IDLE_WRITER_RANK = "IdleH5Writer";

#define CHECK_H5_RETURN(returnCode, reason)                                                        \
    {                                                                                              \
        if (returnCode < 0)                                                                        \
            helper::Throw<std::runtime_error>("Toolkit", "interop::hdf5::HDF5Common",              \
                                              "CHECK_H5_RETURN", reason);                          \
    }
/*
   //need to know ndim before defining this.
   //inconvenient
class HDF5BlockStat
{
public:
  size_t m_offset[ndim];
  size_t m_Count[ndim];
  double m_Min;
  double m_Max;
};
*/

HDF5Common::HDF5Common()
{
    m_DefH5TypeComplexFloat = H5Tcreate(H5T_COMPOUND, sizeof(std::complex<float>));
    H5Tinsert(m_DefH5TypeComplexFloat, "r", 0, H5T_NATIVE_FLOAT);
    H5Tinsert(m_DefH5TypeComplexFloat, "i", H5Tget_size(H5T_NATIVE_FLOAT), H5T_NATIVE_FLOAT);

    m_DefH5TypeComplexDouble = H5Tcreate(H5T_COMPOUND, sizeof(std::complex<double>));
    H5Tinsert(m_DefH5TypeComplexDouble, "r", 0, H5T_NATIVE_DOUBLE);
    H5Tinsert(m_DefH5TypeComplexDouble, "i", H5Tget_size(H5T_NATIVE_DOUBLE), H5T_NATIVE_DOUBLE);

    m_DefH5TypeComplexLongDouble = H5Tcreate(H5T_COMPOUND, sizeof(std::complex<long double>));
    H5Tinsert(m_DefH5TypeComplexLongDouble, "r", 0, H5T_NATIVE_LDOUBLE);
    H5Tinsert(m_DefH5TypeComplexLongDouble, "i", H5Tget_size(H5T_NATIVE_LDOUBLE),
              H5T_NATIVE_LDOUBLE);

    m_PropertyTxfID = H5Pcreate(H5P_DATASET_XFER);

    size_t size_t_max = (size_t)-1;
    if (size_t_max <= std::numeric_limits<unsigned int>::max())
        m_TimeStepH5T = H5T_NATIVE_UINT;
    else if (size_t_max <= std::numeric_limits<unsigned long>::max())
        m_TimeStepH5T = H5T_NATIVE_ULONG;

    // default is m_TimeStepH5T =  H5T_NATIVE_ULLONG;
}

HDF5Common::~HDF5Common() { Close(); }

void HDF5Common::ParseParameters(core::IO &io)
{
    if (m_MPI)
    {
        m_MPI->set_dxpl_mpio(m_PropertyTxfID,
                             H5FD_MPIO_INDEPENDENT); // explicit
        auto itKey = io.m_Parameters.find(PARAMETER_COLLECTIVE);
        if (itKey != io.m_Parameters.end())
        {
            if (itKey->second == "yes" || itKey->second == "true")
                m_MPI->set_dxpl_mpio(m_PropertyTxfID, H5FD_MPIO_COLLECTIVE);
        }

        itKey = io.m_Parameters.find(PARAMETER_HAS_IDLE_WRITER_RANK);
        if (itKey != io.m_Parameters.end())
        {
            if (itKey->second == "yes" || itKey->second == "true")
                m_IdleWriterOn = true;
        }
    }

    m_ChunkVarNames.clear();
    m_ChunkPID = -1;
    m_ChunkDim = 0;

    {
        std::vector<hsize_t> chunkDim;
        auto chunkFlagKey = io.m_Parameters.find(PARAMETER_CHUNK_FLAG);
        if (chunkFlagKey != io.m_Parameters.end())
        { // note space is the delimiter
            std::stringstream ss(chunkFlagKey->second);
            int i;
            while (ss >> i)
                chunkDim.push_back(i);

            m_ChunkPID = H5Pcreate(H5P_DATASET_CREATE);
            m_ChunkDim = chunkDim.size();

            size_t max_int = static_cast<size_t>(std::numeric_limits<int>::max());
            if (m_ChunkDim > max_int)
            {
                helper::Throw<std::overflow_error>("Toolkit", "interop::hdf5::HDF5Common",
                                                   "ParseParameters",
                                                   "chunkDim.size() is "
                                                   "too large to be represented by an int");
            }

            if (m_ChunkDim > 0)
                H5Pset_chunk(m_ChunkPID, static_cast<int>(m_ChunkDim), chunkDim.data());
        }
    }

    //
    // if no chunk dim specified, then ignore this parameter
    //
    if (-1 != m_ChunkPID)
    {
        auto chunkVarKey = io.m_Parameters.find(PARAMETER_CHUNK_VARS);
        if (chunkVarKey != io.m_Parameters.end())
        {
            std::stringstream ss(chunkVarKey->second);
            std::string token;
            while (ss >> token)
                m_ChunkVarNames.insert(token);
        }
    }

    m_OrderByC = (io.m_ArrayOrder == ArrayOrdering::RowMajor);
}

void HDF5Common::Append(const std::string &name, helper::Comm const &comm)
{
    m_PropertyListId = H5Pcreate(H5P_FILE_ACCESS);

    if (MPI_API const *mpi = GetHDF5Common_MPI_API())
    {
        if (mpi && mpi->init(comm, m_PropertyListId, &m_CommRank, &m_CommSize))
        {
            m_MPI = mpi;
        }
    }

    m_FileId = H5Fopen(name.c_str(), H5F_ACC_RDWR, m_PropertyListId);
    H5Pclose(m_PropertyListId);

    std::string ts0;
    StaticGetAdiosStepString(ts0, 0);

    if (m_FileId >= 0)
    {
        if (H5Lexists(m_FileId, ts0.c_str(), H5P_DEFAULT) != 0)
        {
            m_IsGeneratedByAdios = true;
        }
        if (!m_IsGeneratedByAdios)
            helper::Throw<std::ios_base::failure>("Toolkit", "interop::hdf5::HDF5Common", "Append",
                                                  "Likely no such file." + name);

        GetNumAdiosSteps(); // read how many steps exists in this file

        if (0 == m_NumAdiosSteps)
            helper::Throw<std::ios_base::failure>("Toolkit", "interop::hdf5::HDF5Common", "Append",
                                                  "No valid steps found in " + name);
        if (1 == m_NumAdiosSteps)
            m_GroupId = H5Gopen(m_FileId, ts0.c_str(), H5P_DEFAULT);
        else
            SetAdiosStep(m_NumAdiosSteps - 1);

        m_WriteMode = true;
        Advance();
    }
    else
        helper::Throw<std::ios_base::failure>("Toolkit", "interop::hdf5::HDF5Common", "Append",
                                              "Likely no such file." + name);
}

void HDF5Common::Init(const std::string &name, helper::Comm const &comm, bool toWrite)
{
    m_WriteMode = toWrite;
    m_PropertyListId = H5Pcreate(H5P_FILE_ACCESS);

    if (MPI_API const *mpi = GetHDF5Common_MPI_API())
    {
        if (mpi && mpi->init(comm, m_PropertyListId, &m_CommRank, &m_CommSize))
        {
            m_MPI = mpi;
        }
    }

    // std::string ts0 = "/AdiosStep0";
    std::string ts0;
    StaticGetAdiosStepString(ts0, 0);

#ifdef H5_HAVE_SUBFILING_VFD
    bool useMPI = false;
    const char *temp = getenv("H5FD_SUBFILING_IOC_PER_NODE");

    if (NULL != temp)
    {
        int itemp = -1;
        sscanf(temp, "%d", &itemp);
        if (0 == itemp)
            useMPI = true;
    }

    if (!useMPI)
        H5Pset_fapl_subfiling(m_PropertyListId, NULL);
#endif
    if (toWrite)
    {
        /*
         * Create a new file collectively and release property list identifier.
         */
        m_FileId = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, m_PropertyListId);
        if (m_FileId >= 0)
        {
            m_GroupId = H5Gcreate2(m_FileId, ts0.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            if (m_GroupId < 0)
            {
                helper::Throw<std::ios_base::failure>(
                    "Toolkit", "interop::hdf5::HDF5Common", "Init",
                    "Unable to create HDF5 group " + ts0 + " in call to Open");
            }
        }
    }
    else
    {
        // read a file collectively
        m_FileId = H5Fopen(name.c_str(), H5F_ACC_RDONLY, m_PropertyListId);
        if (m_FileId >= 0)
        {
            if (H5Lexists(m_FileId, ts0.c_str(), H5P_DEFAULT) != 0)
            {
                m_GroupId = H5Gopen(m_FileId, ts0.c_str(), H5P_DEFAULT);
                m_IsGeneratedByAdios = true;
            }
        }
    }

    H5Pclose(m_PropertyListId);
}

void HDF5Common::WriteAdiosSteps()
{
    if (m_FileId < 0)
    {
        helper::Throw<std::invalid_argument>("Toolkit", "interop::hdf5::HDF5Common",
                                             "WriteAdiosSteps",
                                             "invalid HDF5 file to record "
                                             "steps, in call to Write");
    }

    if (!m_WriteMode)
    {
        return;
    }

    hid_t s = H5Screate(H5S_SCALAR);
    hid_t attr = H5Aexists(m_FileId, ATTRNAME_NUM_STEPS.c_str());

    if (0 == attr)
        attr = H5Acreate(m_FileId, ATTRNAME_NUM_STEPS.c_str(), m_TimeStepH5T, s, H5P_DEFAULT,
                         H5P_DEFAULT);
    else
        attr = H5Aopen(m_FileId, ATTRNAME_NUM_STEPS.c_str(), H5P_DEFAULT);

    size_t totalAdiosSteps = m_CurrentAdiosStep + 1;

    if (m_GroupId < 0)
    {
        totalAdiosSteps = m_CurrentAdiosStep;
    }

    H5Awrite(attr, m_TimeStepH5T, &totalAdiosSteps);

    H5Sclose(s);
    H5Aclose(attr);
}

size_t HDF5Common::GetAdiosStep() const { return m_NumAdiosSteps; }

// A valid file should have at least 1 step.
size_t HDF5Common::GetNumAdiosSteps()
{
    if (m_WriteMode)
    {
        // return static_cast<unsigned int>(-1);
        return 0;
    }

    if (m_FileId < 0)
    {
        helper::Throw<std::invalid_argument>("Toolkit", "interop::hdf5::HDF5Common",
                                             "GetNumAdiosSteps",
                                             "invalid HDF5 file to read step attribute");
    }

    if (!m_IsGeneratedByAdios)
    {
        return 1;
    }

    if (m_NumAdiosSteps <= 0)
    {
        hsize_t numobj;
        H5Gget_num_objs(m_FileId, &numobj);
        m_NumAdiosSteps = numobj;

        if (H5Aexists(m_FileId, ATTRNAME_NUM_STEPS.c_str()))
        {
            hid_t attr = H5Aopen(m_FileId, ATTRNAME_NUM_STEPS.c_str(), H5P_DEFAULT);

            H5Aread(attr, m_TimeStepH5T, &m_NumAdiosSteps);
            H5Aclose(attr);
        }
    }

    return m_NumAdiosSteps;
}

// read from all time steps
void HDF5Common::ReadAllVariables(core::IO &io)
{
    if (!m_IsGeneratedByAdios)
    {
        FindVarsFromH5(io, m_FileId, "/", "", 0);
        return;
    }

    GetNumAdiosSteps();

    for (size_t i = 0; i < m_NumAdiosSteps; i++)
    {
        ReadVariables(i, io);
    }
}

void HDF5Common::FindVarsFromH5(core::IO &io, hid_t top_id, const char *gname, const char *heritage,
                                size_t ts)
{
    // int i = 0;
    // std::string stepStr;
    hsize_t numObj;

    hid_t gid = H5Gopen2(top_id, gname, H5P_DEFAULT);
    HDF5TypeGuard g(gid, E_H5_GROUP);
    ///    if (gid > 0) {
    herr_t ret = H5Gget_num_objs(gid, &numObj);
    if (ret >= 0)
    {
        hsize_t k = 0;
        char name[100];
        for (k = 0; k < numObj; k++)
        {
            ssize_t result = H5Gget_objname_by_idx(gid, (hsize_t)k, name, sizeof(name));
            if (result >= 0)
            {
                int currType = H5Gget_objtype_by_idx(gid, k);
                if ((currType == H5G_DATASET) || (currType == H5G_TYPE))
                {
                    std::string nameStr = name;
                    // if (!(0 == nameStr.find(PREFIX_BLOCKINFO)) &&
                    //  !(0 == nameStr.find(PREFIX_STAT)))
                    // cave in to cadacity requirement to pass pull request.
                    // This is odd
                    if (std::string::npos == nameStr.find(PREFIX_BLOCKINFO) &&
                        std::string::npos == nameStr.find(PREFIX_STAT))
                    {
                        hid_t datasetId = H5Dopen(gid, name, H5P_DEFAULT);
                        HDF5TypeGuard d(datasetId, E_H5_DATASET);

                        std::string longName;

                        if (strcmp(gname, "/") == 0)
                        {
                            longName = std::string("/") + name;
                        }
                        else
                        {
                            longName = std::string(heritage) + "/" + gname + "/" + name;
                        }
                        // CreateVar(io, datasetId, name);
                        ReadNativeAttrToIO(io, datasetId, longName);
                        CreateVar(io, datasetId, longName, ts);
                    }
                }
                else if (currType == H5G_GROUP)
                {
                    std::string heritageNext = heritage;
                    if (top_id != m_FileId)
                    {
                        heritageNext += "/";
                        heritageNext += gname;
                    }
                    FindVarsFromH5(io, gid, name, heritageNext.c_str(), ts);
                }
            }
        }
    }
}

// read variables from the input timestep
void HDF5Common::ReadVariables(size_t ts, core::IO &io)
{
    std::string stepStr;
    hsize_t numObj;

    StaticGetAdiosStepString(stepStr, ts);
    hid_t gid = H5Gopen2(m_FileId, stepStr.c_str(), H5P_DEFAULT);
    HDF5TypeGuard g(gid, E_H5_GROUP);

    herr_t ret = H5Gget_num_objs(gid, &numObj);
    if (ret >= 0)
    {
        hsize_t k = 0;
        char name[50];
        for (k = 0; k < numObj; k++)
        {
            ssize_t result = H5Gget_objname_by_idx(gid, (hsize_t)k, name, sizeof(name));
            if (result >= 0)
            {
                int currType = H5Gget_objtype_by_idx(gid, k);
                if (currType == H5G_GROUP)
                {
                    FindVarsFromH5(io, gid, name, "", ts);
                }
                else if ((currType == H5G_DATASET) || (currType == H5G_TYPE))
                {
                    std::string nameStr = name;
                    // if (!(0 == nameStr.find(PREFIX_BLOCKINFO)) &&
                    //  !(0 == nameStr.find(PREFIX_STAT)))
                    // cave in to cadacity requirement to pass pull request.
                    // This is odd
                    if (std::string::npos == nameStr.find(PREFIX_BLOCKINFO) &&
                        std::string::npos == nameStr.find(PREFIX_STAT))

                    {
                        hid_t datasetId = H5Dopen(gid, name, H5P_DEFAULT);
                        HDF5TypeGuard d(datasetId, E_H5_DATASET);
                        ReadNativeAttrToIO(io, datasetId, name);
                        CreateVar(io, datasetId, name, ts);
                    }
                }
            }
        }
    }
}

void HDF5Common::AddSingleString(core::IO &io, std::string const &name, hid_t datasetId, size_t ts)
{
    try
    {
        auto &foo = io.DefineVariable<std::string>(name);
        // 0 is a dummy holder. Just to make sure the ts entry is in there
        foo.m_AvailableStepBlockIndexOffsets[ts + 1] = std::vector<size_t>({0});
        foo.m_AvailableStepsStart = ts;
        // default was set to 0 while m_AvailabelStepsStart is 1.
        // correcting

        if (0 == foo.m_AvailableStepsCount)
        {
            foo.m_AvailableStepsCount++;
        }
    }
    catch (std::exception &e)
    {
        // invalid variable, do not define
        printf("WARNING: IO is not accepting definition of variable: %s. "
               "\nException: %s. \nSkipping. \n",
               name.c_str(), e.what());
    }
}

void HDF5Common::AddVarString(core::IO &io, std::string const &name, hid_t datasetId, size_t ts)
{
    core::Variable<std::string> *v = io.InquireVariable<std::string>(name);
    if (v != NULL)
    {
        v->m_AvailableStepsCount++;
        v->m_AvailableStepBlockIndexOffsets[ts + 1] = std::vector<size_t>({0});
        return;
    }

    // create a new string var if it is a single value
    hid_t dspace = H5Dget_space(datasetId);
    const int ndims = H5Sget_simple_extent_ndims(dspace);
    std::vector<hsize_t> dims(ndims);
    H5Sget_simple_extent_dims(dspace, dims.data(), NULL);
    H5Sclose(dspace);

    if ((ndims > 0))
    {
        bool isSingleElement = true;
        for (int i = 0; i < ndims; i++)
        {
            if (dims[i] > 1)
            {
                isSingleElement = false;
                break;
            }
        } // for

        if (!isSingleElement)
        {
            if (ndims != 2)
            {
                printf("WARNING: IO is not accepting definition of array "
                       "string variable: %s. dim: %d "
                       "Skipping! \n",
                       name.c_str(), ndims);
                return;
            }

            // adding special case from a dataset requirement
            for (unsigned long i = 0; i < dims[0]; i++)
                for (unsigned long j = 0; j < dims[1]; j++)
                {
                    std::string curr = name + "__" + std::to_string(i) + "_" + std::to_string(j);
                    // AddSingleString(io, curr, datasetId, ts);
                    try
                    {
                        auto &foo = io.DefineVariable<std::string>(curr);
                        // 0 is a dummy holder. Just to make sure the ts entry
                        // is in there
                        foo.m_AvailableStepBlockIndexOffsets[ts + 1] = std::vector<size_t>({0});
                        foo.m_AvailableStepsStart = ts;
                        // default was set to 0 while m_AvailabelStepsStart
                        // is 1. correcting
                        foo.m_Start = {i, j};
                        foo.m_Count = {1, 1};
                        if (0 == foo.m_AvailableStepsCount)
                        {
                            foo.m_AvailableStepsCount++;
                        }
                    }
                    catch (std::exception &e)
                    {
                        // invalid variable, do not define
                        printf("WARNING: IO is not accepting definition of "
                               "variable: %s. \nException: %s. \n"
                               "Skipping. \n",
                               name.c_str(), e.what());
                    }
                }
            return;
        }
    }

    // Dims shape;
    // shape.resize(ndims);
    AddSingleString(io, name, datasetId, ts);
}

template <class T>
void HDF5Common::AddVar(core::IO &io, std::string const &name, hid_t datasetId, size_t ts)
{
    core::Variable<T> *v = io.InquireVariable<T>(name);
    if (NULL == v)
    {
        hid_t dspace = H5Dget_space(datasetId);
        const int ndims = H5Sget_simple_extent_ndims(dspace);
        std::vector<hsize_t> dims(ndims);
        H5Sget_simple_extent_dims(dspace, dims.data(), NULL);
        H5Sclose(dspace);

        Dims shape;
        shape.resize(ndims);
        if (ndims > 0)
        {
            bool isOrderC = (io.m_ArrayOrder == ArrayOrdering::RowMajor);
            for (int i = 0; i < ndims; i++)
            {
                if (isOrderC)
                {
                    shape[i] = dims[i];
                }
                else
                {
                    shape[i] = dims[ndims - 1 - i];
                }
            }
        }

        Dims zeros(shape.size(), 0);

        try
        {
            auto &foo = io.DefineVariable<T>(name, shape, zeros, shape);
            // 0 is a dummy holder. Just to make sure the ts entry is in there
            foo.m_AvailableStepBlockIndexOffsets[ts + 1] = std::vector<size_t>({0});
            foo.m_AvailableStepsStart = ts;
            // default was set to 0 while m_AvailabelStepsStart is 1.
            // correcting

            if (0 == foo.m_AvailableStepsCount)
            {
                foo.m_AvailableStepsCount++;
            }
        }
        catch (std::exception &e)
        {
            // invalid variable, do not define
            printf("WARNING: IO is not accepting definition of variable: %s. "
                   "\nException: %s. \nSkipping. \n",
                   name.c_str(), e.what());
        }
    }
    else
    {
        /*    if (0 == v->m_AvailableStepsCount) { // default was set to 0 while
        m_AvailabelStepsStart is 1. v->m_AvailableStepsCount ++;
        }
        */
        v->m_AvailableStepsCount++;
        v->m_AvailableStepBlockIndexOffsets[ts + 1] = std::vector<size_t>({0});
    }
}

void HDF5Common::CreateVar(core::IO &io, hid_t datasetId, std::string const &nameSuggested,
                           size_t ts)
{
    std::string name;
    ReadADIOSName(datasetId, name);
    if (name.size() == 0)
    {
        name = nameSuggested;
    }

    hid_t h5Type = H5Dget_type(datasetId);
    HDF5TypeGuard t(h5Type, E_H5_DATATYPE);

    // int8 is mapped to "signed char" by IO.
    // so signed needs to be checked before char, otherwise type is char instead
    // of "signed char". Inqvar considers "signed char" and "char" different
    // types and returns error

    if (H5Tget_class(h5Type) == H5T_STRING)
    // if (H5Tequal(H5T_STRING, h5Type))
    {
        // AddVar<std::string>(io, name, datasetId, ts);
        // in ADIOS, a string var has to be a single value (not arrays)
        // so we treat 1D-array of 1 string in hdf5 as a single value in adios
        // otherwise direct mapping between hdf5 dataspace and adios global
        // dimension will fail for string type
        AddVarString(io, name, datasetId, ts);
        return;
    }

    if (H5Tequal(H5T_NATIVE_INT8, h5Type))
    {
        AddVar<int8_t>(io, name, datasetId, ts);
    }
    else if (H5Tequal(H5T_NATIVE_UINT8, h5Type))
    {
        AddVar<uint8_t>(io, name, datasetId, ts);
    }
    else if (H5Tequal(H5T_NATIVE_INT16, h5Type))
    {
        AddVar<int16_t>(io, name, datasetId, ts);
    }
    else if (H5Tequal(H5T_NATIVE_UINT16, h5Type))
    {
        AddVar<uint16_t>(io, name, datasetId, ts);
    }
    else if (H5Tequal(H5T_NATIVE_INT32, h5Type))
    {
        AddVar<int32_t>(io, name, datasetId, ts);
    }
    else if (H5Tequal(H5T_NATIVE_UINT32, h5Type))
    {
        AddVar<uint32_t>(io, name, datasetId, ts);
    }
    else if (H5Tequal(H5T_NATIVE_INT64, h5Type))
    {
        AddVar<int64_t>(io, name, datasetId, ts);
    }
    else if (H5Tequal(H5T_NATIVE_UINT64, h5Type))
    {
        AddVar<uint64_t>(io, name, datasetId, ts);
    }
    else if (H5Tequal(H5T_NATIVE_FLOAT, h5Type))
    {
        AddVar<float>(io, name, datasetId, ts);
    }
    else if (H5Tequal(H5T_NATIVE_DOUBLE, h5Type))
    {
        AddVar<double>(io, name, datasetId, ts);
    }
    else if (H5Tequal(H5T_NATIVE_LDOUBLE, h5Type))
    {
        AddVar<long double>(io, name, datasetId, ts);
    }
    else if (H5Tequal(m_DefH5TypeComplexFloat, h5Type))
    {
        AddVar<std::complex<float>>(io, name, datasetId, ts);
    }
    else if (H5Tequal(m_DefH5TypeComplexDouble, h5Type))
    {
        AddVar<std::complex<double>>(io, name, datasetId, ts);
    }
    else if (H5Tequal(m_DefH5TypeComplexLongDouble, h5Type))
    {
        // TODO:AddVar<std::complex<long double>>(io, name, datasetId, ts);
    }

    // H5Tclose(h5Type);
}

void HDF5Common::Close()
{
    if (m_FileId < 0)
    {
        return;
    }

    WriteAdiosSteps();

    if (m_GroupId >= 0)
    {
        H5Gclose(m_GroupId);
    }

    // close defined types
    // although H5Fclose will clean them anyways.
    H5Tclose(m_DefH5TypeComplexLongDouble);
    H5Tclose(m_DefH5TypeComplexDouble);
    H5Tclose(m_DefH5TypeComplexFloat);

    H5Pclose(m_PropertyTxfID);

    if (-1 != m_ChunkPID)
        H5Pclose(m_ChunkPID);

    H5Fclose(m_FileId);

    m_FileId = -1;
    m_GroupId = -1;
}

void HDF5Common::SetAdiosStep(size_t step)
{
    if (m_WriteMode)
        helper::Throw<std::ios_base::failure>("Toolkit", "interop::hdf5::HDF5Common",
                                              "SetAdiosStep",
                                              "unable to change step at Write MODE");

    if (step < 0)
        helper::Throw<std::ios_base::failure>("Toolkit", "interop::hdf5::HDF5Common",
                                              "SetAdiosStep", "unable to change to negative step");

    GetNumAdiosSteps();

    if (step >= m_NumAdiosSteps)
        helper::Throw<std::ios_base::failure>("Toolkit", "interop::hdf5::HDF5Common",
                                              "SetAdiosStep",
                                              "given time step is more than actual known steps");

    if (m_CurrentAdiosStep == step)
    {
        return;
    }

    if (m_GroupId >= 0)
        H5Gclose(m_GroupId);

    std::string stepName;
    StaticGetAdiosStepString(stepName, step);
    m_GroupId = H5Gopen(m_FileId, stepName.c_str(), H5P_DEFAULT);
    if (m_GroupId < 0)
    {
        helper::Throw<std::ios_base::failure>(
            "Toolkit", "interop::hdf5::HDF5Common", "SetAdiosStep",
            "ERROR: unable to open HDF5 group " + stepName + ", in call to Open");
    }

    m_CurrentAdiosStep = step;
}

//
// This function is intend to pair with CreateVarsFromIO().
//
// Because of the collective call requirement, we creata
// all variables through CreateVarFromIO() at BeginStep().
// At EndStep(), this function is called to remove unwritten variables
// to comply with ADIOS custom that define all vars, use a few per step
// and only see these few at the step.
//
// note that this works with 1 rank currently.
// need to find a general way in parallel HDF5 to
// detect whether a dataset called H5Dwrite  or not
//
void HDF5Common::CleanUpNullVars(core::IO &io)
{
    if (!m_WriteMode)
        return;

    if (m_CommSize != 1) // because the H5D_storage_size > 0 after H5Dcreate
                         // when there are multiple processors
        return;

    const core::VarMap &variables = io.GetVariables();
    for (const auto &vpair : variables)
    {
        const std::string &varName = vpair.first;
        const DataType varType = vpair.second->m_Type;
#define declare_template_instantiation(T)                                                          \
    if (varType == helper::GetDataType<T>())                                                       \
    {                                                                                              \
        core::Variable<T> *v = io.InquireVariable<T>(varName);                                     \
        if (!v)                                                                                    \
            return;                                                                                \
        RemoveEmptyDataset(varName);                                                               \
    }
        ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation
    }
}

void HDF5Common::Advance()
{
    if (m_WriteMode)
        CheckWriteGroup();

    if (m_GroupId >= 0)
    {
        H5Gclose(m_GroupId);
        m_GroupId = -1;
    }

    if (m_WriteMode)
    {
        // m_GroupId = H5Gcreate2(m_FileId, tsname.c_str(), H5P_DEFAULT,
        //                       H5P_DEFAULT, H5P_DEFAULT);
    }
    else
    {
        if (0 == m_NumAdiosSteps)
        {
            GetNumAdiosSteps();
        }
        if ((m_CurrentAdiosStep + 1) >= m_NumAdiosSteps)
        {
            return;
        }

        // std::string stepName =
        //    "/AdiosStep" + std::to_string(m_CurrentAdiosStep + 1);
        std::string stepName;
        StaticGetAdiosStepString(stepName, m_CurrentAdiosStep + 1);
        m_GroupId = H5Gopen(m_FileId, stepName.c_str(), H5P_DEFAULT);
        if (m_GroupId < 0)
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "interop::hdf5::HDF5Common", "Advance",
                                                  "unable to open HDF5 group " + stepName +
                                                      ", in call to Open");
        }
    }
    ++m_CurrentAdiosStep;
}

void HDF5Common::CheckWriteGroup()
{
    if (!m_WriteMode)
    {
        return;
    }
    if (m_GroupId >= 0)
    {
        return;
    }

    // std::string stepName = "/AdiosStep" +
    // std::to_string(m_CurrentAdiosStep);
    std::string stepName;
    StaticGetAdiosStepString(stepName, m_CurrentAdiosStep);
    m_GroupId = H5Gcreate2(m_FileId, stepName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (m_GroupId < 0)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "interop::hdf5::HDF5Common",
                                              "CheckWriteGroup",
                                              "Unable to create group " + stepName);
    }
}

void HDF5Common::ReadStringScalarDataset(hid_t dataSetId, std::string &result)
{
    hid_t h5Type = H5Dget_type(dataSetId); // get actual type;
    size_t typesize = H5Tget_size(h5Type);

    if (H5Tis_variable_str(h5Type))
    {
        hid_t d_space = H5Dget_space(dataSetId);
        std::vector<char> vc(typesize); // byte buffer to vlen strings
        auto status = H5Dread(dataSetId, h5Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, vc.data());
        CHECK_H5_RETURN(status, "ReadStringScalar_variable_str");
        result.assign(*((char **)vc.data()));

        // free dynamically allocated vlen memory from H5Dread
        // recent versions shift to use H5Treclaim()
        H5Dvlen_reclaim(h5Type, d_space, H5P_DEFAULT, vc.data());

        // H5Treclaim(attr_type, attr_space, H5P_DEFAULT, vc.data());
    }
    else
    {
        char *val = (char *)(calloc(typesize, sizeof(char)));
        hid_t ret2 = H5Dread(dataSetId, h5Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, val);
        CHECK_H5_RETURN(ret2, "ReadStringScalarDataset");

        result.assign(val, typesize);
        free(val);
    }
    H5Tclose(h5Type);
}

hid_t HDF5Common::GetTypeStringScalar(const std::string &input)
{
    /* Create a datatype to refer to. */
    hid_t type = H5Tcopy(H5T_C_S1);
    hid_t ret = H5Tset_size(type, input.size());
    CHECK_H5_RETURN(ret, "GetTypeStringScalar, H5Tset_size");
    ret = H5Tset_strpad(type, H5T_STR_NULLTERM);
    CHECK_H5_RETURN(ret, "GetTypeStringScalar, H5Tset_strpad");
    return type;
}

void HDF5Common::CreateDataset(const std::string &varName, hid_t h5Type, hid_t filespaceID,
                               std::vector<hid_t> &datasetChain)
{
    std::vector<std::string> list;
    char delimiter = '/';
    int delimiterLength = 1;
    std::string s = std::string(varName);
    size_t pos = 0;
    std::string token;
    while ((pos = s.find(delimiter)) != std::string::npos)
    {
        if (pos > 0)
        { // "///a/b/c" == "a/b/c"
            token = s.substr(0, pos);
            list.push_back(token);
        }
        s.erase(0, pos + delimiterLength);
    }
    list.push_back(s);

    hid_t topId = m_GroupId;
    if (list.size() > 1)
    {
        for (size_t i = 0; i < list.size() - 1; i++)
        {
            if (H5Lexists(topId, list[i].c_str(), H5P_DEFAULT) == 0)
            { // does not exist, so create
                topId = H5Gcreate2(topId, list[i].c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            }
            else
            {
                topId = H5Gopen(topId, list[i].c_str(), H5P_DEFAULT);
            }
            datasetChain.push_back(topId);
        }
    }

    hid_t varCreateProperty = H5P_DEFAULT;
    if (-1 != m_ChunkPID)
    {
        if (m_ChunkVarNames.size() == 0) // applies to all var
            varCreateProperty = m_ChunkPID;
        else if (m_ChunkVarNames.find(varName) != m_ChunkVarNames.end())
            varCreateProperty = m_ChunkPID;
    }

    /*
    hid_t dsetID = H5Dcreate(topId, list.back().c_str(), h5Type, filespaceID,
                             H5P_DEFAULT, varCreateProperty, H5P_DEFAULT);
    */

    hid_t dsetID = -1;
    if (H5Lexists(topId, list.back().c_str(), H5P_DEFAULT) == 0)
    {
        dsetID = H5Dcreate(topId, list.back().c_str(), h5Type, filespaceID, H5P_DEFAULT,
                           varCreateProperty, H5P_DEFAULT);
        if (list.back().compare(varName) != 0)
        {
            StoreADIOSName(varName, dsetID); // only stores when not the same
        }
    }
    else
        dsetID = H5Dopen(topId, list.back().c_str(), H5P_DEFAULT);

    datasetChain.push_back(dsetID);
}

void HDF5Common::StoreADIOSName(const std::string adiosName, hid_t dsetID)
{
    hid_t attrSpace = H5Screate(H5S_SCALAR);
    hid_t atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, adiosName.size());
    H5Tset_strpad(atype, H5T_STR_NULLTERM);
    hid_t attr = H5Acreate2(dsetID, ATTRNAME_GIVEN_ADIOSNAME.c_str(), atype, attrSpace, H5P_DEFAULT,
                            H5P_DEFAULT);
    H5Awrite(attr, atype, adiosName.c_str());

    H5Sclose(attrSpace);
    H5Tclose(atype);
    H5Aclose(attr);
}

void HDF5Common::ReadADIOSName(hid_t dsetID, std::string &adiosName)
{
    // htri_t H5Lexists( hid_t loc_id, const char *name, hid_t lapl_id )
    if (H5Aexists(dsetID, ATTRNAME_GIVEN_ADIOSNAME.c_str()) <= 0)
    {
        return;
    }

    hid_t attrID = H5Aopen(dsetID, ATTRNAME_GIVEN_ADIOSNAME.c_str(), H5P_DEFAULT);
    if (attrID < 0)
    {
        return;
    }

    hid_t attrType = H5Aget_type(attrID);
    size_t typeSize = H5Tget_size(attrType);

    char *val = (char *)(calloc(typeSize, sizeof(char)));
    hid_t ret2 = H5Aread(attrID, attrType, val);

    H5Tclose(attrType);
    H5Aclose(attrID);

    CHECK_H5_RETURN(ret2, "ReadADIOSName");
    adiosName.assign(val, typeSize);
    free(val);
}

bool HDF5Common::OpenDataset(const std::string &varName, std::vector<hid_t> &datasetChain)
{
    std::vector<std::string> list;
    char delimiter = '/';
    int delimiterLength = 1;
    std::string s = std::string(varName);
    size_t pos = 0;
    std::string token;
    while ((pos = s.find(delimiter)) != std::string::npos)
    {
        if (pos > 0)
        { // "///a/b/c" == "a/b/c"
            token = s.substr(0, pos);
            list.push_back(token);
        }
        s.erase(0, pos + delimiterLength);
    }
    list.push_back(s);

    if (list.size() == 1)
    {
        if (H5Lexists(m_GroupId, list[0].c_str(), H5P_DEFAULT) == 0)
        {
            datasetChain.push_back(-1);
            return false;
        }
        else
        {
            hid_t dsetID = H5Dopen(m_GroupId, list[0].c_str(), H5P_DEFAULT);
            datasetChain.push_back(dsetID);
            return true;
        }
    }

    hid_t topId = m_GroupId;

    for (size_t i = 0; i < list.size() - 1; i++)
    {
        if (H5Lexists(topId, list[i].c_str(), H5P_DEFAULT) == 0)
        { // does not exist, err
            // topId = H5Gcreate2(topId, list[i].c_str(), H5P_DEFAULT,
            // H5P_DEFAULT,H5P_DEFAULT);
            printf("Unable to open HDF5 group: %s for %s. Quit. \n", list[i].c_str(),
                   varName.c_str());
            return false;
        }
        else
        {
            topId = H5Gopen(topId, list[i].c_str(), H5P_DEFAULT);
        }
        datasetChain.push_back(topId);
    }

    hid_t dsetID = H5Dopen(topId, list.back().c_str(), H5P_DEFAULT);

    datasetChain.push_back(dsetID);
    return true;
}

//
// We use H5Dget_storage_size to see whether H5Dwrite has been called.
// and looks like when there are multiple processors, H5Dcreate causes
// storage allcoation already. So we limit this function when there is only
// one rank.
//
void HDF5Common::RemoveEmptyDataset(const std::string &varName)
{
    if (m_CommSize > 1)
        return;

    std::vector<std::string> list;
    char delimiter = '/';
    int delimiterLength = 1;
    std::string s = std::string(varName);
    size_t pos = 0;
    std::string token;
    while ((pos = s.find(delimiter)) != std::string::npos)
    {
        if (pos > 0)
        { // "///a/b/c" == "a/b/c"
            token = s.substr(0, pos);
            list.push_back(token);
        }
        s.erase(0, pos + delimiterLength);
    }
    list.push_back(s);

    if (list.size() == 1)
    {
        if (H5Lexists(m_GroupId, list[0].c_str(), H5P_DEFAULT) != 0)
        {
            hid_t dsetID = H5Dopen(m_GroupId, list[0].c_str(), H5P_DEFAULT);
            HDF5TypeGuard d(dsetID, E_H5_DATASET);

            H5D_space_status_t status;
            herr_t s1 = H5Dget_space_status(dsetID, &status);
            CHECK_H5_RETURN(s1, "RemoveEmptyDataset");
            if (0 == H5Dget_storage_size(dsetID)) /*nothing is written */
                H5Ldelete(m_GroupId, list[0].c_str(), H5P_DEFAULT);
        }
        return;
    }

    hid_t topId = m_GroupId;
    std::vector<hid_t> datasetChain;

    for (size_t i = 0; i < list.size() - 1; i++)
    {
        if (H5Lexists(topId, list[i].c_str(), H5P_DEFAULT) == 0)
            break;
        else
            topId = H5Gopen(topId, list[i].c_str(), H5P_DEFAULT);

        datasetChain.push_back(topId);
    }
    hid_t dsetID = H5Dopen(topId, list.back().c_str(), H5P_DEFAULT);
    datasetChain.push_back(dsetID);

    HDF5DatasetGuard g(datasetChain);

    if (H5Lexists(topId, list.back().c_str(), H5P_DEFAULT) != 0)
    {
        if (0 == H5Dget_storage_size(dsetID)) // nothing is written
            H5Ldelete(topId, list.back().c_str(), H5P_DEFAULT);
    }
}

// trim from right
inline std::string &rtrim(std::string &s, const char *t = " \t\n\r\f\v")
{
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

// trim from left
inline std::string &ltrim(std::string &s, const char *t = " \t\n\r\f\v")
{
    s.erase(0, s.find_first_not_of(t));
    return s;
}

// trim from left & right
inline std::string &trim(std::string &s, const char *t = " \t\n\r\f\v")
{
    return ltrim(rtrim(s, t), t);
}

void HDF5Common::ReadInStringAttr(core::IO &io, const std::string &attrName, hid_t attrId,
                                  hid_t h5Type, hid_t sid)
{
    hsize_t typeSize = H5Tget_size(h5Type);
    H5S_class_t stype = H5Sget_simple_extent_type(sid);

    if (H5Tis_variable_str(h5Type))
    {
        hid_t attr_space = H5Aget_space(attrId);
        if (H5S_SCALAR == stype)
        {
            std::vector<char> vc(typeSize); // byte buffer to vlen strings
            auto status = H5Aread(attrId, h5Type, vc.data());
            CHECK_H5_RETURN(status, "ReadInStringAttr_scalar")

            std::string c_str(*((char **)vc.data()));

            // free dynamically allocated vlen memory from H5Aread
            // later versions use H5Treclaim() instead.
            H5Dvlen_reclaim(h5Type, attr_space, H5P_DEFAULT, vc.data());

            io.DefineAttribute<std::string>(attrName, c_str);
        }
        else
        {
            hsize_t ndims = H5Sget_simple_extent_ndims(sid);
            if (ndims != 1)
                CHECK_H5_RETURN(-1, "Only handles 1-D string array");

            // ndims must be 1
            hsize_t dims[1];
            hid_t ret = H5Sget_simple_extent_dims(sid, dims, NULL);
            CHECK_H5_RETURN(ret, "ReadInStringAttr");

            std::vector<char *> vc(dims[0]);
            auto status = H5Aread(attrId, h5Type, vc.data());
            CHECK_H5_RETURN(status, "ReadInStringAttr");

            std::vector<std::string> stringArray;
            for (auto const &val : vc)
                // stringArray.push_back(auxiliary::strip(std::string(val), {'\0'}));
                stringArray.push_back(std::string(val));

            status = H5Dvlen_reclaim(h5Type, attr_space, H5P_DEFAULT, vc.data());
            io.DefineAttribute<std::string>(attrName, stringArray.data(), dims[0]);
        }

        return;
    }
    //
    // regular string, not variables
    //
    if (H5S_SCALAR == stype)
    {
        auto val = std::unique_ptr<char[]>(new char[typeSize]);
        H5Aread(attrId, h5Type, &val[0]);

        auto strValue = std::string(&val[0], typeSize);
        io.DefineAttribute<std::string>(attrName, strValue);
    }
    else
    { // array
        hsize_t ndims = H5Sget_simple_extent_ndims(sid);
        if (ndims != 1)
        {
            return; // so far io can only handle 1-D array
        }
        // ndims must be 1
        hsize_t dims[1];
        hid_t ret = H5Sget_simple_extent_dims(sid, dims, NULL);
        CHECK_H5_RETURN(ret, "ReadInStringAttr");
        auto val = std::unique_ptr<char[]>(new char[typeSize * dims[0]]);
        H5Aread(attrId, h5Type, val.get());

        std::vector<std::string> stringArray;
        for (hsize_t i = 0; i < dims[0]; i++)
        {
            auto input = std::string(&val[i * typeSize], typeSize);
            // remove the padded empty space;
            rtrim(input);
            stringArray.push_back(input);
        }

        io.DefineAttribute<std::string>(attrName, stringArray.data(), dims[0]);
    }
}

template <class T>
void HDF5Common::AddNonStringAttribute(core::IO &io, std::string const &attrName, hid_t attrId,
                                       hid_t h5Type, hsize_t arraySize)
{
    if (arraySize == 0)
    { // SCALAR
        T val;
        H5Aread(attrId, h5Type, &val);
        io.DefineAttribute(attrName, val);
    }
    else
    {
        std::vector<T> val(arraySize);
        H5Aread(attrId, h5Type, val.data());
        io.DefineAttribute(attrName, val.data(), arraySize);
    }
}

void HDF5Common::ReadInNonStringAttr(core::IO &io, const std::string &attrName, hid_t attrId,
                                     hid_t h5Type, hid_t sid)
{
    hsize_t ndims = H5Sget_simple_extent_ndims(sid);

    if (ndims > 1)
    {
        return; // so far adios2 io can only handle 1-D array
    }

    hsize_t dims[1];
    dims[0] = 0;
    if (ndims == 1)
    {
        hid_t ret = H5Sget_simple_extent_dims(sid, dims, NULL);
        CHECK_H5_RETURN(ret, "ReadInNonStringAttr");
    }

    if (H5Tequal(H5T_NATIVE_INT8, h5Type))
    {
        AddNonStringAttribute<int8_t>(io, attrName, attrId, h5Type, dims[0]);
    }
    else if (H5Tequal(H5T_NATIVE_UINT8, h5Type))
    {
        AddNonStringAttribute<uint8_t>(io, attrName, attrId, h5Type, dims[0]);
    }
    else if (H5Tequal(H5T_NATIVE_INT16, h5Type))
    {
        AddNonStringAttribute<int16_t>(io, attrName, attrId, h5Type, dims[0]);
    }
    else if (H5Tequal(H5T_NATIVE_UINT16, h5Type))
    {
        AddNonStringAttribute<uint16_t>(io, attrName, attrId, h5Type, dims[0]);
    }
    else if (H5Tequal(H5T_NATIVE_INT32, h5Type))
    {
        AddNonStringAttribute<int32_t>(io, attrName, attrId, h5Type, dims[0]);
    }
    else if (H5Tequal(H5T_NATIVE_UINT32, h5Type))
    {
        AddNonStringAttribute<uint32_t>(io, attrName, attrId, h5Type, dims[0]);
    }
    else if (H5Tequal(H5T_NATIVE_INT64, h5Type))
    {
        AddNonStringAttribute<int64_t>(io, attrName, attrId, h5Type, dims[0]);
    }
    else if (H5Tequal(H5T_NATIVE_UINT64, h5Type))
    {
        AddNonStringAttribute<uint64_t>(io, attrName, attrId, h5Type, dims[0]);
    }
    else if (H5Tequal(H5T_NATIVE_FLOAT, h5Type))
    {
        AddNonStringAttribute<float>(io, attrName, attrId, h5Type, dims[0]);
    }
    else if (H5Tequal(H5T_NATIVE_DOUBLE, h5Type))
    {
        AddNonStringAttribute<double>(io, attrName, attrId, h5Type, dims[0]);
    }
    else if (H5Tequal(H5T_NATIVE_LDOUBLE, h5Type))
    {
        AddNonStringAttribute<long double>(io, attrName, attrId, h5Type, dims[0]);
    }
}

void HDF5Common::WriteStringAttr(core::IO &io, core::Attribute<std::string> *adiosAttr,
                                 const std::string &attrName, hid_t parentID)
{
    // core::Attribute<std::string> *adiosAttr =
    // io.InquireAttribute<std::string>(attrName);

    if (adiosAttr == NULL)
    {
        return;
    }

    if (adiosAttr->m_IsSingleValue)
    {
        hid_t h5Type = GetTypeStringScalar(adiosAttr->m_DataSingleValue.data());
        hid_t s = H5Screate(H5S_SCALAR);
        hid_t attr = H5Acreate2(parentID, attrName.c_str(), h5Type, s, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, h5Type, (adiosAttr->m_DataSingleValue.data()));
        H5Sclose(s);
        H5Tclose(h5Type);
        H5Aclose(attr);
    }
    else if (adiosAttr->m_Elements >= 1)
    {
        // is array
        size_t max = 0;
        size_t idxWithMax = 0;
        for (size_t i = 0; i < adiosAttr->m_Elements; i++)
        {
            size_t curr = adiosAttr->m_DataArray[i].size();
            if (max < curr)
            {
                max = curr;
                idxWithMax = i;
            }
        }

        hid_t h5Type = GetTypeStringScalar(adiosAttr->m_DataArray[idxWithMax]);
        // std::vector<char> temp;
        std::string all;
        for (size_t i = 0; i < adiosAttr->m_Elements; i++)
        {
            std::string curr = adiosAttr->m_DataArray[i];
            curr.resize(max, ' ');
            all.append(curr);
        }

        hsize_t onedim[1] = {adiosAttr->m_Elements};
        hid_t s = H5Screate_simple(1, onedim, NULL);
        hid_t attr = H5Acreate2(parentID, attrName.c_str(), h5Type, s, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, h5Type, all.c_str());
        H5Sclose(s);
        H5Aclose(attr);
        H5Tclose(h5Type);
    }
}

template <class T>
void HDF5Common::WriteNonStringAttr(core::IO &io, core::Attribute<T> *adiosAttr, hid_t parentID,
                                    const char *h5AttrName)
{
    if (adiosAttr == NULL)
    {
        return;
    }
    hid_t h5Type = GetHDF5Type<T>();
    if (adiosAttr->m_IsSingleValue)
    {
        hid_t s = H5Screate(H5S_SCALAR);
        // hid_t attr = H5Acreate2(parentID, adiosAttr->m_Name.c_str(), h5Type,
        // s,
        hid_t attr = H5Acreate2(parentID, h5AttrName, h5Type, s, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, h5Type, &(adiosAttr->m_DataSingleValue));
        H5Sclose(s);
        H5Aclose(attr);
    }
    else if (adiosAttr->m_Elements >= 1)
    {
        hsize_t onedim[1] = {adiosAttr->m_Elements};
        hid_t s = H5Screate_simple(1, onedim, NULL);
        // hid_t attr = H5Acreate2(parentID, adiosAttr->m_Name.c_str(), h5Type,
        // s,
        hid_t attr = H5Acreate2(parentID, h5AttrName, h5Type, s, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, h5Type, adiosAttr->m_DataArray.data());
        H5Sclose(s);
        H5Aclose(attr);
    }
}

void HDF5Common::LocateAttrParent(const std::string &attrName, std::vector<std::string> &list,
                                  std::vector<hid_t> &parentChain)
{
    char delimiter = '/';
    int delimiterLength = 1;
    std::string s = std::string(attrName);
    size_t pos = 0;
    std::string token;
    while ((pos = s.find(delimiter)) != std::string::npos)
    {
        if (pos > 0)
        { // "///a/b/c" == "a/b/c"
            token = s.substr(0, pos);
            list.push_back(token);
        }
        s.erase(0, pos + delimiterLength);
    }
    list.push_back(s);

    if (list.size() == 1)
    {
        return;
    }

    hid_t topId = m_FileId;
    if (list.size() >= 1)
    {
        std::string ts;
        for (size_t i = 0; i < m_CurrentAdiosStep; i++)
        {
            StaticGetAdiosStepString(ts, i);
            for (size_t j = 0; j < list.size() - 1; j++)
            {
                ts += delimiter;
                ts += list[j].c_str();
            }
            if (H5Lexists(m_FileId, ts.c_str(), H5P_DEFAULT) <= 0)
                continue;
            else
            {
                topId = H5Dopen(m_FileId, ts.c_str(), H5P_DEFAULT);
                break;
            }
        } // for

        if (topId != m_FileId)
            parentChain.push_back(topId);
        return;
    } // if

    // hid_t dsetID = H5Dopen(topId, list.back().c_str(), H5P_DEFAULT);

    // parentChain.push_back(dsetID);
    // return dsetID;
}

//
// Create Variables from IO
//
void HDF5Common::CreateVarsFromIO(core::IO &io)
{
    if (!m_WriteMode)
        return;

    CheckWriteGroup(); // making sure all processors are creating new step

    if (!m_IdleWriterOn)
        return;

    const core::VarMap &variables = io.GetVariables();
    for (const auto &vpair : variables)
    {
        const std::string &varName = vpair.first;
        const DataType varType = vpair.second->m_Type;
#define declare_template_instantiation(T)                                                          \
    if (varType == helper::GetDataType<T>())                                                       \
    {                                                                                              \
        core::Variable<T> *v = io.InquireVariable<T>(varName);                                     \
        if (!v)                                                                                    \
            return;                                                                                \
        DefineDataset(*v);                                                                         \
    }
        ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation
    }
}
//
// write attr from io to hdf5
// right now adios only support global attr
// var does not have attr
//
void HDF5Common::WriteAttrFromIO(core::IO &io)
{
    if (m_FileId < 0)
    {
        return;
    }
    if (!m_WriteMode)
    {
        return;
    }

    const std::map<std::string, Params> &attributesInfo = io.GetAvailableAttributes();

    for (const auto &apair : attributesInfo)
    {
        std::string attrName = apair.first;
        Params temp = apair.second;
        DataType attrType = helper::GetDataTypeFromString(temp["Type"]);

        hid_t parentID = m_FileId;
#ifdef NO_ATTR_VAR_ASSOC
        std::vector<hid_t> chain;
        std::vector<std::string> list;
        LocateAttrParent(attrName, list, chain);
        HDF5DatasetGuard g(chain);

        if (chain.size() > 0)
        {
            parentID = chain.back();
        }
#else
        // will list out all attr at root level
        // to make it easy to be consistant with ADIOS2 attr symantic
        std::vector<std::string> list;
        list.push_back(attrName);
#endif
        // if (H5Aexists(parentID, attrName.c_str()) > 0)
        if (H5Aexists(parentID, list.back().c_str()) > 0)
        {
            continue;
        }

        if (attrType == DataType::Struct)
        {
            // not supported
        }
        else if (attrType == helper::GetDataType<std::string>())
        {
            // WriteStringAttr(io, attrName, parentID);
            core::Attribute<std::string> *adiosAttr = io.InquireAttribute<std::string>(attrName);
            WriteStringAttr(io, adiosAttr, list.back(), parentID);
        }
//
// note no std::complext attr types
//
#define declare_template_instantiation(T)                                                          \
    else if (attrType == helper::GetDataType<T>())                                                 \
    {                                                                                              \
        core::Attribute<T> *adiosAttr = io.InquireAttribute<T>(attrName);                          \
        WriteNonStringAttr(io, adiosAttr, parentID, list.back().c_str());                          \
    }
        ADIOS2_FOREACH_ATTRIBUTE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation
    }
}

//
// read attr from hdf5 to IO
// only the global ones are retrieved for now,
// until adios2 starts to support var level attrs
//
void HDF5Common::ReadAttrToIO(core::IO &io)
{
    hsize_t numAttrs;
    // herr_t ret = H5Gget_num_objs(m_FileId, &numObj);
    H5O_info_t oinfo;

#if H5_VERSION_GE(1, 11, 0)
#if (defined H5_USE_16_API || defined H5_USE_18_API)
    herr_t ret = H5Oget_info(m_FileId, &oinfo);
#else
    herr_t ret = H5Oget_info(m_FileId, &oinfo, H5O_INFO_ALL);
#endif
#else
    herr_t ret = H5Oget_info(m_FileId, &oinfo);
#endif
    if (ret >= 0)
    {
        numAttrs = oinfo.num_attrs;
        hsize_t k = 0;
        const int MAX_ATTR_NAME_SIZE = 100;
        for (k = 0; k < numAttrs; k++)
        {
            char attrName[MAX_ATTR_NAME_SIZE];
            ret = (herr_t)H5Aget_name_by_idx(m_FileId, ".", H5_INDEX_CRT_ORDER, H5_ITER_DEC,
                                             (hsize_t)k, attrName, (size_t)MAX_ATTR_NAME_SIZE,
                                             H5P_DEFAULT);
            if (ret >= 0)
            {
                // if (strcmp(attrName, ATTRNAME_NUM_STEPS.c_str()) == 0) {
                if (ATTRNAME_NUM_STEPS.compare(attrName) == 0)
                {
                    continue;
                }

                hid_t attrId = H5Aopen(m_FileId, attrName, H5P_DEFAULT);
                if (attrId < 0)
                {
                    continue;
                }
                HDF5TypeGuard ag(attrId, E_H5_ATTRIBUTE);

                hid_t sid = H5Aget_space(attrId);
                HDF5TypeGuard sg(sid, E_H5_SPACE);
                // H5S_class_t stype = H5Sget_simple_extent_type(sid);

                hid_t attrType = H5Aget_type(attrId);
                bool isString = (H5Tget_class(attrType) == H5T_STRING);
                if (isString)
                {
                    ReadInStringAttr(io, attrName, attrId, attrType, sid);
                }
                else
                {
                    ReadInNonStringAttr(io, attrName, attrId, attrType, sid);
                }
            }
        }
    }
}

void HDF5Common::ReadNativeAttrToIO(core::IO &io, hid_t datasetId, std::string const &pathFromRoot)
{
    hsize_t numAttrs;

    H5O_info_t oinfo;
#if H5_VERSION_GE(1, 11, 0)
#if (defined H5_USE_16_API || defined H5_USE_18_API)
    herr_t ret = H5Oget_info(datasetId, &oinfo);
#else
    herr_t ret = H5Oget_info(datasetId, &oinfo, H5O_INFO_ALL);
#endif
#else
    herr_t ret = H5Oget_info(datasetId, &oinfo);
#endif
    if (ret >= 0)
    {
        numAttrs = oinfo.num_attrs;

        if (numAttrs <= 0)
        {
            return; // warning: reading attrs at every var can be very time
                    // consuimg
        }
        hsize_t k = 0;
        const int MAX_ATTR_NAME_SIZE = 100;
        for (k = 0; k < numAttrs; k++)
        {
            char attrName[MAX_ATTR_NAME_SIZE];
            ret = (herr_t)H5Aget_name_by_idx(datasetId, ".", H5_INDEX_CRT_ORDER, H5_ITER_DEC,
                                             (hsize_t)k, attrName, (size_t)MAX_ATTR_NAME_SIZE,
                                             H5P_DEFAULT);
            if (ret >= 0)
            {
                hid_t attrId = H5Aopen(datasetId, attrName, H5P_DEFAULT);
                if (attrId < 0)
                {
                    continue;
                }
                HDF5TypeGuard ag(attrId, E_H5_ATTRIBUTE);
                if (ATTRNAME_GIVEN_ADIOSNAME.compare(attrName) == 0)
                {
                    continue;
                }

                hid_t sid = H5Aget_space(attrId);
                HDF5TypeGuard sg(sid, E_H5_SPACE);

                hid_t attrType = H5Aget_type(attrId);
                bool isString = (H5Tget_class(attrType) == H5T_STRING);

                std::string attrNameInAdios = pathFromRoot + "/" + attrName;
                if (isString)
                {
                    ReadInStringAttr(io, attrNameInAdios, attrId, attrType, sid);
                }
                else
                {
                    ReadInNonStringAttr(io, attrNameInAdios, attrId, attrType, sid);
                }
            }
        }
    }
}

void HDF5Common::StaticGetAdiosStepString(std::string &stepName, size_t ts)
{
    stepName = "/Step" + std::to_string(ts);
}

void HDF5Common::CheckVariableOperations(const core::VariableBase &variable) const
{
    if (!variable.m_Operations.empty())
    {
        helper::Throw<std::runtime_error>("Toolkit", "interop::hdf5::HDF5Common",
                                          "CheckVariableOperations",
                                          "ADIOS2 Operators are not supported for HDF5 engine");
    }
}

#define declare_template_instantiation(T)                                                          \
    template void HDF5Common::Write(core::Variable<T> &, const T *);

ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

} // end namespace interop
} // end namespace adios

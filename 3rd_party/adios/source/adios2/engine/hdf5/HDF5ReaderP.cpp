/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * HDF5ReaderP.cpp
 *
 *  Created on: March 20, 2017
 *      Author: Junmin
 */

#include "HDF5ReaderP.h"
#include "HDF5ReaderP.tcc"

#include "adios2/helper/adiosFunctions.h" //CSVToVector
#include "adios2/helper/adiosFunctions.h" //IsHDF5
#include <limits>
#include <stdexcept>
#include <vector>

namespace adios2
{
namespace core
{
namespace engine
{

#define CHECK_H5_RETURN(returnCode, reason)                                                        \
    {                                                                                              \
        if (returnCode < 0)                                                                        \
        {                                                                                          \
            helper::Throw<std::runtime_error>("Engine", "HDF5ReaderP", "CHECK_H5_RETURN", reason); \
        }                                                                                          \
    }

HDF5ReaderP::HDF5ReaderP(IO &io, const std::string &name, const Mode openMode, helper::Comm comm)
: Engine("HDF5Reader", io, name, openMode, std::move(comm))
{
    if (!helper::IsHDF5File(name, io, m_Comm, {}))
    {
        helper::Throw<std::invalid_argument>("Engine", "HDF5ReaderP", "HDF5ReaderP",
                                             "Invalid HDF5 file found");
    }

    Init();
    m_IsOpen = true;
}

HDF5ReaderP::~HDF5ReaderP()
{
    if (IsValid())
        DoClose();
}

bool HDF5ReaderP::IsValid()
{
    bool isValid = false;

    if (m_OpenMode != Mode::Read)
    {
        return isValid;
    }
    if (m_H5File.m_FileId >= 0)
    {
        isValid = true;
    }
    return isValid;
}

void HDF5ReaderP::Init()
{
    if (m_OpenMode != Mode::Read)
    {
        helper::Throw<std::invalid_argument>("Engine", "HDF5ReaderP", "Init",
                                             "HDF5Reader only supports OpenMode::Read "
                                             ", in call to Open");
    }

    m_H5File.Init(m_Name, m_Comm, false);
    m_H5File.ParseParameters(m_IO);

    /*
     */
    m_H5File.ReadAttrToIO(m_IO);
    if (!m_InStreamMode)
    {
        // m_H5File.ReadVariables(0, m_IO);
        m_H5File.ReadAllVariables(m_IO);
    }
    else
    {
        // m_H5File.ReadVariables(0, m_IO);
        m_H5File.ReadAllVariables(m_IO);
    }
}

// returns slab size (>0) (datasize read in)
// returns -1 to advise do not continue
template <class T>
size_t HDF5ReaderP::ReadDataset(hid_t dataSetId, hid_t h5Type, Variable<T> &variable, T *values)
{
    hid_t fileSpace = H5Dget_space(dataSetId);
    interop::HDF5TypeGuard g_fs(fileSpace, interop::E_H5_SPACE);

    if (fileSpace < 0)
    {
        return 0;
    }

    size_t slabsize = 1u;

    size_t ndims = std::max(variable.m_Shape.size(), variable.m_Count.size());
    if (0u == ndims)
    { // is scalar
      // hid_t myclass = H5Tget_class(h5Type);
        if (H5Tget_class(h5Type) == H5T_STRING)
        {
            m_H5File.ReadStringScalarDataset(dataSetId, *(std::string *)values);
        }
        else
        {
            hid_t ret = H5Dread(dataSetId, h5Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values);
            CHECK_H5_RETURN(ret, "ReadDataset, H5Dread");
        }
    }
    else
    {
        std::vector<hsize_t> start(ndims), count(ndims), stride(ndims);
        bool isOrderC = (m_IO.m_ArrayOrder == ArrayOrdering::RowMajor);

        for (size_t i = 0u; i < ndims; i++)
        {
            if (isOrderC)
            {
                count[i] = variable.m_Count[i];
                start[i] = variable.m_Start[i];
            }
            else
            {
                count[i] = variable.m_Count[ndims - 1u - i];
                start[i] = variable.m_Start[ndims - 1u - i];
            }
            slabsize *= count[i];
            stride[i] = 1;
        }
        hid_t ret = H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, start.data(), stride.data(),
                                        count.data(), NULL);
        if (ret < 0)
            return 0;

        size_t max_int = static_cast<size_t>(std::numeric_limits<int>::max());
        if (ndims > max_int)
        {
            helper::Throw<std::overflow_error>("Engine", "HDF5ReaderP", "ReadDataset",
                                               "Number of dimensions is too large to be "
                                               "represented by an int");
        }

        hid_t memDataSpace = H5Screate_simple(static_cast<int>(ndims), count.data(), NULL);
        interop::HDF5TypeGuard g_mds(memDataSpace, interop::E_H5_SPACE);

        size_t elementsRead = 1;
        for (size_t i = 0u; i < ndims; i++)
        {
            elementsRead *= count[i];
        }

        if (H5Tget_class(h5Type) != H5T_STRING)
        {
            ret = H5Dread(dataSetId, h5Type, memDataSpace, fileSpace, H5P_DEFAULT, values);
            // H5Sclose(memDataSpace);
        }
        else if (1 == elementsRead)
        {
            hid_t h5Type = H5Dget_type(dataSetId); // get actual type;
            size_t typesize = H5Tget_size(h5Type);

            char *val = (char *)(calloc(typesize, sizeof(char)));
            H5Dread(dataSetId, h5Type, memDataSpace, fileSpace, H5P_DEFAULT, val);

            ((std::string *)values)->assign(val, typesize);
            free(val);

            H5Tclose(h5Type);
        }
    }

    return slabsize;
}

template <class T>
void HDF5ReaderP::UseHDFRead(Variable<T> &variable, T *data, hid_t h5Type)
{

    if (!m_H5File.m_IsGeneratedByAdios)
    {
        std::size_t found = variable.m_Name.find("__");
        if (found != std::string::npos)
        {
            std::string h5name = variable.m_Name.substr(0, found);
            hid_t dataSetId = H5Dopen(m_H5File.m_FileId, h5name.c_str(), H5P_DEFAULT);
            if (dataSetId < 0)
                return;
            interop::HDF5TypeGuard g_ds(dataSetId, interop::E_H5_DATASET);
            ReadDataset(dataSetId, h5Type, variable, data);
            return;
        }
        else
        {
            // UseHDFReadNativeFile(variable, data, h5Type);
            hid_t dataSetId = H5Dopen(m_H5File.m_FileId, variable.m_Name.c_str(), H5P_DEFAULT);
            if (dataSetId < 0)
            {
                return;
            }

            interop::HDF5TypeGuard g_ds(dataSetId, interop::E_H5_DATASET);
            ReadDataset(dataSetId, h5Type, variable, data);
            return;
        }
    }

    T *values = data;

    size_t ts = 0;
    size_t variableStart = variable.m_StepsStart;
    /*
      // looks like m_StepsStart is defaulted to be 0 now.
    if (!m_InStreamMode && (variableStart == 1))
    { // variableBase::m_StepsStart min=1
        variableStart = 0;
    }
    */

    while (ts < variable.m_StepsCount)
    {
        try
        {
            m_H5File.SetAdiosStep(variableStart + ts);
        }
        catch (std::exception &e)
        {
            printf("[Not fatal] %s\n", e.what());
            break;
        }

        std::vector<hid_t> chain;
        if (!m_H5File.OpenDataset(variable.m_Name, chain))
        {
            return;
        }
        hid_t dataSetId = chain.back();
        interop::HDF5DatasetGuard g(chain);

        if (dataSetId < 0)
        {
            return;
        }

        size_t slabsize = ReadDataset(dataSetId, h5Type, variable, values);

        if (slabsize == 0)
        {
            break;
        }

        ts++;
        values += slabsize;
    } // while
}

/*
 */

StepStatus HDF5ReaderP::BeginStep(StepMode mode, const float timeoutSeconds)
{
    const size_t ts = m_H5File.GetNumAdiosSteps();

    if (m_StreamAt >= ts)
    {
        m_IO.m_ReadStreaming = false;
        return StepStatus::EndOfStream;
    }

    if (m_DeferredStack.size() > 0)
    {
        // EndStep was not called!
        return StepStatus::NotReady;
    }

    if (m_InStreamMode && (m_StreamAt == m_IO.m_EngineStep))
    {
        // EndStep() was not called after BeginStep()
        // otherwise m_StreamAt would have been increated by 1
        return StepStatus::OtherError;
    }

    m_InStreamMode = true;

    m_IO.m_ReadStreaming = true;
    m_IO.m_EngineStep = m_StreamAt;

    return StepStatus::OK;
}

size_t HDF5ReaderP::CurrentStep() const { return m_StreamAt; }

void HDF5ReaderP::EndStep()
{
    if (m_DeferredStack.size() > 0)
    {
        PerformGets();
    }
    m_StreamAt++;
    m_H5File.Advance();
}

void HDF5ReaderP::PerformGets()
{

#define declare_type(T)                                                                            \
    for (std::string variableName : m_DeferredStack)                                               \
    {                                                                                              \
        const DataType type = m_IO.InquireVariableType(variableName);                              \
        if (type == helper::GetDataType<T>())                                                      \
        {                                                                                          \
            Variable<T> *var = m_IO.InquireVariable<T>(variableName);                              \
            if (var != nullptr)                                                                    \
            {                                                                                      \
                if (m_InStreamMode)                                                                \
                {                                                                                  \
                    var->m_StepsStart = m_StreamAt;                                                \
                    var->m_StepsCount = 1;                                                         \
                }                                                                                  \
                hid_t h5Type = m_H5File.GetHDF5Type<T>();                                          \
                UseHDFRead(*var, var->GetData(), h5Type);                                          \
            }                                                                                      \
        }                                                                                          \
    }
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    //
    // clear the list here so only read the variables once
    // (other functions like EndStep() will call PerformToGet() also)
    //
    m_DeferredStack.clear();
}

#define declare_type(T)                                                                            \
    void HDF5ReaderP::DoGetSync(Variable<T> &variable, T *data) { GetSyncCommon(variable, data); } \
                                                                                                   \
    void HDF5ReaderP::DoGetDeferred(Variable<T> &variable, T *data)                                \
    {                                                                                              \
        GetDeferredCommon(variable, data);                                                         \
    }                                                                                              \
                                                                                                   \
    std::map<size_t, std::vector<typename Variable<T>::BPInfo>> HDF5ReaderP::DoAllStepsBlocksInfo( \
        const Variable<T> &variable) const                                                         \
    {                                                                                              \
        return GetAllStepsBlocksInfo(variable);                                                    \
    }                                                                                              \
                                                                                                   \
    std::vector<typename Variable<T>::BPInfo> HDF5ReaderP::DoBlocksInfo(                           \
        const Variable<T> &variable, const size_t step) const                                      \
    {                                                                                              \
        return GetBlocksInfo(variable, step);                                                      \
    }

ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

void HDF5ReaderP::DoClose(const int transportIndex)
{
    /*
     */
    EndStep();
    /*
     */
    m_H5File.Close();
}

size_t HDF5ReaderP::DoSteps() const { return m_H5File.GetAdiosStep(); }
} // end namespace engine
} // end namespace core
} // end namespace adios2

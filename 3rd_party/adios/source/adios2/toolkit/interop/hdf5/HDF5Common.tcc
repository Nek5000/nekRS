/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * HDF5Common.tcc
 *
 *  Created on: Jun 1, 2017
 *      Author: Junmin
 */

#ifndef ADIOS2_TOOLKIT_INTEROP_HDF5_HDF5COMMON_TCC_
#define ADIOS2_TOOLKIT_INTEROP_HDF5_HDF5COMMON_TCC_

#include "HDF5Common.h"
#include <iostream>
#include <vector>

#include "adios2/helper/adiosFunctions.h"

namespace adios2
{
namespace interop
{

template <class T>
void HDF5Common::DefineDataset(core::Variable<T> &variable)
{
    CheckVariableOperations(variable);
    size_t dimSize = std::max(variable.m_Shape.size(), variable.m_Count.size());
    hid_t h5Type = GetHDF5Type<T>();

    if (dimSize == 0)
    {
        // write scalar
        hid_t filespaceID = H5Screate(H5S_SCALAR);
        HDF5TypeGuard g0(filespaceID, E_H5_SPACE);

        std::vector<hid_t> chain;
        CreateDataset(variable.m_Name, h5Type, filespaceID, chain);
        HDF5DatasetGuard g(chain);
        return;
    }
    // dimSize > 0
    std::vector<hsize_t> dimsf, count, offset;
    GetHDF5SpaceSpec(variable, dimsf, count, offset);

    size_t max_int = static_cast<size_t>(std::numeric_limits<int>::max());
    if (dimSize > max_int)
    {
        helper::Throw<std::overflow_error>("Toolkit", "interop::hdf5::HDF5Common", "DefineDataset",
                                           "dimSize is too large "
                                           "to be represented by an int");
    }

    hid_t fileSpace = H5Screate_simple(static_cast<int>(dimSize), dimsf.data(), NULL);
    HDF5TypeGuard fs(fileSpace, E_H5_SPACE);

    std::vector<hid_t> chain;
    CreateDataset(variable.m_Name, h5Type, fileSpace, chain);
    HDF5DatasetGuard g(chain);
}

template <>
void HDF5Common::DefineDataset<std::string>(core::Variable<std::string> &variable)
{
    std::cout << "...Needs actual string size, so defer to later? var name=" << variable.m_Name
              << std::endl;
}

template <class T>
void HDF5Common::GetHDF5SpaceSpec(const core::Variable<T> &variable, std::vector<hsize_t> &dimsf,
                                  std::vector<hsize_t> &count, std::vector<hsize_t> &offset)
{
    size_t dimSize = std::max(variable.m_Shape.size(), variable.m_Count.size());
    for (size_t i = 0; i < dimSize; ++i)
    {
        if (variable.m_Shape.size() == dimSize)
        {
            dimsf.push_back(variable.m_Shape[i]);
        }
        else
        {
            dimsf.push_back(variable.m_Count[i]);
        }

        if (variable.m_Count.size() == dimSize)
        {
            count.push_back(variable.m_Count[i]);
            if (variable.m_Start.size() == dimSize)
            {
                offset.push_back(variable.m_Start[i]);
            }
            else
            {
                offset.push_back(0);
            }
        }
        else
        {
            count.push_back(variable.m_Shape[i]);
            offset.push_back(0);
        }
    }

    if (dimSize <= 1)
        return;
    if (m_OrderByC)
        return;

    for (size_t i = 0; i < dimSize / 2; ++i)
    {
        std::swap(dimsf[i], dimsf[dimSize - 1 - i]);
        std::swap(count[i], count[dimSize - 1 - i]);
        std::swap(offset[i], offset[dimSize - 1 - i]);
    }
}

template <class T>
void HDF5Common::Write(core::Variable<T> &variable, const T *values)
{
    CheckWriteGroup();
    CheckVariableOperations(variable);
    size_t dimSize = std::max(variable.m_Shape.size(), variable.m_Count.size());
    hid_t h5Type = GetHDF5Type<T>();

    if (std::is_same<T, std::string>::value)
    {
        h5Type = GetTypeStringScalar(*(std::string *)values);
    }

    if (dimSize == 0)
    {
#ifndef RELAY_DEFINE_TO_HDF5 // RELAY_DEFINE_TO_HDF5 = variables in io are
                             // created at begin_step
        // write scalar
        hid_t filespaceID = H5Screate(H5S_SCALAR);

        std::vector<hid_t> chain;
        CreateDataset(variable.m_Name, h5Type, filespaceID, chain);
        HDF5DatasetGuard g(chain);
        hid_t dsetID = chain.back();

        if (std::is_same<T, std::string>::value)
        {
            H5Dwrite(dsetID, h5Type, H5S_ALL, H5S_ALL, m_PropertyTxfID,
                     ((std::string *)values)->data());
            H5Tclose(h5Type);
        }
        else
        {
            H5Dwrite(dsetID, h5Type, H5S_ALL, H5S_ALL, m_PropertyTxfID, values);
        }
        H5Sclose(filespaceID);
#else

        if (std::is_same<T, std::string>::value)
        {
            hid_t filespaceID = H5Screate(H5S_SCALAR);
            std::vector<hid_t> chain;
            CreateDataset(variable.m_Name, h5Type, filespaceID, chain);
            HDF5DatasetGuard g(chain);
            H5Sclose(filespaceID);
        }
        hid_t dsetID = H5Dopen(m_GroupId, variable.m_Name.c_str(), H5P_DEFAULT);
        if (-1 != dsetID)
        { // exists
            if (std::is_same<T, std::string>::value)
            {
                H5Dwrite(dsetID, h5Type, H5S_ALL, H5S_ALL, m_PropertyTxfID,
                         ((std::string *)values)->data());
                H5Tclose(h5Type);
            }
            else
            {
                H5Dwrite(dsetID, h5Type, H5S_ALL, H5S_ALL, m_PropertyTxfID, values);
            }
        }
#endif

        return;
    }

    std::vector<hsize_t> dimsf, count, offset;
    GetHDF5SpaceSpec(variable, dimsf, count, offset);

    size_t max_int = static_cast<size_t>(std::numeric_limits<int>::max());
    if (dimSize > max_int)
    {
        helper::Throw<std::overflow_error>("Toolkit", "interop::hdf5::HDF5Common", "Write",
                                           "dimSize is too large "
                                           "to be represented by an int");
    }

    hid_t fileSpace = H5Screate_simple(static_cast<int>(dimSize), dimsf.data(), NULL);
#ifndef RELAY_DEFINE_TO_HDF5 // RELAY_DEFINE_TO_HDF5 = variables in io are
                             // created at begin_step
    std::vector<hid_t> chain;
    CreateDataset(variable.m_Name, h5Type, fileSpace, chain);
    hid_t dsetID = chain.back();
    HDF5DatasetGuard g(chain);
#else
    hid_t dsetID = H5Dopen(m_GroupId, variable.m_Name.c_str(), H5P_DEFAULT);
    HDF5TypeGuard k(dsetID, E_H5_DATASET);
    if (-1 == dsetID)
    {
    }
#endif
    hid_t memSpace = H5Screate_simple(static_cast<int>(dimSize), count.data(), NULL);

    // Select hyperslab
    fileSpace = H5Dget_space(dsetID);
    H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, offset.data(), NULL, count.data(), NULL);

    herr_t status;

    if (!variable.m_MemoryStart.empty())
    {
        auto blockSize = helper::GetTotalSize(variable.m_Count);
        T *k = reinterpret_cast<T *>(calloc(blockSize, sizeof(T)));

        adios2::Dims zero(variable.m_Start.size(), 0);
        helper::CopyMemoryBlock(k, zero, variable.m_Count, true, values,
                                zero, // variable.m_Start,
                                variable.m_Count, true, false, Dims(), Dims(),
                                variable.m_MemoryStart, variable.m_MemoryCount);

        status = H5Dwrite(dsetID, h5Type, memSpace, fileSpace, m_PropertyTxfID, k);
        free(k);
    }
    else
    {
        status = H5Dwrite(dsetID, h5Type, memSpace, fileSpace, m_PropertyTxfID, values);
    }

    if (status < 0)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "interop::hdf5::HDF5Common", "Write",
                                              "HDF5 file Write failed");
    }

#ifdef NO_STAT
    size_t valuesSize = adios2::helper::GetTotalSize(variable.m_Count);
    T min, max;
    auto memSpace = variable.GetMemorySpace(values);
    adios2::helper::GetMinMaxThreads(values, valuesSize, min, max, 1, memSpace);

    int chainSize = chain.size();
    hid_t parentId = m_GroupId;
    if (chainSize > 1)
    {
        parentId = chain[chainSize - 2];
    }
    AddBlockInfo(variable, parentId);

    std::vector<T> stats = {min, max};
    AddStats(variable, parentId, stats);
#endif
    H5Sclose(fileSpace);
    H5Sclose(memSpace);
}

template <class T>
void HDF5Common::AddStats(const core::Variable<T> &variable, hid_t parentId, std::vector<T> &stats)
{

    hid_t h5Type = GetHDF5Type<T>();

    std::string statInfo_name = PREFIX_STAT + variable.m_Name;
    hsize_t numStat = stats.size(); // min, max etc
    hsize_t statDim[2] = {(hsize_t)m_CommSize, numStat};
    hid_t statSpace_id = H5Screate_simple(numStat, statDim, NULL);
    hid_t statId = H5Dcreate(parentId, statInfo_name.c_str(), h5Type, statSpace_id, H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
    hsize_t statLocalDim[1] = {numStat};
    hid_t statLocal_id = H5Screate_simple(1, statLocalDim, NULL);

    hsize_t statOffset[2] = {(hsize_t)m_CommRank, 0};
    hsize_t statCount[2] = {1, numStat};
    H5Sselect_hyperslab(statSpace_id, H5S_SELECT_SET, statOffset, NULL, statCount, NULL);

    H5Dwrite(statId, h5Type, statLocal_id, statSpace_id, m_PropertyTxfID, stats.data());

    H5Sclose(statLocal_id);
    H5Sclose(statSpace_id);
    H5Dclose(statId);
    //     H5Pclose(plistID);
}

template <class T>
void HDF5Common::AddBlockInfo(const core::Variable<T> &variable, hid_t parentId)
{
    int dimSize = std::max(variable.m_Shape.size(), variable.m_Count.size());
    hsize_t metaDim[2];
    metaDim[1] = dimSize * 2;
    metaDim[0] = m_CommSize;
    hid_t metaSpace_id = H5Screate_simple(2, metaDim, NULL);
    std::string blockInfo_name = PREFIX_BLOCKINFO + variable.m_Name;
    hid_t metaId = H5Dcreate(parentId, blockInfo_name.c_str(), H5T_NATIVE_HSIZE, metaSpace_id,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    std::vector<size_t> blocks(dimSize * 2);
    for (int i = 0; i < dimSize; i++)
    {
        blocks[i + dimSize] = variable.m_Count[i];
        blocks[i] = variable.m_Start[i];
    }
    hsize_t blockDim[1] = {(hsize_t)(dimSize * 2)};
    hid_t metaLocal_id = H5Screate_simple(1, blockDim, NULL);

    hsize_t metaOffset[2] = {(hsize_t)m_CommRank, 0};
    hsize_t metaCount[2] = {1, (hsize_t)(dimSize * 2)};
    H5Sselect_hyperslab(metaSpace_id, H5S_SELECT_SET, metaOffset, NULL, metaCount, NULL);

    H5Dwrite(metaId, H5T_NATIVE_HSIZE, metaLocal_id, metaSpace_id, m_PropertyTxfID, blocks.data());

    H5Sclose(metaLocal_id);
    H5Sclose(metaSpace_id);
    H5Dclose(metaId);

    //      H5Pclose(plistID);
}
//
//
template <>
hid_t HDF5Common::GetHDF5Type<std::string>()
{
    return H5T_C_S1;
}

template <>
hid_t HDF5Common::GetHDF5Type<char>()
{
    return H5T_NATIVE_UCHAR;
}

template <>
hid_t HDF5Common::GetHDF5Type<int8_t>()
{
    return H5T_NATIVE_INT8;
}

template <>
hid_t HDF5Common::GetHDF5Type<uint8_t>()
{
    return H5T_NATIVE_UINT8;
}

template <>
hid_t HDF5Common::GetHDF5Type<int16_t>()
{
    return H5T_NATIVE_INT16;
}
template <>
hid_t HDF5Common::GetHDF5Type<uint16_t>()
{
    return H5T_NATIVE_UINT16;
}
template <>
hid_t HDF5Common::GetHDF5Type<int32_t>()
{
    return H5T_NATIVE_INT32;
}
template <>
hid_t HDF5Common::GetHDF5Type<uint32_t>()
{
    return H5T_NATIVE_UINT32;
}
template <>
hid_t HDF5Common::GetHDF5Type<int64_t>()
{
    return H5T_NATIVE_INT64;
}
template <>
hid_t HDF5Common::GetHDF5Type<uint64_t>()
{
    return H5T_NATIVE_UINT64;
}
template <>
hid_t HDF5Common::GetHDF5Type<float>()
{
    return H5T_NATIVE_FLOAT;
}
template <>
hid_t HDF5Common::GetHDF5Type<double>()
{
    return H5T_NATIVE_DOUBLE;
}
template <>
hid_t HDF5Common::GetHDF5Type<long double>()
{
    return H5T_NATIVE_LDOUBLE;
}
template <>
hid_t HDF5Common::GetHDF5Type<std::complex<float>>()
{
    return m_DefH5TypeComplexFloat;
}
template <>
hid_t HDF5Common::GetHDF5Type<std::complex<double>>()
{
    return m_DefH5TypeComplexDouble;
}
template <>
hid_t HDF5Common::GetHDF5Type<std::complex<long double>>()
{
    return m_DefH5TypeComplexLongDouble;
}

} // end namespace interop
} // end namespace adios

#endif /* ADIOS2_TOOLKIT_INTEROP_HDF5_HDF5COMMON_H_ */

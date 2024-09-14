/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * RefactorMDR.cpp :
 *
 *  Created on: Sep 14, 2023
 *      Author: Norbert Podhorszki <pnorbert@ornl.gov>
 */

#include "RefactorMDR.h"
#include "adios2/helper/adiosFunctions.h"
#include "adios2/operator/compress/CompressNull.h"
#include <assert.h>
#include <cstring>

#include <stdio.h>

#include <mgard/compress_x.hpp>

namespace adios2
{
namespace core
{
namespace refactor
{

RefactorMDR::RefactorMDR(const Params &parameters)
: Operator("mdr", REFACTOR_MDR, "refactor", parameters)
{
    config.normalize_coordinates = false;
    config.log_level = 1;
    config.decomposition = mgard_x::decomposition_type::MultiDim;
    config.domain_decomposition = mgard_x::domain_decomposition_type::MaxDim;
    // config.domain_decomposition = mgard_x::domain_decomposition_type::Block;
    // config.block_size = 64;

    config.dev_type = mgard_x::device_type::AUTO;
    config.prefetch = false;
    // config.max_memory_footprint = max_memory_footprint;

    config.lossless = mgard_x::lossless_type::Huffman_Zstd;
    // mgard_x::lossless_type::Huffman, mgard_x::lossless_type::Huffman_Zstd,
    // mgard_x::lossless_type::CPU_Lossless

    // This should be changed to 4 for float later
    config.total_num_bitplanes = 64;
}

size_t RefactorMDR::GetEstimatedSize(const size_t ElemCount, const size_t ElemSize,
                                     const size_t ndims, const size_t *dims) const
{
    std::cout << "RefactorMDR::GetEstimatedSize() called \n";
    mgard_x::data_type mtype =
        (ElemSize == 8 ? mgard_x::data_type::Double : mgard_x::data_type::Float);
    DataType datatype = (ElemSize == 8 ? DataType::Double : DataType::Float);
    Dims dimsV(ndims);
    for (size_t i = 0; i < ndims; ++i)
    {
        dimsV[i] = dims[i];
    }
    Dims convertedDims = ConvertDims(dimsV, datatype, 3);
    mgard_x::DIM mgardDim = ndims;
    std::vector<mgard_x::SIZE> mgardCount;
    for (const auto &c : convertedDims)
    {
        mgardCount.push_back(c);
    }

    mgard_x::Config cfg(config); // copy of const config
    if (mtype == mgard_x::data_type::Float)
    {
        cfg.total_num_bitplanes = 32;
    }

    auto s = mgard_x::MDR::MDRMaxOutputDataSize(mgardDim, mtype, mgardCount, config);

    size_t sizeIn = helper::GetTotalSize(convertedDims, ElemSize);
    std::cout << "RefactorMDR Estimated Max output size = " << s << " for input size = " << sizeIn
              << std::endl;
    return (size_t)s + 128; // count in the header
};

struct RefactorMDRHeader
{
    uint64_t ndims;
    std::vector<uint64_t> dims; // ndims values
    adios2::DataType type;
    uint8_t mgard_version_major;
    uint8_t mgard_version_minor;
    uint8_t mgard_version_patch;
    bool wasRefactored;
    uint64_t metadataHeaderSize;
    uint64_t metadataSize;
    uint8_t nSubdomains;
    uint8_t nLevels;
    uint8_t nBitPlanes;
};

size_t RefactorMDR::Operate(const char *dataIn, const Dims &blockStart, const Dims &blockCount,
                            const DataType type, char *bufferOut)
{
    const uint8_t bufferVersion = 1;
    size_t bufferOutOffset = 0;

    MakeCommonHeader(bufferOut, bufferOutOffset, bufferVersion);

    Dims convertedDims = ConvertDims(blockCount, type, 3);

    const size_t ndims = convertedDims.size();
    if (ndims > 5)
    {
        helper::Throw<std::invalid_argument>("Operator", "RefactorMDR", "Operate",
                                             "MGARD does not support data in " +
                                                 std::to_string(ndims) + " dimensions");
    }

    // mgard V1 metadata
    PutParameter(bufferOut, bufferOutOffset, ndims);
    for (const auto &d : convertedDims)
    {
        PutParameter(bufferOut, bufferOutOffset, d);
    }
    PutParameter(bufferOut, bufferOutOffset, type);
    PutParameter(bufferOut, bufferOutOffset, static_cast<uint8_t>(MGARD_VERSION_MAJOR));
    PutParameter(bufferOut, bufferOutOffset, static_cast<uint8_t>(MGARD_VERSION_MINOR));
    PutParameter(bufferOut, bufferOutOffset, static_cast<uint8_t>(MGARD_VERSION_PATCH));
    // mgard V1 metadata end

    // set type
    mgard_x::data_type mgardType;
    if (type == helper::GetDataType<float>())
    {
        mgardType = mgard_x::data_type::Float;
        config.total_num_bitplanes = 32;
    }
    else if (type == helper::GetDataType<double>())
    {
        mgardType = mgard_x::data_type::Double;
        config.total_num_bitplanes = 64;
    }
    else if (type == helper::GetDataType<std::complex<float>>())
    {
        mgardType = mgard_x::data_type::Float;
        config.total_num_bitplanes = 32;
    }
    else if (type == helper::GetDataType<std::complex<double>>())
    {
        mgardType = mgard_x::data_type::Double;
        config.total_num_bitplanes = 64;
    }
    else
    {
        helper::Throw<std::invalid_argument>("Operator", "RefactorMDR", "Operate",
                                             "MGARD only supports float and double types");
    }
    // set type end

    // set mgard style dim info
    mgard_x::DIM mgardDim = ndims;
    std::vector<mgard_x::SIZE> mgardCount;
    for (const auto &c : convertedDims)
    {
        mgardCount.push_back(c);
    }
    // set mgard style dim info end

    // input under this size will not be refactored
    const size_t thresholdSize = 100000;

    size_t sizeIn = helper::GetTotalSize(blockCount, helper::GetDataTypeSize(type));

    if (sizeIn < thresholdSize)
    {
        /* disable compression and add marker in the header*/
        PutParameter(bufferOut, bufferOutOffset, false);
        headerSize = bufferOutOffset;
        return 0;
    }
    PutParameter(bufferOut, bufferOutOffset, true);

    mgard_x::MDR::RefactoredMetadata refactored_metadata;
    mgard_x::MDR::RefactoredData refactored_data;
    char *dataptr = const_cast<char *>(dataIn);
    mgard_x::pin_memory(dataptr, sizeIn, config);
    mgard_x::MDR::MDRefactor(mgardDim, mgardType, mgardCount, dataIn, refactored_metadata,
                             refactored_data, config, false);
    mgard_x::unpin_memory(dataptr, config);

    size_t nbytes = SerializeRefactoredData(refactored_metadata, refactored_data,
                                            bufferOut + bufferOutOffset, SIZE_MAX);

    bufferOutOffset += nbytes;

    return bufferOutOffset;
}

// return number of bytes written
size_t RefactorMDR::SerializeRefactoredData(mgard_x::MDR::RefactoredMetadata &refactored_metadata,
                                            mgard_x::MDR::RefactoredData &refactored_data,
                                            char *buffer, size_t maxsize)
{
    size_t offset = 0;
    size_t MDRHeaderSize = 0;

    /* Metadata header */
    const uint64_t metadata_header_size = refactored_metadata.header.size();
    {
        PutParameter<uint64_t>(buffer, offset, metadata_header_size);
        std::memcpy(buffer + offset, refactored_metadata.header.data(), metadata_header_size);
        offset += metadata_header_size;
    }

    /* Metadata */
    std::vector<mgard_x::Byte> serialized_metadata = refactored_metadata.Serialize();
    const uint64_t metadata_size = (uint64_t)serialized_metadata.size();
    {
        PutParameter<uint64_t>(buffer, offset, metadata_size);
        std::memcpy(buffer + offset, serialized_metadata.data(), metadata_size);
        offset += metadata_size;
    }

    std::cout << "MDR metadata seralized " << offset << " bytes. header = " << metadata_header_size
              << " metadata = " << metadata_size << "\n";

    /* 3D table of subdomain X level X bitplane offsets, not all of them has data */
    uint8_t nSubdomains = refactored_metadata.metadata.size();
    uint8_t nLevels = 0;
    uint8_t nBitPlanes = 0;

    for (size_t subdomain_id = 0; subdomain_id < refactored_metadata.metadata.size();
         subdomain_id++)
    {
        nLevels = std::max(nLevels,
                           (uint8_t)refactored_metadata.metadata[subdomain_id].level_sizes.size());
        for (size_t level_idx = 0;
             level_idx < refactored_metadata.metadata[subdomain_id].level_sizes.size(); level_idx++)
        {
            nBitPlanes = std::max(
                nBitPlanes,
                (uint8_t)refactored_metadata.metadata[subdomain_id].level_sizes[level_idx].size());
        }
    }

    PutParameter<uint8_t>(buffer, offset, nSubdomains);
    PutParameter<uint8_t>(buffer, offset, nLevels);
    PutParameter<uint8_t>(buffer, offset, nBitPlanes);
    uint64_t tableSize = (nSubdomains * nLevels * nBitPlanes);
    uint64_t *table = (uint64_t *)(buffer + offset);
    std::fill(table, table + tableSize, 0ULL);
    offset += tableSize * sizeof(uint64_t);
    MDRHeaderSize = offset;

    /* Individual components of refactored data */
    size_t nBlocks = 0;
    uint64_t tableIdx = 0;
    for (size_t subdomain_id = 0; subdomain_id < refactored_metadata.metadata.size();
         subdomain_id++)
    {
        for (size_t level_idx = 0;
             level_idx < refactored_metadata.metadata[subdomain_id].level_sizes.size(); level_idx++)
        {
            tableIdx = subdomain_id * nLevels * nBitPlanes + level_idx * nBitPlanes;
            for (size_t bitplane_idx = 0;
                 bitplane_idx <
                 refactored_metadata.metadata[subdomain_id].level_sizes[level_idx].size();
                 bitplane_idx++)
            {
                std::memcpy(buffer + offset,
                            refactored_data.data[subdomain_id][level_idx][bitplane_idx],
                            refactored_metadata.metadata[subdomain_id]
                                .level_sizes[level_idx][bitplane_idx]);
                table[tableIdx] = offset - MDRHeaderSize;
                offset +=
                    refactored_metadata.metadata[subdomain_id].level_sizes[level_idx][bitplane_idx];
                ++nBlocks;
                ++tableIdx;
            }
        }
    }

    std::cout << "MDR serialized " << offset << " bytes, MDR header size = " << MDRHeaderSize
              << " subdomains = " << (size_t)nSubdomains << " levels = " << (size_t)nLevels
              << " bitplanes = " << (size_t)nBitPlanes << " blocks = " << nBlocks << "\n";
    return offset;
}

size_t RefactorMDR::GetHeaderSize() const { return headerSize; }

size_t RefactorMDR::ReconstructV1(const char *bufferIn, const size_t sizeIn, char *dataOut)
{
    // Do NOT remove even if the buffer version is updated. Data might be still
    // in lagacy formats. This function must be kept for backward compatibility.
    // If a newer buffer format is implemented, create another function, e.g.
    // ReconstructV1 and keep this function for reconstructing legacy data.

    config.log_level = 1;
    // double s = std::numeric_limits<double>::infinity();

    size_t bufferInOffset = 0;

    const size_t ndims = GetParameter<size_t, size_t>(bufferIn, bufferInOffset);
    Dims blockCount(ndims);
    for (size_t i = 0; i < ndims; ++i)
    {
        blockCount[i] = GetParameter<size_t, size_t>(bufferIn, bufferInOffset);
    }
    const DataType type = GetParameter<DataType>(bufferIn, bufferInOffset);
    m_VersionInfo = " Data is compressed using MGARD Version " +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) + "." +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) + "." +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) +
                    ". Please make sure a compatible version is used for decompression.";

    const bool isRefactored = GetParameter<bool>(bufferIn, bufferInOffset);

    if (!isRefactored)
    {
        // data was copied as is from this offset
        headerSize += bufferInOffset;
        m_AccuracyProvided = m_AccuracyRequested;
        return 0;
    }

    size_t sizeOut = helper::GetTotalSize(blockCount, helper::GetDataTypeSize(type));

    mgard_x::MDR::RefactoredMetadata refactored_metadata;

    /* Metadata header */
    {
        const uint64_t metadata_header_size = GetParameter<uint64_t>(bufferIn, bufferInOffset);
        refactored_metadata.header.resize(metadata_header_size);
        std::memcpy(refactored_metadata.header.data(), bufferIn + bufferInOffset,
                    metadata_header_size);
        bufferInOffset += metadata_header_size;
    }

    /* Metadata */
    {
        const uint64_t metadata_size = GetParameter<uint64_t>(bufferIn, bufferInOffset);
        std::vector<mgard_x::Byte> serialized_metadata;
        serialized_metadata.resize(metadata_size);
        std::memcpy(serialized_metadata.data(), bufferIn + bufferInOffset, metadata_size);
        bufferInOffset += metadata_size;
        refactored_metadata.Deserialize(serialized_metadata);
        refactored_metadata.InitializeForReconstruction();
    }

    mgard_x::MDR::RefactoredData refactored_data;
    refactored_data.InitializeForReconstruction(refactored_metadata);

    const uint8_t nSubdomains = GetParameter<uint8_t>(bufferIn, bufferInOffset);
    const uint8_t nLevels = GetParameter<uint8_t>(bufferIn, bufferInOffset);
    const uint8_t nBitPlanes = GetParameter<uint8_t>(bufferIn, bufferInOffset);

    const uint64_t tableSize = (nSubdomains * nLevels * nBitPlanes);
    const uint64_t *table = (uint64_t *)(bufferIn + bufferInOffset);
    bufferInOffset += tableSize * sizeof(uint64_t);

    /*
        Reconstruction
    */
    const char *componentData = bufferIn + bufferInOffset; // data pieces are here in memory

    mgard_x::MDR::ReconstructedData reconstructed_data;
    for (auto &metadata : refactored_metadata.metadata)
    {
        metadata.requested_tol = m_AccuracyRequested.error;
        metadata.requested_s = m_AccuracyRequested.norm;
    }
    mgard_x::MDR::MDRequest(refactored_metadata, config);
    /*for (auto &metadata : refactored_metadata.metadata)
    {
        metadata.PrintStatus();
    }*/

    bool first_reconstruction = true; // only will be needed with progressive reconstruction

    // Assemble data pieces from buffer
    {
        uint64_t tableIdx;
        int num_subdomains = refactored_metadata.metadata.size();
        assert(nSubdomains == refactored_metadata.metadata.size());
        for (int subdomain_id = 0; subdomain_id < num_subdomains; subdomain_id++)
        {
            mgard_x::MDR::MDRMetadata metadata = refactored_metadata.metadata[subdomain_id];
            int num_levels = metadata.level_sizes.size();
            for (int level_idx = 0; level_idx < num_levels; level_idx++)
            {
                assert(nLevels >= metadata.level_sizes.size());
                tableIdx = subdomain_id * nLevels * nBitPlanes + level_idx * nBitPlanes;
                int num_bitplanes = metadata.level_sizes[level_idx].size();
                int loaded_bitplanes = metadata.loaded_level_num_bitplanes[level_idx];
                int reqested_bitplanes = metadata.requested_level_num_bitplanes[level_idx];
                assert(nBitPlanes >= metadata.requested_level_num_bitplanes[level_idx]);
                for (int bitplane_idx = loaded_bitplanes; bitplane_idx < reqested_bitplanes;
                     bitplane_idx++)
                {

                    uint64_t componentSize = refactored_metadata.metadata[subdomain_id]
                                                 .level_sizes[level_idx][bitplane_idx];
                    const mgard_x::Byte *cdata =
                        reinterpret_cast<const mgard_x::Byte *>(componentData + table[tableIdx]);

                    refactored_data.data[subdomain_id][level_idx][bitplane_idx] =
                        const_cast<mgard_x::Byte *>(cdata);

                    ++tableIdx;
                }
                if (first_reconstruction)
                {
                    // initialize level signs
                    refactored_data.level_signs[subdomain_id][level_idx] =
                        (bool *)malloc(sizeof(bool) * metadata.level_num_elems[level_idx]);
                    std::memset(refactored_data.level_signs[subdomain_id][level_idx], 0,
                                sizeof(bool) * metadata.level_num_elems[level_idx]);
                }
            }
        }
    }

    /* Initialize reconstructed data manually here to force using
        user allocated memory for the final result
    */
    if (false)
    {
        reconstructed_data.Initialize(1);
        reconstructed_data.data[0] = reinterpret_cast<mgard_x::Byte *>(dataOut);
        std::memset(reconstructed_data.data[0], 0, sizeOut);
        std::vector<mgard_x::SIZE> offsets = std::vector<mgard_x::SIZE>(ndims, 0);
        std::vector<mgard_x::SIZE> shape = std::vector<mgard_x::SIZE>();
        for (const auto &c : blockCount)
        {
            shape.push_back(static_cast<mgard_x::SIZE>(c));
        }
        reconstructed_data.offset[0] = offsets;
        reconstructed_data.shape[0] = shape;
    }

    mgard_x::MDR::MDReconstruct(refactored_metadata, refactored_data, reconstructed_data, config,
                                false);

    size_t refactoredSize =
        helper::GetTotalSize(reconstructed_data.shape[0], helper::GetDataTypeSize(type));
    assert(sizeOut == refactoredSize);
    std::memcpy(dataOut, reconstructed_data.data[0], sizeOut);

    first_reconstruction = false;

    /*for (int subdomain_id = 0; subdomain_id < reconstructed_data.data.size(); subdomain_id++)
    {
        std::cout << "reconstructed_data " << subdomain_id << " : offset(";
        for (auto n : reconstructed_data.offset[subdomain_id])
        {
            std::cout << n << " ";
        }
        std::cout << ") shape(";
        for (auto n : reconstructed_data.shape[subdomain_id])
        {
            std::cout << n << " ";
        }
        std::cout << ")\n";
    }*/

    if (type == DataType::FloatComplex || type == DataType::Float)
    {
        m_AccuracyProvided.error =
            std::max(m_AccuracyRequested.error, adios2::ops::mdr::FLOAT_ROUNDING_ERROR_LIMIT);
    }
    else
    {
        m_AccuracyProvided.error =
            std::max(m_AccuracyRequested.error, adios2::ops::mdr::DOUBLE_ROUNDING_ERROR_LIMIT);
    }
    m_AccuracyProvided.norm = m_AccuracyRequested.norm;
    m_AccuracyProvided.relative = false; // should be m_AccuracyRequested.relative

    if (type == DataType::FloatComplex || type == DataType::DoubleComplex)
    {
        sizeOut /= 2;
    }
    return sizeOut;
}

size_t RefactorMDR::InverseOperate(const char *bufferIn, const size_t sizeIn, char *dataOut)
{

    for (auto &p : m_Parameters)
    {
        std::cout << "User parameter " << p.first << " = " << p.second << std::endl;
        const std::string key = helper::LowerCase(p.first);
        // const std::string value = helper::LowerCase(p.second);
        if (key == "accuracy")
        {
            m_AccuracyRequested.error =
                helper::StringTo<double>(p.second, " in Parameter key=" + key);
            std::cout << "Accuracy error set from Parameter to " << m_AccuracyRequested.error;
        }
    }

    size_t bufferInOffset = 1; // skip operator type
    const uint8_t bufferVersion = GetParameter<uint8_t>(bufferIn, bufferInOffset);
    bufferInOffset += 2; // skip two reserved bytes
    headerSize = bufferInOffset;

    if (bufferVersion == 1)
    {
        return ReconstructV1(bufferIn + bufferInOffset, sizeIn - bufferInOffset, dataOut);
    }
    /*else if (bufferVersion == 2)
    {
        // TODO: if a Version 2 mgard buffer is being implemented, put it here
        // and keep the ReconstructV1 routine for backward compatibility
    }*/
    else
    {
        helper::Throw<std::runtime_error>("Operator", "RefactorMDR", "InverseOperate",
                                          "invalid mgard buffer version" + bufferVersion);
    }

    return 0;
}

bool RefactorMDR::IsDataTypeValid(const DataType type) const
{
    if (type == DataType::Double || type == DataType::Float || type == DataType::DoubleComplex ||
        type == DataType::FloatComplex)
    {
        return true;
    }
    return false;
}

} // end namespace compress
} // end namespace core
} // end namespace adios2

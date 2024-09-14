/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IO_ph5.cpp
 *
 *  Created on: Feb 2017
 *      Author: Norbert Podhorszki
 */

#include "IO.h"

#include <hdf5.h>
#include <ios>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

class HDF5NativeWriter
{

public:
    HDF5NativeWriter(const std::string &fileName);
    ~HDF5NativeWriter();

    bool Advance();
    void Close();
    void CheckWriteGroup();

    void WriteScalar(const std::string &varName, const void *data, hid_t h5Type);
    void WriteSimple(const std::string &varName, int dimSize, const void *data, hid_t h5Type,
                     const hsize_t *shape, const hsize_t *offset, const hsize_t *count);

    void applyMetadataCacheEviction();
    void WriteSimpleWithChunking(const std::string &varName, int dimSize, const void *data,
                                 hid_t h5Type, const hsize_t *shape, const hsize_t *offset,
                                 const hsize_t *count);
    int m_CurrentTimeStep;
    unsigned int m_TotalTimeSteps;
    bool m_Chunking;

private:
    hid_t m_FilePropertyListId;
    hid_t m_FileId;
    hid_t m_GroupId;
};

HDF5NativeWriter::HDF5NativeWriter(const std::string &fileName)
: m_CurrentTimeStep(0), m_TotalTimeSteps(0), m_Chunking(false)
{
    m_FilePropertyListId = H5Pcreate(H5P_FILE_ACCESS);

    // read a file collectively
    H5Pset_fapl_mpio(m_FilePropertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

    m_FileId = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, m_FilePropertyListId);

    if (m_FileId < 0)
    {
        throw std::runtime_error("Unable to open " + fileName + " for reading");
    }

    applyMetadataCacheEviction();

    std::string ts0 = "/Step0";

    m_GroupId = H5Gcreate2(m_FileId, ts0.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (m_GroupId < 0)
    {
        throw std::runtime_error("HDF5: Unable to create group " + ts0);
    }

    std::string prefix = "chunking";
    if (fileName.substr(0, prefix.size()) == prefix)
    {
        m_Chunking = true;
    }
}

HDF5NativeWriter::~HDF5NativeWriter() { Close(); }

void HDF5NativeWriter::applyMetadataCacheEviction()
{
#ifdef NEVER
    /*
      see
      https://lists.hdfgroup.org/pipermail/hdf-forum_lists.hdfgroup.org/2011-February/004201.html
      John said the code below worked for the paper but not anymore after
      updates
    */
    H5AC_cache_config_t mdc_config;
    mdc_config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
    H5Pget_mdc_config(m_FilePropertyListId, &mdc_config);
    mdc_config.evictions_enabled = 0; // FALSE
    mdc_config.incr_mode = H5C_incr__off;
    mdc_config.decr_mode = H5C_decr__off;
    H5Pset_mdc_config(m_FilePropertyListId, &mdc_config);
#endif
}

void HDF5NativeWriter::Close()
{
    if (m_FileId < 0)
        return;

    hid_t s = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate(m_FileId, "NumSteps", H5T_NATIVE_UINT, s, H5P_DEFAULT, H5P_DEFAULT);
    unsigned int totalTimeSteps = m_CurrentTimeStep + 1;

    if (m_GroupId < 0)
    {
        totalTimeSteps = m_CurrentTimeStep;
    }
    H5Awrite(attr, H5T_NATIVE_UINT, &totalTimeSteps);
    H5Sclose(s);
    H5Aclose(attr);

    if (m_GroupId >= 0)
    {
        H5Gclose(m_GroupId);
        m_GroupId = -1;
    }

    H5Fclose(m_FileId);
    m_FileId = -1;
    H5Pclose(m_FilePropertyListId);
}

bool HDF5NativeWriter::Advance()
{
    if (m_GroupId >= 0)
    {
        H5Gclose(m_GroupId);
        m_GroupId = -1;
    }

    ++m_CurrentTimeStep;

    return true;
}

void HDF5NativeWriter::CheckWriteGroup()
{
    if (m_GroupId >= 0)
    {
        return;
    }

    std::string timeStepName = "/Step" + std::to_string(m_CurrentTimeStep);
    m_GroupId = H5Gcreate2(m_FileId, timeStepName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (m_GroupId < 0)
    {
        throw std::runtime_error("HDF5: Unable to create group " + timeStepName);
    }
}

void HDF5NativeWriter::WriteScalar(const std::string &varName, const void *data, hid_t h5Type)
{
    CheckWriteGroup();
    // scalar
    hid_t filespaceID = H5Screate(H5S_SCALAR);
    hid_t dsetID = H5Dcreate(m_GroupId, varName.c_str(), h5Type, filespaceID, H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
    hid_t plistID = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plistID, H5FD_MPIO_COLLECTIVE);

    herr_t status = H5Dwrite(dsetID, h5Type, H5S_ALL, H5S_ALL, plistID, data);

    if (status < 0)
    {
        std::cerr << " Write failed. " << std::endl;
    }
    H5Sclose(filespaceID);
    H5Dclose(dsetID);
}

void HDF5NativeWriter::WriteSimple(const std::string &varName, int dimSize, const void *data,
                                   hid_t h5Type, const hsize_t *shape, const hsize_t *offset,
                                   const hsize_t *count)
{
    CheckWriteGroup();
    hid_t fileSpace = H5Screate_simple(dimSize, shape, NULL);

    hid_t dsetID = H5Dcreate(m_GroupId, varName.c_str(), h5Type, fileSpace, H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);

    hid_t memSpace = H5Screate_simple(dimSize, count, NULL);

    // Select hyperslab
    fileSpace = H5Dget_space(dsetID);
    H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, offset, NULL, count, NULL);

    //  Create property list for collective dataset write.

    hid_t plistID = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plistID, H5FD_MPIO_COLLECTIVE);

    herr_t status;

    status = H5Dwrite(dsetID, h5Type, memSpace, fileSpace, plistID, data);

    if (status < 0)
    {
        // error
        std::cerr << " Write failed. " << std::endl;
    }

    H5Dclose(dsetID);
    H5Sclose(fileSpace);
    H5Sclose(memSpace);
    H5Pclose(plistID);
}

void HDF5NativeWriter::WriteSimpleWithChunking(const std::string &varName, int dimSize,
                                               const void *data, hid_t h5Type, const hsize_t *shape,
                                               const hsize_t *offset, const hsize_t *count)
{
    CheckWriteGroup();
    hid_t fileSpace = H5Screate_simple(dimSize, shape, NULL);

    hid_t dsetPid = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dsetPid, dimSize, count);

    size_t bytes = H5Tget_size(h5Type);
    for (int i = 0; i < dimSize; i++)
    {
        bytes *= count[i];
    }
    hid_t access_plistid = H5Pcreate(H5P_DATASET_ACCESS);
    H5Pset_chunk_cache(access_plistid, 101, bytes, 1);

    hid_t dsetID = H5Dcreate(m_GroupId, varName.c_str(), h5Type, fileSpace, H5P_DEFAULT, dsetPid,
                             access_plistid);

    hid_t memSpace = H5Screate_simple(dimSize, count, NULL);

    // Select hyperslab
    fileSpace = H5Dget_space(dsetID);
    H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, offset, NULL, count, NULL);

    //  Create property list for collective dataset write.

    hid_t plistID = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plistID, H5FD_MPIO_COLLECTIVE);

    herr_t status;

    status = H5Dwrite(dsetID, h5Type, memSpace, fileSpace, plistID, data);

    if (status < 0)
    {
        // error
        std::cerr << " Write failed. " << std::endl;
    }

    H5Dclose(dsetID);
    H5Sclose(fileSpace);
    H5Sclose(memSpace);
    H5Pclose(plistID);

    H5Pclose(dsetPid);
    H5Pclose(access_plistid);
}

//
//
std::shared_ptr<HDF5NativeWriter> h5writer;
// HDF5NativeWriter* h5writer;

IO::IO(const Settings &s, MPI_Comm comm)
{
    m_outputfilename = MakeFilename(s.outputfile, ".h5");

    if (s.outputfile[0] == '0')
    {
        std::cout << " no writer. " << std::endl;
        h5writer = nullptr;
        return;
    }
    h5writer = std::make_shared<HDF5NativeWriter>(m_outputfilename);

    if (h5writer == nullptr)
        throw std::ios_base::failure("ERROR: failed to open ADIOS h5writer\n");
}

IO::~IO()
{
    if (h5writer != nullptr)
    {
        h5writer->Close();
    }
    // delete h5writer;
}

void IO::write(int step, const HeatTransfer &ht, const Settings &s, MPI_Comm comm)
{
    if (h5writer == nullptr)
    {
        return;
    }
    std::vector<hsize_t> dims = {s.gndx, s.gndy};
    std::vector<hsize_t> offset = {s.offsx, s.offsy};
    std::vector<hsize_t> count = {s.ndx, s.ndy};

    if (h5writer->m_Chunking)
    {
        h5writer->WriteSimpleWithChunking("T", 2, ht.data_noghost().data(), H5T_NATIVE_DOUBLE,
                                          dims.data(), offset.data(), count.data());
    }
    else
    {
        h5writer->WriteSimple("T", 2, ht.data_noghost().data(), H5T_NATIVE_DOUBLE, dims.data(),
                              offset.data(), count.data());
    }

    h5writer->WriteScalar("gndy", &(s.gndy), H5T_NATIVE_UINT);
    h5writer->WriteScalar("gndx", &(s.gndx), H5T_NATIVE_UINT);

    h5writer->Advance();
}

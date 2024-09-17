/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <iostream>
#include <stdexcept>

#include <adios2.h>

#include <hdf5.h>

#include <gtest/gtest.h>

#include "../adios2/engine/SmallTestData.h"
#include <h5vol/H5Vol.h>

class H5VolWriteReadTest : public ::testing::Test
{
public:
    H5VolWriteReadTest() = default;

    SmallTestData m_TestData;
};

class HDF5NativeReader
{

public:
    HDF5NativeReader(const std::string fileName);
    ~HDF5NativeReader();

    bool Advance();

    void GetVarInfo(const std::string varName, std::vector<hsize_t> &dims, hid_t &h5Type);
    // If offset, count and memspaceSize are provided, then variable would be
    // read by selection
    void ReadString(const std::string varName, std::string &result);
    void ReadVar(const std::string varName, void *dataArray, hsize_t *offset = nullptr,
                 hsize_t *count = nullptr, const size_t memsspaceSize = 0);

    int m_CurrentTimeStep;
    unsigned int m_TotalTimeSteps;

private:
    hid_t m_FilePropertyListId;
    hid_t m_FileId;
    hid_t m_GroupId;
};

class HDF5NativeWriter
{
public:
#ifdef TEST_HDF5_MPI
    HDF5NativeWriter(const std::string &fileName, MPI_Comm comm);
#else
    HDF5NativeWriter(const std::string &fileName);
#endif
    ~HDF5NativeWriter();

    void Advance();

    void CreateAndStoreScalar(std::string const &variableName, hid_t h5Type, const void *values);
    void CreateAndStoreVar(std::string const &variableName, int dimSize, hid_t h5Type,
                           const hsize_t *global_dims, const hsize_t *offsets,
                           const hsize_t *counts, const void *values);

    int m_CurrentTimeStep;
    unsigned int m_TotalTimeSteps;

private:
    void CheckWriteGroup();

    hid_t m_FilePropertyListId;
    hid_t m_FileId;
    hid_t m_GroupId;
};
#ifdef TEST_HDF5_MPI
HDF5NativeWriter::HDF5NativeWriter(const std::string &fileName, MPI_Comm comm)
#else
HDF5NativeWriter::HDF5NativeWriter(const std::string &fileName)
#endif
: m_CurrentTimeStep(0), m_TotalTimeSteps(0)
{
    m_FilePropertyListId = H5Pcreate(H5P_FILE_ACCESS);
#ifdef TEST_HDF5_MPI
    H5Pset_fapl_mpio(m_FilePropertyListId, comm, MPI_INFO_NULL);
#endif

    H5VL_ADIOS2_set(m_FilePropertyListId);

    // std::string ts0 = "/AdiosStep0";
    // stepName = "/Step" + std::to_string(ts);
    std::string ts0 = "/Step0";

    /*
     * Create a new file collectively and release property list identifier.
     */
    m_FileId = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, m_FilePropertyListId);
    if (m_FileId < 0)
    {
        throw std::runtime_error("Unable to create file: " + fileName);
    }

    m_GroupId = H5Gcreate2(m_FileId, ts0.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (m_GroupId < 0)
    {
        throw std::runtime_error("ERROR: Unable to create HDF5 group " + ts0);
    }
}

HDF5NativeWriter::~HDF5NativeWriter()
{
    if (m_FileId < 0)
    {
        return;
    }

    // write NumStep attr
    hid_t s = H5Screate(H5S_SCALAR);

    hid_t attr = H5Acreate(m_FileId, "NumSteps", H5T_NATIVE_UINT, s, H5P_DEFAULT, H5P_DEFAULT);
    unsigned int totalAdiosSteps = m_CurrentTimeStep + 1;

    if (m_GroupId < 0)
    {
        totalAdiosSteps = m_CurrentTimeStep;
    }

    H5Awrite(attr, H5T_NATIVE_UINT, &totalAdiosSteps);

    H5Sclose(s);
    H5Aclose(attr);

    // now close necessary ids
    if (m_GroupId >= 0)
    {
        H5Gclose(m_GroupId);
    }

    H5Fclose(m_FileId);
    H5Pclose(m_FilePropertyListId);

    H5VL_ADIOS2_unset();
}

void HDF5NativeWriter::CheckWriteGroup()
{
    if (m_GroupId >= 0)
    {
        return;
    }

    std::string stepName = "/Step" + std::to_string(m_CurrentTimeStep);

    m_GroupId = H5Gcreate2(m_FileId, stepName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (m_GroupId < 0)
    {
        throw std::runtime_error("ERROR: Unable to create HDF5 group " + stepName);
    }
}

void HDF5NativeWriter::CreateAndStoreScalar(std::string const &variableName, hid_t h5Type,
                                            const void *values)
{
    CheckWriteGroup();

    // write scalar
    hid_t filespaceID = H5Screate(H5S_SCALAR);
    hid_t plistID = H5Pcreate(H5P_DATASET_XFER);
#ifdef TEST_HDF5_MPI
    H5Pset_dxpl_mpio(plistID, H5FD_MPIO_COLLECTIVE);
#endif

    hid_t dsetID;

    if (h5Type != H5T_STRING)
    {
        dsetID = H5Dcreate(m_GroupId, variableName.c_str(), h5Type, filespaceID, H5P_DEFAULT,
                           H5P_DEFAULT, H5P_DEFAULT);
        herr_t status = H5Dwrite(dsetID, h5Type, H5S_ALL, H5S_ALL, plistID, values);
        EXPECT_TRUE(status > 0);
    }
    else
    {
        /* Create a datatype to refer to. */
        hid_t type = H5Tcopy(H5T_C_S1);
        char *strval = (char *)values;
        hid_t ret = H5Tset_size(type, strlen(strval));

        ret = H5Tset_strpad(type, H5T_STR_NULLTERM);

        /* Test creating a "normal" sized string attribute */
        dsetID = H5Dcreate(m_GroupId, variableName.c_str(), type, filespaceID, H5P_DEFAULT,
                           H5P_DEFAULT, H5P_DEFAULT);

        ret = H5Dwrite(dsetID, type, H5S_ALL, H5S_ALL, plistID, values);
        EXPECT_GE(ret, 0);

#ifdef DOUBLECHECK
        size_t typesize = H5Tget_size(type);
        char *val = (char *)(calloc(typesize, sizeof(char)));

        hid_t ret2 = H5Dread(dsetID, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, val);
        EXPECT_GE(ret2, 0);
        std::cerr << "        ....  typesize=" << typesize << "  val=" << val << std::endl;
        free val;
#endif
    }

    H5Sclose(filespaceID);
    H5Dclose(dsetID);
}

void HDF5NativeWriter::CreateAndStoreVar(std::string const &variableName, int dimSize, hid_t h5Type,
                                         const hsize_t *global_dims, const hsize_t *offsets,
                                         const hsize_t *counts, const void *values)
{
    if (h5Type == H5T_STRING)
    {
        throw std::runtime_error("Sync with ADIOS2. It does not store string "
                                 "var with dimensions yet!");
    }

    CheckWriteGroup();
    hid_t fileSpace = H5Screate_simple(dimSize, global_dims, NULL);

    hid_t dsetID = H5Dcreate(m_GroupId, variableName.c_str(), h5Type, fileSpace, H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
    hid_t memSpace = H5Screate_simple(dimSize, counts, NULL);

    // Select hyperslab
    fileSpace = H5Dget_space(dsetID);
    H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, offsets, NULL, counts, NULL);

    //  Create property list for collective dataset write.

    hid_t plistID = H5Pcreate(H5P_DATASET_XFER);
#ifdef TEST_HDF5_MPI
    H5Pset_dxpl_mpio(plistID, H5FD_MPIO_COLLECTIVE);
#endif
    herr_t status = H5Dwrite(dsetID, h5Type, memSpace, fileSpace, plistID, values);

    if (status < 0)
    {
        throw std::runtime_error("ERROR: HDF5 file Write failed, in call to Write\n");
    }

    H5Dclose(dsetID);
    H5Sclose(fileSpace);
    H5Sclose(memSpace);
    H5Pclose(plistID);
}

void HDF5NativeWriter::Advance()
{
    if (m_GroupId >= 0)
    {
        H5Gclose(m_GroupId);
        m_GroupId = -1;
    }
    ++m_CurrentTimeStep;
}

//
//
//
HDF5NativeReader::HDF5NativeReader(const std::string fileName)
: m_CurrentTimeStep(0), m_TotalTimeSteps(0)
{
    m_FilePropertyListId = H5Pcreate(H5P_FILE_ACCESS);
#ifdef TEST_HDF5_MPI
    // read a file collectively
    H5Pset_fapl_mpio(m_FilePropertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif

    H5VL_ADIOS2_set(m_FilePropertyListId);

    m_FileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, m_FilePropertyListId);
    if (m_FileId < 0)
    {
        throw std::runtime_error("Unable to open " + fileName + " for reading");
    }

    std::string ts0 = "/Step0";
    m_GroupId = H5Gopen(m_FileId, ts0.c_str(), H5P_DEFAULT);
    if (m_GroupId < 0)
    {
        throw std::runtime_error("Unable to open group " + ts0 + " for reading");
    }

    hid_t attrId = H5Aopen(m_FileId, "NumSteps", H5P_DEFAULT);
    if (attrId < 0)
    {
        throw std::runtime_error("Unable to open attribute NumSteps");
    }
    H5Aread(attrId, H5T_NATIVE_UINT, &m_TotalTimeSteps);
    H5Aclose(attrId);
}

HDF5NativeReader::~HDF5NativeReader()
{
    if (m_GroupId >= 0)
    {
        H5Gclose(m_GroupId);
    }

    H5Fclose(m_FileId);
    H5Pclose(m_FilePropertyListId);
    H5VL_ADIOS2_unset();
}

void HDF5NativeReader::GetVarInfo(const std::string varName, std::vector<hsize_t> &dims,
                                  hid_t &h5Type)
{
    hid_t dataSetId = H5Dopen(m_GroupId, varName.c_str(), H5P_DEFAULT);
    if (dataSetId < 0)
    {
        throw std::runtime_error("Unable to open dataset " + varName + " when getVarInfo");
    }

    hid_t fileSpaceId = H5Dget_space(dataSetId);
    if (fileSpaceId < 0)
    {
        throw std::runtime_error("Unable to get filespace for dataset " + varName);
    }

    const int ndims = H5Sget_simple_extent_ndims(fileSpaceId);
    if (ndims < 0)
    {
        throw std::runtime_error("Unable to get number of dimensions for dataset " + varName);
    }

    dims.resize(ndims);
    if (H5Sget_simple_extent_dims(fileSpaceId, dims.data(), NULL) != ndims)
    {
        throw std::runtime_error("Unable to get dimensions for dataset " + varName);
    }

    h5Type = H5Dget_type(dataSetId);

    H5Sclose(fileSpaceId);
    H5Dclose(dataSetId);
}

bool HDF5NativeReader::Advance()
{
    if (m_GroupId >= 0)
    {
        H5Gclose(m_GroupId);
        m_GroupId = -1;
    }

    if (m_CurrentTimeStep + 1 >= static_cast<int>(m_TotalTimeSteps))
    {
        return false;
    }

    const std::string tsName = "/Step" + std::to_string(m_CurrentTimeStep + 1);
    m_GroupId = H5Gopen(m_FileId, tsName.c_str(), H5P_DEFAULT);
    if (m_GroupId < 0)
    {
        throw std::runtime_error("Unable to open group " + tsName + " for reading");
    }
    ++m_CurrentTimeStep;

    return true;
}

void HDF5NativeReader::ReadString(const std::string varName, std::string &result)
{
    if (m_GroupId < 0)
    {
        throw std::runtime_error("Can't read variable " + varName +
                                 " since a group is not currently open");
    }

    hid_t dataSetId = H5Dopen(m_GroupId, varName.c_str(), H5P_DEFAULT);
    if (dataSetId < 0)
    {
        throw std::runtime_error("Unable to open dataset " + varName + "when ReadVar");
    }

    hid_t h5Type = H5Dget_type(dataSetId);
    size_t typesize = H5Tget_size(h5Type); // returns a fix number, 30
    char *val = (char *)(calloc(typesize, sizeof(char)));
    hid_t ret2 = H5Dread(dataSetId, h5Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, val);
    EXPECT_TRUE(ret2 >= 0);
    // result.assign(val, typesize);
    result.assign(val, strlen(val));
    free(val);

    H5Dclose(dataSetId);
}

void HDF5NativeReader::ReadVar(const std::string varName, void *dataArray, hsize_t *offset,
                               hsize_t *count, const size_t memspaceSize)
{
    if (m_GroupId < 0)
    {
        throw std::runtime_error("Can't read variable " + varName +
                                 " since a group is not currently open");
    }

    hid_t dataSetId = H5Dopen(m_GroupId, varName.c_str(), H5P_DEFAULT);
    if (dataSetId < 0)
    {
        throw std::runtime_error("Unable to open dataset " + varName + "when ReadVar");
    }
    hid_t fileSpace = H5Dget_space(dataSetId);
    if (fileSpace < 0)
    {
        throw std::runtime_error("Unable to get filespace for dataset " + varName);
    }

    hid_t h5type = H5Dget_type(dataSetId);

    // Extend reader to support read by hyperslab selection
    // Reference link: https://support.hdfgroup.org/HDF5/Tutor/select.html
    // Check if hyperspace is provided
    if (offset && count)
    {
        // Get the dataspace
        hid_t dataspace = H5Dget_space(dataSetId);
        // Define hyperslab in the dataset
        hid_t status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
        if (status < 0)
        {
            throw std::runtime_error("Unable to create a selection for dataset" + varName);
        }

        /*
        hsize_t dimsm[1];
        dimsm[0] = memspaceSize;
        hid_t memspace = H5Screate_simple(1, dimsm, NULL);
        */
        hid_t memspace = H5S_ALL;

        hid_t ret = H5Dread(dataSetId, h5type, memspace, dataspace, H5P_DEFAULT, dataArray);
        EXPECT_TRUE(ret >= 0);
    }
    else
    {
        hid_t ret = H5Dread(dataSetId, h5type, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataArray);
        EXPECT_TRUE(ret >= 0);
    }

    H5Sclose(fileSpace);
    H5Dclose(dataSetId);
}

//******************************************************************************
// 1D 1x8 test data
//******************************************************************************

// HDF5 API write, HDF5 API read
TEST_F(H5VolWriteReadTest, H5VolWriteHDF5Read1D8)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

#ifdef TEST_HDF5_MPI
    const std::string fname = "H5VolTest1D8_MPI.bp";
#else
    const std::string fname = "H5VolTest1D8.bp";
#endif
    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const std::size_t Nx = 8;

    // Number of steps
    const std::size_t NSteps = 3;

    {
#ifdef TEST_HDF5_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

        HDF5NativeWriter h5writer(fname, MPI_COMM_WORLD);
#else
        HDF5NativeWriter h5writer(fname);
#endif

        int dimSize = 1;
        hsize_t global_dims[1] = {Nx * mpiSize};
        hsize_t count[1] = {Nx};
        hsize_t offset[1] = {Nx * mpiRank};

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, step, mpiRank, mpiSize);

            h5writer.CreateAndStoreScalar("iString", H5T_STRING, currentTestData.S1.data());
            h5writer.CreateAndStoreVar("i8", dimSize, H5T_NATIVE_INT8, global_dims, offset, count,
                                       currentTestData.I8.data());
            h5writer.CreateAndStoreVar("i16", dimSize, H5T_NATIVE_SHORT, global_dims, offset, count,
                                       currentTestData.I16.data());
            h5writer.CreateAndStoreVar("i32", dimSize, H5T_NATIVE_INT, global_dims, offset, count,
                                       currentTestData.I32.data());
            h5writer.CreateAndStoreVar("i64", dimSize, H5T_NATIVE_LONG, global_dims, offset, count,
                                       currentTestData.I64.data());
            h5writer.CreateAndStoreVar("u8", dimSize, H5T_NATIVE_UCHAR, global_dims, offset, count,
                                       currentTestData.U8.data());
            h5writer.CreateAndStoreVar("u16", dimSize, H5T_NATIVE_USHORT, global_dims, offset,
                                       count, currentTestData.U16.data());
            h5writer.CreateAndStoreVar("u32", dimSize, H5T_NATIVE_UINT, global_dims, offset, count,
                                       currentTestData.U32.data());
            h5writer.CreateAndStoreVar("u64", dimSize, H5T_NATIVE_ULONG, global_dims, offset, count,
                                       currentTestData.U64.data());
            h5writer.CreateAndStoreVar("r32", dimSize, H5T_NATIVE_FLOAT, global_dims, offset, count,
                                       currentTestData.R32.data());
            h5writer.CreateAndStoreVar("r64", dimSize, H5T_NATIVE_DOUBLE, global_dims, offset,
                                       count, currentTestData.R64.data());
            h5writer.Advance();
        }
    }

    {
        const size_t arraySize = Nx;
        std::string IString;
        std::array<int8_t, arraySize> I8;
        std::array<int16_t, arraySize> I16;
        std::array<int32_t, arraySize> I32;
        std::array<int64_t, arraySize> I64;
        std::array<uint8_t, arraySize> U8;
        std::array<uint16_t, arraySize> U16;
        std::array<uint32_t, arraySize> U32;
        std::array<uint64_t, arraySize> U64;
        std::array<float, arraySize> R32;
        std::array<double, arraySize> R64;

        HDF5NativeReader hdf5Reader(fname);
        // 1D
        hsize_t count[1], offset[1];
        offset[0] = mpiRank * Nx;
        count[0] = Nx;
        size_t globalArraySize = Nx * mpiSize;

        // For each variable, we would verify its global size and type.
        // Then we would retrieve the data back which is written by the
        // current process and validate the value
        for (size_t t = 0; t < NSteps; ++t)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, t, mpiRank, mpiSize);

            std::vector<hsize_t> gDims;
            hid_t h5Type;

            hdf5Reader.GetVarInfo("iString", gDims, h5Type);
            ASSERT_EQ(gDims.size(), 0);
            hdf5Reader.ReadString("iString", IString);

            hdf5Reader.GetVarInfo("i8", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_INT8), 1);
            ASSERT_EQ(gDims.size(), 1);
            ASSERT_EQ(gDims[0], globalArraySize);
            hdf5Reader.ReadVar("i8", I8.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("i16", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_SHORT), 1);
            ASSERT_EQ(gDims.size(), 1);
            ASSERT_EQ(gDims[0], globalArraySize);
            hdf5Reader.ReadVar("i16", I16.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("i32", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_INT), 1);
            ASSERT_EQ(gDims.size(), 1);
            ASSERT_EQ(gDims[0], globalArraySize);
            hdf5Reader.ReadVar("i32", I32.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("i64", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_LONG), 1);
            ASSERT_EQ(gDims.size(), 1);
            ASSERT_EQ(gDims[0], globalArraySize);
            hdf5Reader.ReadVar("i64", I64.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("u8", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_UCHAR), 1);
            ASSERT_EQ(gDims.size(), 1);
            ASSERT_EQ(gDims[0], globalArraySize);
            hdf5Reader.ReadVar("u8", U8.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("u16", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_USHORT), 1);
            ASSERT_EQ(gDims.size(), 1);
            ASSERT_EQ(gDims[0], globalArraySize);
            hdf5Reader.ReadVar("u16", U16.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("u32", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_UINT), 1);
            ASSERT_EQ(gDims.size(), 1);
            ASSERT_EQ(gDims[0], globalArraySize);
            hdf5Reader.ReadVar("u32", U32.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("u64", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_ULONG), 1);
            ASSERT_EQ(gDims.size(), 1);
            ASSERT_EQ(gDims[0], globalArraySize);
            hdf5Reader.ReadVar("u64", U64.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("r32", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_FLOAT), 1);
            ASSERT_EQ(gDims.size(), 1);
            ASSERT_EQ(gDims[0], globalArraySize);
            hdf5Reader.ReadVar("r32", R32.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("r64", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_DOUBLE), 1);
            ASSERT_EQ(gDims.size(), 1);
            ASSERT_EQ(gDims[0], globalArraySize);
            hdf5Reader.ReadVar("r64", R64.data(), offset, count, arraySize);

            // #EXPECT_EQ(IString, currentTestData.S1);

            // Check if it's correct
            for (size_t i = 0; i < Nx; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                EXPECT_EQ(I8[i], currentTestData.I8[i]) << msg;
                EXPECT_EQ(I16[i], currentTestData.I16[i]) << msg;
                EXPECT_EQ(I32[i], currentTestData.I32[i]) << msg;
                EXPECT_EQ(I64[i], currentTestData.I64[i]) << msg;
                EXPECT_EQ(U8[i], currentTestData.U8[i]) << msg;
                EXPECT_EQ(U16[i], currentTestData.U16[i]) << msg;
                EXPECT_EQ(U32[i], currentTestData.U32[i]) << msg;
                EXPECT_EQ(U64[i], currentTestData.U64[i]) << msg;
                EXPECT_EQ(R32[i], currentTestData.R32[i]) << msg;
                EXPECT_EQ(R64[i], currentTestData.R64[i]) << msg;
            }
            hdf5Reader.Advance();
        }
    }
}

// Native HDF5 write, ADIOS2 read

//******************************************************************************
// 2D 2x4 test data
//******************************************************************************

// HDF5 API write, HDF5 API read
TEST_F(H5VolWriteReadTest, H5VolWriteHDF5Read2D2x4)
{
    // Each process would write a 2x4 array and all processes would
    // form a 2D 2 * (numberOfProcess*Nx) matrix where Nx is 4 here
#ifdef TEST_HDF5_MPI
    std::string fname = "H5VolTest2D2x4_MPI.bp";
#else
    std::string fname = "H5VolTest2D2x4.bp";
#endif
    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const std::size_t Nx = 4;

    // Number of rows
    const std::size_t Ny = 2;

    // Number of steps
    const std::size_t NSteps = 3;

    {
#ifdef TEST_HDF5_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

        HDF5NativeWriter h5writer(fname, MPI_COMM_WORLD);
#else
        HDF5NativeWriter h5writer(fname);
#endif

        int dimSize = 2;
        hsize_t global_dims[2] = {Ny, Nx * mpiSize};
        hsize_t count[2] = {Ny, Nx};
        hsize_t offset[2] = {0, Nx * mpiRank};

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, step, mpiRank, mpiSize);

            h5writer.CreateAndStoreScalar("iString", H5T_STRING, currentTestData.S1.data());
            h5writer.CreateAndStoreVar("i8", dimSize, H5T_NATIVE_INT8, global_dims, offset, count,
                                       currentTestData.I8.data());
            h5writer.CreateAndStoreVar("i16", dimSize, H5T_NATIVE_SHORT, global_dims, offset, count,
                                       currentTestData.I16.data());
            h5writer.CreateAndStoreVar("i32", dimSize, H5T_NATIVE_INT, global_dims, offset, count,
                                       currentTestData.I32.data());
            h5writer.CreateAndStoreVar("i64", dimSize, H5T_NATIVE_LONG, global_dims, offset, count,
                                       currentTestData.I64.data());
            h5writer.CreateAndStoreVar("u8", dimSize, H5T_NATIVE_UCHAR, global_dims, offset, count,
                                       currentTestData.U8.data());
            h5writer.CreateAndStoreVar("u16", dimSize, H5T_NATIVE_USHORT, global_dims, offset,
                                       count, currentTestData.U16.data());
            h5writer.CreateAndStoreVar("u32", dimSize, H5T_NATIVE_UINT, global_dims, offset, count,
                                       currentTestData.U32.data());
            h5writer.CreateAndStoreVar("u64", dimSize, H5T_NATIVE_ULONG, global_dims, offset, count,
                                       currentTestData.U64.data());
            h5writer.CreateAndStoreVar("r32", dimSize, H5T_NATIVE_FLOAT, global_dims, offset, count,
                                       currentTestData.R32.data());
            h5writer.CreateAndStoreVar("r64", dimSize, H5T_NATIVE_DOUBLE, global_dims, offset,
                                       count, currentTestData.R64.data());
            h5writer.Advance();
        }
    }

    {
        HDF5NativeReader hdf5Reader(fname);

        std::string IString;
        const size_t arraySize = Nx * Ny;
        std::array<int8_t, arraySize> I8;
        std::array<int16_t, arraySize> I16;
        std::array<int32_t, arraySize> I32;
        std::array<int64_t, arraySize> I64;
        std::array<uint8_t, arraySize> U8;
        std::array<uint16_t, arraySize> U16;
        std::array<uint32_t, arraySize> U32;
        std::array<uint64_t, arraySize> U64;
        std::array<float, arraySize> R32;
        std::array<double, arraySize> R64;
        // 2D
        hsize_t count[2], offset[2];

        offset[0] = 0;
        offset[1] = mpiRank * Nx;
        count[0] = Ny;
        count[1] = Nx;

        size_t globalArraySize = Nx * mpiSize;

        // For each variable, we would verify its global size and type.
        // Then we would retrieve the data back which is written by the
        // current process and validate the value
        for (size_t t = 0; t < NSteps; ++t)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, t, mpiRank, mpiSize);

            std::vector<hsize_t> gDims;
            hid_t h5Type;

            hdf5Reader.GetVarInfo("iString", gDims, h5Type);
            // ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_IN), 1);
            ASSERT_EQ(gDims.size(), 0);
            hdf5Reader.ReadString("iString", IString);

            hdf5Reader.GetVarInfo("i8", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_INT8), 1);
            ASSERT_EQ(gDims.size(), 2);
            ASSERT_EQ(gDims[0], 2);
            ASSERT_EQ(gDims[1], globalArraySize);
            hdf5Reader.ReadVar("i8", I8.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("i16", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_SHORT), 1);
            ASSERT_EQ(gDims.size(), 2);
            ASSERT_EQ(gDims[0], 2);
            ASSERT_EQ(gDims[1], globalArraySize);
            hdf5Reader.ReadVar("i16", I16.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("i32", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_INT), 1);
            ASSERT_EQ(gDims.size(), 2);
            ASSERT_EQ(gDims[0], 2);
            ASSERT_EQ(gDims[1], globalArraySize);
            hdf5Reader.ReadVar("i32", I32.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("i64", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_LONG), 1);
            ASSERT_EQ(gDims.size(), 2);
            ASSERT_EQ(gDims[0], 2);
            ASSERT_EQ(gDims[1], globalArraySize);
            hdf5Reader.ReadVar("i64", I64.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("u8", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_UCHAR), 1);
            ASSERT_EQ(gDims.size(), 2);
            ASSERT_EQ(gDims[0], 2);
            ASSERT_EQ(gDims[1], globalArraySize);
            hdf5Reader.ReadVar("u8", U8.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("u16", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_USHORT), 1);
            ASSERT_EQ(gDims.size(), 2);
            ASSERT_EQ(gDims[0], 2);
            ASSERT_EQ(gDims[1], globalArraySize);
            hdf5Reader.ReadVar("u16", U16.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("u32", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_UINT), 1);
            ASSERT_EQ(gDims.size(), 2);
            ASSERT_EQ(gDims[0], 2);
            ASSERT_EQ(gDims[1], globalArraySize);
            hdf5Reader.ReadVar("u32", U32.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("u64", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_ULONG), 1);
            ASSERT_EQ(gDims.size(), 2);
            ASSERT_EQ(gDims[0], 2);
            ASSERT_EQ(gDims[1], globalArraySize);
            hdf5Reader.ReadVar("u64", U64.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("r32", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_FLOAT), 1);
            ASSERT_EQ(gDims.size(), 2);
            ASSERT_EQ(gDims[0], 2);
            ASSERT_EQ(gDims[1], globalArraySize);
            hdf5Reader.ReadVar("r32", R32.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("r64", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_DOUBLE), 1);
            ASSERT_EQ(gDims.size(), 2);
            ASSERT_EQ(gDims[0], 2);
            ASSERT_EQ(gDims[1], globalArraySize);
            hdf5Reader.ReadVar("r64", R64.data(), offset, count, arraySize);

            EXPECT_EQ(IString, currentTestData.S1);

            // Check if it's correct
            for (size_t i = 0; i < Nx; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                EXPECT_EQ(I8[i], currentTestData.I8[i]) << msg;
                EXPECT_EQ(I16[i], currentTestData.I16[i]) << msg;
                EXPECT_EQ(I32[i], currentTestData.I32[i]) << msg;
                EXPECT_EQ(I64[i], currentTestData.I64[i]) << msg;
                EXPECT_EQ(U8[i], currentTestData.U8[i]) << msg;
                EXPECT_EQ(U16[i], currentTestData.U16[i]) << msg;
                EXPECT_EQ(U32[i], currentTestData.U32[i]) << msg;
                EXPECT_EQ(U64[i], currentTestData.U64[i]) << msg;
                EXPECT_EQ(R32[i], currentTestData.R32[i]) << msg;
                EXPECT_EQ(R64[i], currentTestData.R64[i]) << msg;
            }
            hdf5Reader.Advance();
        }
    }
}

//******************************************************************************
// 2D 4x2 test data
//******************************************************************************

// HDF5 API write, HDF5 API read
TEST_F(H5VolWriteReadTest, H5VolWriteHDF5Read2D4x2)
{

    // Each process would write a 4x2 array and all processes would
    // form a 2D 4 * (NumberOfProcess * Nx) matrix where Nx is 2 here
    std::string fname = "H5VolTest2D4x2.bp";

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const std::size_t Nx = 2;
    // Number of cols
    const std::size_t Ny = 4;

    // Number of steps
    const std::size_t NSteps = 3;

    {
#ifdef TEST_HDF5_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

        HDF5NativeWriter h5writer(fname, MPI_COMM_WORLD);
#else
        HDF5NativeWriter h5writer(fname);
#endif

        int dimSize = 2;
        hsize_t global_dims[2] = {Ny, Nx * mpiSize};
        hsize_t count[2] = {Ny, Nx};
        hsize_t offset[2] = {0, Nx * mpiRank};

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, step, mpiRank, mpiSize);

            h5writer.CreateAndStoreScalar("iString", H5T_STRING, currentTestData.S1.data());
            h5writer.CreateAndStoreVar("i8", dimSize, H5T_NATIVE_INT8, global_dims, offset, count,
                                       currentTestData.I8.data());
            h5writer.CreateAndStoreVar("i16", dimSize, H5T_NATIVE_SHORT, global_dims, offset, count,
                                       currentTestData.I16.data());
            h5writer.CreateAndStoreVar("i32", dimSize, H5T_NATIVE_INT, global_dims, offset, count,
                                       currentTestData.I32.data());
            h5writer.CreateAndStoreVar("i64", dimSize, H5T_NATIVE_LONG, global_dims, offset, count,
                                       currentTestData.I64.data());
            h5writer.CreateAndStoreVar("u8", dimSize, H5T_NATIVE_UCHAR, global_dims, offset, count,
                                       currentTestData.U8.data());
            h5writer.CreateAndStoreVar("u16", dimSize, H5T_NATIVE_USHORT, global_dims, offset,
                                       count, currentTestData.U16.data());
            h5writer.CreateAndStoreVar("u32", dimSize, H5T_NATIVE_UINT, global_dims, offset, count,
                                       currentTestData.U32.data());
            h5writer.CreateAndStoreVar("u64", dimSize, H5T_NATIVE_ULONG, global_dims, offset, count,
                                       currentTestData.U64.data());
            h5writer.CreateAndStoreVar("r32", dimSize, H5T_NATIVE_FLOAT, global_dims, offset, count,
                                       currentTestData.R32.data());
            h5writer.CreateAndStoreVar("r64", dimSize, H5T_NATIVE_DOUBLE, global_dims, offset,
                                       count, currentTestData.R64.data());
            h5writer.Advance();
        }
    }

    {

        HDF5NativeReader hdf5Reader(fname);

        std::string IString;
        const size_t arraySize = Nx * Ny;
        std::array<int8_t, arraySize> I8;
        std::array<int16_t, arraySize> I16;
        std::array<int32_t, arraySize> I32;
        std::array<int64_t, arraySize> I64;
        std::array<uint8_t, arraySize> U8;
        std::array<uint16_t, arraySize> U16;
        std::array<uint32_t, arraySize> U32;
        std::array<uint64_t, arraySize> U64;
        std::array<float, arraySize> R32;
        std::array<double, arraySize> R64;
        // 2D
        hsize_t count[2], offset[2];
        offset[0] = 0;
        offset[1] = mpiRank * Nx;
        count[0] = Ny;
        count[1] = Nx;
        size_t globalArraySize = Nx * mpiSize;

        // For each variable, we would verify its global size and type.
        // Then we would retrieve the data back which is written by the
        // current process and validate the value
        for (size_t t = 0; t < NSteps; ++t)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, t, mpiRank, mpiSize);

            std::vector<hsize_t> gDims;
            hid_t h5Type;

            hdf5Reader.GetVarInfo("iString", gDims, h5Type);
            // ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_IN), 1);
            ASSERT_EQ(gDims.size(), 0);
            hdf5Reader.ReadString("iString", IString);

            hdf5Reader.GetVarInfo("i8", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_INT8), 1);
            ASSERT_EQ(gDims.size(), 2);
            ASSERT_EQ(gDims[0], 4);
            ASSERT_EQ(gDims[1], globalArraySize);
            hdf5Reader.ReadVar("i8", I8.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("i16", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_SHORT), 1);
            ASSERT_EQ(gDims.size(), 2);
            ASSERT_EQ(gDims[0], 4);
            ASSERT_EQ(gDims[1], globalArraySize);
            hdf5Reader.ReadVar("i16", I16.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("i32", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_INT), 1);
            ASSERT_EQ(gDims.size(), 2);
            ASSERT_EQ(gDims[0], 4);
            ASSERT_EQ(gDims[1], globalArraySize);
            hdf5Reader.ReadVar("i32", I32.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("i64", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_LONG), 1);
            ASSERT_EQ(gDims.size(), 2);
            ASSERT_EQ(gDims[0], 4);
            ASSERT_EQ(gDims[1], globalArraySize);
            hdf5Reader.ReadVar("i64", I64.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("u8", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_UCHAR), 1);
            ASSERT_EQ(gDims.size(), 2);
            ASSERT_EQ(gDims[0], 4);
            ASSERT_EQ(gDims[1], globalArraySize);
            hdf5Reader.ReadVar("u8", U8.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("u16", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_USHORT), 1);
            ASSERT_EQ(gDims.size(), 2);
            ASSERT_EQ(gDims[0], 4);
            ASSERT_EQ(gDims[1], globalArraySize);
            hdf5Reader.ReadVar("u16", U16.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("u32", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_UINT), 1);
            ASSERT_EQ(gDims.size(), 2);
            ASSERT_EQ(gDims[0], 4);
            ASSERT_EQ(gDims[1], globalArraySize);
            hdf5Reader.ReadVar("u32", U32.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("u64", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_ULONG), 1);
            ASSERT_EQ(gDims.size(), 2);
            ASSERT_EQ(gDims[0], 4);
            ASSERT_EQ(gDims[1], globalArraySize);
            hdf5Reader.ReadVar("u64", U64.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("r32", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_FLOAT), 1);
            ASSERT_EQ(gDims.size(), 2);
            ASSERT_EQ(gDims[0], 4);
            ASSERT_EQ(gDims[1], globalArraySize);
            hdf5Reader.ReadVar("r32", R32.data(), offset, count, arraySize);

            hdf5Reader.GetVarInfo("r64", gDims, h5Type);
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_DOUBLE), 1);
            ASSERT_EQ(gDims.size(), 2);
            ASSERT_EQ(gDims[0], 4);
            ASSERT_EQ(gDims[1], globalArraySize);
            hdf5Reader.ReadVar("r64", R64.data(), offset, count, arraySize);

            EXPECT_EQ(IString, currentTestData.S1);
            // Check if it's correct
            for (size_t i = 0; i < Nx; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                EXPECT_EQ(I8[i], currentTestData.I8[i]) << msg;
                EXPECT_EQ(I16[i], currentTestData.I16[i]) << msg;
                EXPECT_EQ(I32[i], currentTestData.I32[i]) << msg;
                EXPECT_EQ(I64[i], currentTestData.I64[i]) << msg;
                EXPECT_EQ(U8[i], currentTestData.U8[i]) << msg;
                EXPECT_EQ(U16[i], currentTestData.U16[i]) << msg;
                EXPECT_EQ(U32[i], currentTestData.U32[i]) << msg;
                EXPECT_EQ(U64[i], currentTestData.U64[i]) << msg;
                EXPECT_EQ(R32[i], currentTestData.R32[i]) << msg;
                EXPECT_EQ(R64[i], currentTestData.R64[i]) << msg;
            }
            hdf5Reader.Advance();
        }
    }
}

//******************************************************************************
// main
//******************************************************************************

int main(int argc, char **argv)
{
#ifdef TEST_HDF5_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(nullptr, nullptr, MPI_THREAD_MULTIPLE, &provided);
#endif

    int result;
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();

#ifdef TEST_HDF5_MPI
    MPI_Finalize();
#endif

    return result;
}

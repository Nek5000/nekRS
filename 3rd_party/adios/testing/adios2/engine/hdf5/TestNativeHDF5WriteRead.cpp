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

#include "../SmallTestData.h"

class HDF5WriteReadTest : public ::testing::Test
{
public:
    HDF5WriteReadTest() = default;

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

    /*
      void WriteVar(const std::string varName, void *dataArray,
                    hsize_t *offset = nullptr, hsize_t *count = nullptr,
                    const size_t memsspaceSize = 0);
    */
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
        EXPECT_TRUE(status >= 0);
    }
    else
    {
        /* Create a datatype to refer to. */
        hid_t type = H5Tcopy(H5T_C_S1);
        char *strval = (char *)values;
        hid_t ret = H5Tset_size(type, strlen(strval));
        EXPECT_TRUE(ret >= 0);
        ret = H5Tset_strpad(type, H5T_STR_NULLTERM);
        EXPECT_TRUE(ret >= 0);
        /* Test creating a "normal" sized string attribute */
        dsetID = H5Dcreate(m_GroupId, variableName.c_str(), type, filespaceID, H5P_DEFAULT,
                           H5P_DEFAULT, H5P_DEFAULT);

        ret = H5Dwrite(dsetID, type, H5S_ALL, H5S_ALL, plistID, values);
        EXPECT_TRUE(ret >= 0);
#ifdef DOUBLECHECK
        size_t typesize = H5Tget_size(type);
        char *val = (char *)(calloc(typesize, sizeof(char)));

        hid_t ret2 = H5Dread(dsetID, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, val);
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

    const std::string tsName = "Step" + std::to_string(m_CurrentTimeStep + 1);
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
    size_t typesize = H5Tget_size(h5Type);

    char *val = (char *)(calloc(typesize, sizeof(char)));
    hid_t ret2 = H5Dread(dataSetId, h5Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, val);
    EXPECT_TRUE(ret2 >= 0);
    result.assign(val, typesize);
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

        hsize_t dimsm[1];
        dimsm[0] = memspaceSize;
        hid_t memspace = H5Screate_simple(1, dimsm, NULL);

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

// ADIOS2 write, native HDF5 read
TEST_F(HDF5WriteReadTest, ADIOS2HDF5WriteHDF5Read1D8)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    const std::string fname = "ADIOS2HDF5WriteHDF5Read1D8.h5";

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const std::size_t Nx = 8;

    // Number of steps
    const std::size_t NSteps = 3;

#ifdef TEST_HDF5_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // Write test data using ADIOS2

#ifdef TEST_HDF5_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    adios2::IO io = adios.DeclareIO("TestIO");

    // Declare 1D variables (NumOfProcesses * Nx)
    // The local process' part (start, count) can be defined now or later
    // before Write().
    {
        adios2::Dims shape{static_cast<unsigned int>(Nx * mpiSize)};
        adios2::Dims start{static_cast<unsigned int>(Nx * mpiRank)};
        adios2::Dims count{static_cast<unsigned int>(Nx)};

        io.DefineVariable<std::string>("iString");
        io.DefineVariable<int8_t>("i8", shape, start, count);
        io.DefineVariable<int16_t>("i16", shape, start, count);
        io.DefineVariable<int32_t>("i32", shape, start, count);
        io.DefineVariable<int64_t>("i64", shape, start, count);
        io.DefineVariable<uint8_t>("u8", shape, start, count);
        io.DefineVariable<uint16_t>("u16", shape, start, count);
        io.DefineVariable<uint32_t>("u32", shape, start, count);
        io.DefineVariable<uint64_t>("u64", shape, start, count);
        io.DefineVariable<float>("r32", shape, start, count);
        io.DefineVariable<double>("r64", shape, start, count);
    }

    // Create the HDF5 Engine
    io.SetEngine("HDF5");

    // HDf5 engine calls the HDF5 common object that calls the hDF5 library.
    // The IO functionality, SetParameters and AddTransports will be added
    // in the future. For now `io.AddTransport("file", {
    // "library", "MPI"}});` is omitted.
    // })
    // io.AddTransport("File");

    adios2::Engine engine = io.Open(fname, adios2::Mode::Write);

    for (size_t step = 0; step < NSteps; ++step)
    {
        // Generate test data for each process uniquely
        SmallTestData currentTestData =
            generateNewSmallTestData(m_TestData, step, mpiRank, mpiSize);

        // Retrieve the variables that previously went out of scope
        auto var_iString = io.InquireVariable<std::string>("iString");
        auto var_i8 = io.InquireVariable<int8_t>("i8");
        auto var_i16 = io.InquireVariable<int16_t>("i16");
        auto var_i32 = io.InquireVariable<int32_t>("i32");
        auto var_i64 = io.InquireVariable<int64_t>("i64");
        auto var_u8 = io.InquireVariable<uint8_t>("u8");
        auto var_u16 = io.InquireVariable<uint16_t>("u16");
        auto var_u32 = io.InquireVariable<uint32_t>("u32");
        auto var_u64 = io.InquireVariable<uint64_t>("u64");
        auto var_r32 = io.InquireVariable<float>("r32");
        auto var_r64 = io.InquireVariable<double>("r64");

        // Make a 1D selection to describe the local dimensions of the
        // variable we write and its offsets in the global spaces
        adios2::Box<adios2::Dims> sel({mpiRank * Nx}, {Nx});
        var_i8.SetSelection(sel);
        var_i16.SetSelection(sel);
        var_i32.SetSelection(sel);
        var_i64.SetSelection(sel);
        var_u8.SetSelection(sel);
        var_u16.SetSelection(sel);
        var_u32.SetSelection(sel);
        var_u64.SetSelection(sel);
        var_r32.SetSelection(sel);
        var_r64.SetSelection(sel);

        // Write each one
        // fill in the variable with values from starting index to
        // starting index + count
        engine.BeginStep();
        engine.Put(var_iString, currentTestData.S1);
        engine.Put(var_i8, currentTestData.I8.data());
        engine.Put(var_i16, currentTestData.I16.data());
        engine.Put(var_i32, currentTestData.I32.data());
        engine.Put(var_i64, currentTestData.I64.data());
        engine.Put(var_u8, currentTestData.U8.data());
        engine.Put(var_u16, currentTestData.U16.data());
        engine.Put(var_u32, currentTestData.U32.data());
        engine.Put(var_u64, currentTestData.U64.data());
        engine.Put(var_r32, currentTestData.R32.data());
        engine.Put(var_r64, currentTestData.R64.data());
        // Advance to the next time step
        engine.EndStep();
    }

    // Close the file
    engine.Close();

    bool doRead = true;
    if (doRead)
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

            // auto var_iString = io.InquireVariable<std::string>("iString");
            hdf5Reader.GetVarInfo("iString", gDims, h5Type);
            // ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_IN), 1);
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
#ifdef _WIN32
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_LLONG), 1);
#else
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_LONG), 1);
#endif
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
#ifdef _WIN32
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_ULLONG), 1);
#else
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_ULONG), 1);
#endif
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

// ADIOS2 write, ADIOS2 read
TEST_F(HDF5WriteReadTest, ADIOS2HDF5WriteADIOS2HDF5Read1D8)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    const std::string fname = "ADIOS2HDF5WriteADIOS2HDF5Read1D8.h5";

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const std::size_t Nx = 8;

    // Number of steps
    const std::size_t NSteps = 3;

#ifdef TEST_HDF5_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // Write test data using ADIOS2

#ifdef TEST_HDF5_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    adios2::IO io = adios.DeclareIO("TestIO");

    // Declare 1D variables (NumOfProcesses * Nx)
    // The local process' part (start, count) can be defined now or later
    // before Write().
    {
        adios2::Dims shape{static_cast<unsigned int>(Nx * mpiSize)};
        adios2::Dims start{static_cast<unsigned int>(Nx * mpiRank)};
        adios2::Dims count{static_cast<unsigned int>(Nx)};

        io.DefineVariable<std::string>("iString");
        io.DefineVariable<int8_t>("i8", shape, start, count);
        io.DefineVariable<int16_t>("i16", shape, start, count);
        io.DefineVariable<int32_t>("i32", shape, start, count);
        io.DefineVariable<int64_t>("i64", shape, start, count);
        io.DefineVariable<uint8_t>("u8", shape, start, count);
        io.DefineVariable<uint16_t>("u16", shape, start, count);
        io.DefineVariable<uint32_t>("u32", shape, start, count);
        io.DefineVariable<uint64_t>("u64", shape, start, count);
        io.DefineVariable<float>("r32", shape, start, count);
        io.DefineVariable<double>("r64", shape, start, count);
    }

    // Create the HDF5 Engine
    io.SetEngine("HDF5");

    // HDf5 engine calls the HDF5 common object that calls the hDF5 library.
    // The IO functionality, SetParameters and AddTransports will be added
    // in the future. For now `io.AddTransport("file", {
    // "library", "MPI"}});` is omitted.
    // })
    // io.AddTransport("File");

    adios2::Engine engine = io.Open(fname, adios2::Mode::Write);

    for (size_t step = 0; step < NSteps; ++step)
    {
        // Generate test data for each process uniquely
        SmallTestData currentTestData =
            generateNewSmallTestData(m_TestData, step, mpiRank, mpiSize);

        // Retrieve the variables that previously went out of scope
        auto var_iString = io.InquireVariable<std::string>("iString");
        auto var_i8 = io.InquireVariable<int8_t>("i8");
        auto var_i16 = io.InquireVariable<int16_t>("i16");
        auto var_i32 = io.InquireVariable<int32_t>("i32");
        auto var_i64 = io.InquireVariable<int64_t>("i64");
        auto var_u8 = io.InquireVariable<uint8_t>("u8");
        auto var_u16 = io.InquireVariable<uint16_t>("u16");
        auto var_u32 = io.InquireVariable<uint32_t>("u32");
        auto var_u64 = io.InquireVariable<uint64_t>("u64");
        auto var_r32 = io.InquireVariable<float>("r32");
        auto var_r64 = io.InquireVariable<double>("r64");

        // Make a 1D selection to describe the local dimensions of the
        // variable we write and its offsets in the global spaces

        adios2::Box<adios2::Dims> sel({mpiRank * Nx}, {Nx});

        EXPECT_THROW(var_iString.SetSelection(sel), std::invalid_argument);

        var_i8.SetSelection(sel);
        var_i16.SetSelection(sel);
        var_i32.SetSelection(sel);
        var_i64.SetSelection(sel);
        var_u8.SetSelection(sel);
        var_u16.SetSelection(sel);
        var_u32.SetSelection(sel);
        var_u64.SetSelection(sel);
        var_r32.SetSelection(sel);
        var_r64.SetSelection(sel);

        // Write each one
        // fill in the variable with values from starting index to
        // starting index + count
        engine.BeginStep();
        engine.Put(var_iString, currentTestData.S1);
        engine.Put(var_i8, currentTestData.I8.data());
        engine.Put(var_i16, currentTestData.I16.data());
        engine.Put(var_i32, currentTestData.I32.data());
        engine.Put(var_i64, currentTestData.I64.data());
        engine.Put(var_u8, currentTestData.U8.data());
        engine.Put(var_u16, currentTestData.U16.data());
        engine.Put(var_u32, currentTestData.U32.data());
        engine.Put(var_u64, currentTestData.U64.data());
        engine.Put(var_r32, currentTestData.R32.data());
        engine.Put(var_r64, currentTestData.R64.data());
        // Advance to the next time step
        engine.EndStep();
    }

    // Close the file
    engine.Close();

    {
        adios2::IO io = adios.DeclareIO("HDF5ReadIO");
        io.SetEngine("HDF5");

        adios2::Engine hdf5Reader = io.Open(fname, adios2::Mode::Read);

        auto var_iString = io.InquireVariable<std::string>("iString");
        EXPECT_TRUE(var_iString);
        ASSERT_EQ(var_iString.Shape().size(), 0);
        ASSERT_EQ(var_iString.Steps(), NSteps);

        auto var_i8 = io.InquireVariable<int8_t>("i8");
        EXPECT_TRUE(var_i8);
        ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i8.Steps(), NSteps);
        ASSERT_EQ(var_i8.Shape()[0], mpiSize * Nx);

        auto var_i16 = io.InquireVariable<int16_t>("i16");
        EXPECT_TRUE(var_i16);
        ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i16.Steps(), NSteps);
        ASSERT_EQ(var_i16.Shape()[0], mpiSize * Nx);

        auto var_i32 = io.InquireVariable<int32_t>("i32");
        EXPECT_TRUE(var_i32);
        ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i32.Steps(), NSteps);
        ASSERT_EQ(var_i32.Shape()[0], mpiSize * Nx);

        auto var_i64 = io.InquireVariable<int64_t>("i64");
        EXPECT_TRUE(var_i64);
        ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i64.Steps(), NSteps);
        ASSERT_EQ(var_i64.Shape()[0], mpiSize * Nx);

        auto var_u8 = io.InquireVariable<uint8_t>("u8");
        EXPECT_TRUE(var_u8);
        ASSERT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u8.Steps(), NSteps);
        ASSERT_EQ(var_u8.Shape()[0], mpiSize * Nx);

        auto var_u16 = io.InquireVariable<uint16_t>("u16");
        EXPECT_TRUE(var_u16);
        ASSERT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u16.Steps(), NSteps);
        ASSERT_EQ(var_u16.Shape()[0], mpiSize * Nx);

        auto var_u32 = io.InquireVariable<uint32_t>("u32");
        EXPECT_TRUE(var_u32);
        ASSERT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u32.Steps(), NSteps);
        ASSERT_EQ(var_u32.Shape()[0], mpiSize * Nx);

        auto var_u64 = io.InquireVariable<uint64_t>("u64");
        EXPECT_TRUE(var_u64);
        ASSERT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u64.Steps(), NSteps);
        ASSERT_EQ(var_u64.Shape()[0], mpiSize * Nx);

        auto var_r32 = io.InquireVariable<float>("r32");
        EXPECT_TRUE(var_r32);
        ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r32.Steps(), NSteps);
        ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);

        auto var_r64 = io.InquireVariable<double>("r64");
        EXPECT_TRUE(var_r64);
        ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r64.Steps(), NSteps);
        ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);

        // TODO: other types

        SmallTestData testData;

        std::string IString;
        std::array<int8_t, Nx> I8;
        std::array<int16_t, Nx> I16;
        std::array<int32_t, Nx> I32;
        std::array<int64_t, Nx> I64;
        std::array<uint8_t, Nx> U8;
        std::array<uint16_t, Nx> U16;
        std::array<uint32_t, Nx> U32;
        std::array<uint64_t, Nx> U64;
        std::array<float, Nx> R32;
        std::array<double, Nx> R64;

        const adios2::Dims start{mpiRank * Nx};
        const adios2::Dims count{Nx};

        const adios2::Box<adios2::Dims> sel(start, count);

        var_i8.SetSelection(sel);
        var_i16.SetSelection(sel);
        var_i32.SetSelection(sel);
        var_i64.SetSelection(sel);

        var_u8.SetSelection(sel);
        var_u16.SetSelection(sel);
        var_u32.SetSelection(sel);
        var_u64.SetSelection(sel);

        var_r32.SetSelection(sel);
        var_r64.SetSelection(sel);

        for (size_t t = 0; t < NSteps; ++t)
        {
            var_i8.SetStepSelection({t, 1});
            var_i16.SetStepSelection({t, 1});
            var_i32.SetStepSelection({t, 1});
            var_i64.SetStepSelection({t, 1});

            var_u8.SetStepSelection({t, 1});
            var_u16.SetStepSelection({t, 1});
            var_u32.SetStepSelection({t, 1});
            var_u64.SetStepSelection({t, 1});

            var_r32.SetStepSelection({t, 1});
            var_r64.SetStepSelection({t, 1});

            // Generate test data for each rank uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(t), mpiRank, mpiSize);

            hdf5Reader.Get(var_iString, IString);

            hdf5Reader.Get(var_i8, I8.data());
            hdf5Reader.Get(var_i16, I16.data());
            hdf5Reader.Get(var_i32, I32.data());
            hdf5Reader.Get(var_i64, I64.data());

            hdf5Reader.Get(var_u8, U8.data());
            hdf5Reader.Get(var_u16, U16.data());
            hdf5Reader.Get(var_u32, U32.data());
            hdf5Reader.Get(var_u64, U64.data());

            hdf5Reader.Get(var_r32, R32.data());
            hdf5Reader.Get(var_r64, R64.data());
            hdf5Reader.PerformGets();

            EXPECT_EQ(IString, currentTestData.S1);

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
        }
        hdf5Reader.Close();
    }
}

// Native HDF5 write, ADIOS2 read
TEST_F(HDF5WriteReadTest, HDF5WriteADIOS2HDF5Read1D8)
{
    std::string fname = "HDF5WriteADIOS2HDF5Read1D8.h5";

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
            h5writer.CreateAndStoreVar("ch", dimSize, H5T_NATIVE_UCHAR, global_dims, offset, count,
                                       currentTestData.CHAR.data());
            h5writer.CreateAndStoreVar("i8", dimSize, H5T_NATIVE_INT8, global_dims, offset, count,
                                       currentTestData.I8.data());
            h5writer.CreateAndStoreVar("i16", dimSize, H5T_NATIVE_SHORT, global_dims, offset, count,
                                       currentTestData.I16.data());
            h5writer.CreateAndStoreVar("i32", dimSize, H5T_NATIVE_INT, global_dims, offset, count,
                                       currentTestData.I32.data());
#ifdef _WIN32
            h5writer.CreateAndStoreVar("i64", dimSize, H5T_NATIVE_LLONG, global_dims, offset, count,
                                       currentTestData.I64.data());
#else
            h5writer.CreateAndStoreVar("i64", dimSize, H5T_NATIVE_LONG, global_dims, offset, count,
                                       currentTestData.I64.data());
#endif
            h5writer.CreateAndStoreVar("u8", dimSize, H5T_NATIVE_UCHAR, global_dims, offset, count,
                                       currentTestData.U8.data());
            h5writer.CreateAndStoreVar("u16", dimSize, H5T_NATIVE_USHORT, global_dims, offset,
                                       count, currentTestData.U16.data());
            h5writer.CreateAndStoreVar("u32", dimSize, H5T_NATIVE_UINT, global_dims, offset, count,
                                       currentTestData.U32.data());
#ifdef _WIN32
            h5writer.CreateAndStoreVar("u64", dimSize, H5T_NATIVE_ULLONG, global_dims, offset,
                                       count, currentTestData.U64.data());
#else
            h5writer.CreateAndStoreVar("u64", dimSize, H5T_NATIVE_ULONG, global_dims, offset, count,
                                       currentTestData.U64.data());
#endif
            h5writer.CreateAndStoreVar("r32", dimSize, H5T_NATIVE_FLOAT, global_dims, offset, count,
                                       currentTestData.R32.data());
            h5writer.CreateAndStoreVar("r64", dimSize, H5T_NATIVE_DOUBLE, global_dims, offset,
                                       count, currentTestData.R64.data());
            h5writer.Advance();
        }
    }

    { // ADIOS2 read back

        // Write test data using ADIOS2

#ifdef TEST_HDF5_MPI
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif
        adios2::IO io = adios.DeclareIO("HDF5ReadIO");
        io.SetEngine("HDF5");

        adios2::Engine hdf5Reader = io.Open(fname, adios2::Mode::Read);

        auto var_iString = io.InquireVariable<std::string>("iString");
        EXPECT_TRUE(var_iString);
        ASSERT_EQ(var_iString.Shape().size(), 0);
        ASSERT_EQ(var_iString.Steps(), NSteps);

        auto var_ch = io.InquireVariable<uint8_t>("ch");
        EXPECT_TRUE(var_ch);
        ASSERT_EQ(var_ch.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_ch.Steps(), NSteps);
        ASSERT_EQ(var_ch.Shape()[0], mpiSize * Nx);

        auto var_i8 = io.InquireVariable<int8_t>("i8");
        EXPECT_TRUE(var_i8);
        ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i8.Steps(), NSteps);
        ASSERT_EQ(var_i8.Shape()[0], mpiSize * Nx);

        auto var_i16 = io.InquireVariable<int16_t>("i16");
        EXPECT_TRUE(var_i16);
        ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i16.Steps(), NSteps);
        ASSERT_EQ(var_i16.Shape()[0], mpiSize * Nx);

        auto var_i32 = io.InquireVariable<int32_t>("i32");
        EXPECT_TRUE(var_i32);
        ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i32.Steps(), NSteps);
        ASSERT_EQ(var_i32.Shape()[0], mpiSize * Nx);

        auto var_i64 = io.InquireVariable<int64_t>("i64");
        EXPECT_TRUE(var_i64);
        ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i64.Steps(), NSteps);
        ASSERT_EQ(var_i64.Shape()[0], mpiSize * Nx);

        auto var_u8 = io.InquireVariable<uint8_t>("u8");
        EXPECT_TRUE(var_u8);
        ASSERT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u8.Steps(), NSteps);
        ASSERT_EQ(var_u8.Shape()[0], mpiSize * Nx);

        auto var_u16 = io.InquireVariable<uint16_t>("u16");
        EXPECT_TRUE(var_u16);
        ASSERT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u16.Steps(), NSteps);
        ASSERT_EQ(var_u16.Shape()[0], mpiSize * Nx);

        auto var_u32 = io.InquireVariable<uint32_t>("u32");
        EXPECT_TRUE(var_u32);
        ASSERT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u32.Steps(), NSteps);
        ASSERT_EQ(var_u32.Shape()[0], mpiSize * Nx);

        auto var_u64 = io.InquireVariable<uint64_t>("u64");
        EXPECT_TRUE(var_u64);
        ASSERT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u64.Steps(), NSteps);
        ASSERT_EQ(var_u64.Shape()[0], mpiSize * Nx);

        auto var_r32 = io.InquireVariable<float>("r32");
        EXPECT_TRUE(var_r32);
        ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r32.Steps(), NSteps);
        ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);

        auto var_r64 = io.InquireVariable<double>("r64");
        EXPECT_TRUE(var_r64);
        ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r64.Steps(), NSteps);
        ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);

        // TODO: other types

        SmallTestData testData;

        std::string IString;
        std::array<uint8_t, Nx> CHAR;
        std::array<int8_t, Nx> I8;
        std::array<int16_t, Nx> I16;
        std::array<int32_t, Nx> I32;
        std::array<int64_t, Nx> I64;
        std::array<uint8_t, Nx> U8;
        std::array<uint16_t, Nx> U16;
        std::array<uint32_t, Nx> U32;
        std::array<uint64_t, Nx> U64;
        std::array<float, Nx> R32;
        std::array<double, Nx> R64;

        const adios2::Dims start{mpiRank * Nx};
        const adios2::Dims count{Nx};

        const adios2::Box<adios2::Dims> sel(start, count);

        var_ch.SetSelection(sel);
        var_i8.SetSelection(sel);
        var_i16.SetSelection(sel);
        var_i32.SetSelection(sel);
        var_i64.SetSelection(sel);

        var_u8.SetSelection(sel);
        var_u16.SetSelection(sel);
        var_u32.SetSelection(sel);
        var_u64.SetSelection(sel);

        var_r32.SetSelection(sel);
        var_r64.SetSelection(sel);

        for (size_t t = 0; t < NSteps; ++t)
        {
            var_ch.SetStepSelection({t, 1});
            var_i8.SetStepSelection({t, 1});
            var_i16.SetStepSelection({t, 1});
            var_i32.SetStepSelection({t, 1});
            var_i64.SetStepSelection({t, 1});

            var_u8.SetStepSelection({t, 1});
            var_u16.SetStepSelection({t, 1});
            var_u32.SetStepSelection({t, 1});
            var_u64.SetStepSelection({t, 1});

            var_r32.SetStepSelection({t, 1});
            var_r64.SetStepSelection({t, 1});

            // Generate test data for each rank uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(t), mpiRank, mpiSize);

            hdf5Reader.Get(var_iString, IString);

            hdf5Reader.Get(var_ch, CHAR.data());
            hdf5Reader.Get(var_i8, I8.data());
            hdf5Reader.Get(var_i16, I16.data());
            hdf5Reader.Get(var_i32, I32.data());
            hdf5Reader.Get(var_i64, I64.data());

            hdf5Reader.Get(var_u8, U8.data());
            hdf5Reader.Get(var_u16, U16.data());
            hdf5Reader.Get(var_u32, U32.data());
            hdf5Reader.Get(var_u64, U64.data());

            hdf5Reader.Get(var_r32, R32.data());
            hdf5Reader.Get(var_r64, R64.data());

            hdf5Reader.PerformGets();

            EXPECT_EQ(IString, currentTestData.S1);

            for (size_t i = 0; i < Nx; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                EXPECT_EQ(static_cast<char>(CHAR[i]), currentTestData.CHAR[i]) << msg;
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
        }
        hdf5Reader.Close();
    }
}

//******************************************************************************
// 2D 2x4 test data
//******************************************************************************

// ADIOS2 write, native HDF5 read
TEST_F(HDF5WriteReadTest, ADIOS2HDF5WriteHDF5Read2D2x4)
{
    // Each process would write a 2x4 array and all processes would
    // form a 2D 2 * (numberOfProcess*Nx) matrix where Nx is 4 here
    std::string fname = "ADIOS2HDF5WriteHDF5Read2D2x4Test.h5";

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const std::size_t Nx = 4;

    // Number of rows
    const std::size_t Ny = 2;

    // Number of steps
    const std::size_t NSteps = 3;

#ifdef TEST_HDF5_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // Write test data using ADIOS2
    {
#ifdef TEST_HDF5_MPI
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif
        adios2::IO io = adios.DeclareIO("TestIO");

        // Declare 2D variables (Ny * (NumOfProcesses * Nx))
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            adios2::Dims shape{static_cast<unsigned int>(Ny),
                               static_cast<unsigned int>(Nx * mpiSize)};
            adios2::Dims start{static_cast<unsigned int>(0),
                               static_cast<unsigned int>(mpiRank * Nx)};
            adios2::Dims count{static_cast<unsigned int>(Ny), static_cast<unsigned int>(Nx)};
            auto var_iString = io.DefineVariable<std::string>("iString");
            EXPECT_TRUE(var_iString);
            auto var_i8 = io.DefineVariable<int8_t>("i8", shape, start, count);
            EXPECT_TRUE(var_i8);
            auto var_i16 = io.DefineVariable<int16_t>("i16", shape, start, count);
            EXPECT_TRUE(var_i16);
            auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count);
            EXPECT_TRUE(var_i32);
            auto var_i64 = io.DefineVariable<int64_t>("i64", shape, start, count);
            EXPECT_TRUE(var_i64);
            auto var_u8 = io.DefineVariable<uint8_t>("u8", shape, start, count);
            EXPECT_TRUE(var_u8);
            auto var_u16 = io.DefineVariable<uint16_t>("u16", shape, start, count);
            EXPECT_TRUE(var_u16);
            auto var_u32 = io.DefineVariable<uint32_t>("u32", shape, start, count);
            EXPECT_TRUE(var_u32);
            auto var_u64 = io.DefineVariable<uint64_t>("u64", shape, start, count);
            EXPECT_TRUE(var_u64);
            auto var_r32 = io.DefineVariable<float>("r32", shape, start, count);
            EXPECT_TRUE(var_r32);
            auto var_r64 = io.DefineVariable<double>("r64", shape, start, count);
            EXPECT_TRUE(var_r64);
        }

        // Create the HDF5 Engine
        io.SetEngine("HDF5");

        // HDf5 engine calls the HDF5 common object that calls the hDF5 library.
        // The IO functionality, SetParameters and AddTransports will be added
        // in the future. For now `io.AddTransport("file", {
        // "library", "MPI"}});` is omitted.
        // })
        io.AddTransport("file");

        adios2::Engine engine = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, step, mpiRank, mpiSize);

            // Retrieve the variables that previously went out of scope
            auto var_iString = io.InquireVariable<std::string>("iString");
            auto var_i8 = io.InquireVariable<int8_t>("i8");
            auto var_i16 = io.InquireVariable<int16_t>("i16");
            auto var_i32 = io.InquireVariable<int32_t>("i32");
            auto var_i64 = io.InquireVariable<int64_t>("i64");
            auto var_u8 = io.InquireVariable<uint8_t>("u8");
            auto var_u16 = io.InquireVariable<uint16_t>("u16");
            auto var_u32 = io.InquireVariable<uint32_t>("u32");
            auto var_u64 = io.InquireVariable<uint64_t>("u64");
            auto var_r32 = io.InquireVariable<float>("r32");
            auto var_r64 = io.InquireVariable<double>("r64");

            // Make a 2D selection to describe the local dimensions of the
            // variable we write and its offsets in the global spaces
            adios2::Box<adios2::Dims> sel({0, static_cast<unsigned int>(mpiRank * Nx)}, {Ny, Nx});
            var_i8.SetSelection(sel);
            var_i16.SetSelection(sel);
            var_i32.SetSelection(sel);
            var_i64.SetSelection(sel);
            var_u8.SetSelection(sel);
            var_u16.SetSelection(sel);
            var_u32.SetSelection(sel);
            var_u64.SetSelection(sel);
            var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            // Write each one
            // fill in the variable with values from starting index to
            // starting index + count
            engine.Put(var_iString, currentTestData.S1);
            engine.Put(var_i8, currentTestData.I8.data());
            engine.Put(var_i16, currentTestData.I16.data());
            engine.Put(var_i32, currentTestData.I32.data());
            engine.Put(var_i64, currentTestData.I64.data());
            engine.Put(var_u8, currentTestData.U8.data());
            engine.Put(var_u16, currentTestData.U16.data());
            engine.Put(var_u32, currentTestData.U32.data());
            engine.Put(var_u64, currentTestData.U64.data());
            engine.Put(var_r32, currentTestData.R32.data());
            engine.Put(var_r64, currentTestData.R64.data());

            // Advance to the next time step
            engine.EndStep();
        }

        // Close the file
        engine.Close();
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
#ifdef _WIN32
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_LLONG), 1);
#else
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_LONG), 1);
#endif
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
#ifdef _WIN32
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_ULLONG), 1);
#else
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_ULONG), 1);
#endif
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

// ADIOS2 write, ADIOS2 read
TEST_F(HDF5WriteReadTest, ADIOS2HDF5WriteADIOS2HDF5Read2D2x4)
{
    std::string fname = "ADIOS2HDF5WriteADIOS2HDF5Read2D2x4Test.h5";
    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const std::size_t Nx = 4;

    // Number of rows
    const std::size_t Ny = 2;

    // Number of steps
    const std::size_t NSteps = 3;

#ifdef TEST_HDF5_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

#ifdef TEST_HDF5_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    // Write test data using ADIOS2
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        // Declare 2D variables (Ny * (NumOfProcesses * Nx))
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            adios2::Dims shape{static_cast<unsigned int>(Ny),
                               static_cast<unsigned int>(Nx * mpiSize)};
            adios2::Dims start{static_cast<unsigned int>(0),
                               static_cast<unsigned int>(mpiRank * Nx)};
            adios2::Dims count{static_cast<unsigned int>(Ny), static_cast<unsigned int>(Nx)};

            auto var_iString = io.DefineVariable<std::string>("iString");
            EXPECT_TRUE(var_iString);
            auto var_i8 = io.DefineVariable<int8_t>("i8", shape, start, count);
            EXPECT_TRUE(var_i8);
            auto var_i16 = io.DefineVariable<int16_t>("i16", shape, start, count);
            EXPECT_TRUE(var_i16);
            auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count);
            EXPECT_TRUE(var_i32);
            auto var_i64 = io.DefineVariable<int64_t>("i64", shape, start, count);
            EXPECT_TRUE(var_i64);
            auto var_u8 = io.DefineVariable<uint8_t>("u8", shape, start, count);
            EXPECT_TRUE(var_u8);
            auto var_u16 = io.DefineVariable<uint16_t>("u16", shape, start, count);
            EXPECT_TRUE(var_u16);
            auto var_u32 = io.DefineVariable<uint32_t>("u32", shape, start, count);
            EXPECT_TRUE(var_u32);
            auto var_u64 = io.DefineVariable<uint64_t>("u64", shape, start, count);
            EXPECT_TRUE(var_u64);
            auto var_r32 = io.DefineVariable<float>("r32", shape, start, count);
            EXPECT_TRUE(var_r32);
            auto var_r64 = io.DefineVariable<double>("r64", shape, start, count);
            EXPECT_TRUE(var_r64);
        }

        // Create the HDF5 Engine
        io.SetEngine("HDF5");

        // HDf5 engine calls the HDF5 common object that calls the hDF5 library.
        // The IO functionality, SetParameters and AddTransports will be added
        // in the future. For now `io.AddTransport("file", {
        // "library", "MPI"}});` is omitted.
        // })
        io.AddTransport("file");

        adios2::Engine engine = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, step, mpiRank, mpiSize);

            // Retrieve the variables that previously went out of scope
            auto var_iString = io.InquireVariable<std::string>("iString");
            auto var_i8 = io.InquireVariable<int8_t>("i8");
            auto var_i16 = io.InquireVariable<int16_t>("i16");
            auto var_i32 = io.InquireVariable<int32_t>("i32");
            auto var_i64 = io.InquireVariable<int64_t>("i64");
            auto var_u8 = io.InquireVariable<uint8_t>("u8");
            auto var_u16 = io.InquireVariable<uint16_t>("u16");
            auto var_u32 = io.InquireVariable<uint32_t>("u32");
            auto var_u64 = io.InquireVariable<uint64_t>("u64");
            auto var_r32 = io.InquireVariable<float>("r32");
            auto var_r64 = io.InquireVariable<double>("r64");

            // Make a 2D selection to describe the local dimensions of the
            // variable we write and its offsets in the global spaces
            adios2::Box<adios2::Dims> sel({0, static_cast<unsigned int>(mpiRank * Nx)}, {Ny, Nx});
            var_i8.SetSelection(sel);
            var_i16.SetSelection(sel);
            var_i32.SetSelection(sel);
            var_i64.SetSelection(sel);
            var_u8.SetSelection(sel);
            var_u16.SetSelection(sel);
            var_u32.SetSelection(sel);
            var_u64.SetSelection(sel);
            var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            // Write each one
            // fill in the variable with values from starting index to
            // starting index + count
            engine.Put(var_iString, currentTestData.S1);
            engine.Put(var_i8, currentTestData.I8.data());
            engine.Put(var_i16, currentTestData.I16.data());
            engine.Put(var_i32, currentTestData.I32.data());
            engine.Put(var_i64, currentTestData.I64.data());
            engine.Put(var_u8, currentTestData.U8.data());
            engine.Put(var_u16, currentTestData.U16.data());
            engine.Put(var_u32, currentTestData.U32.data());
            engine.Put(var_u64, currentTestData.U64.data());
            engine.Put(var_r32, currentTestData.R32.data());
            engine.Put(var_r64, currentTestData.R64.data());

            // Advance to the next time step
            engine.EndStep();
        }

        // Close the file
        engine.Close();
    }

    { // reading back
        adios2::IO io = adios.DeclareIO("HDF5ReadIO");
        io.SetEngine("HDF5");

        adios2::Engine hdf5Reader = io.Open(fname, adios2::Mode::Read);

        auto var_iString = io.InquireVariable<std::string>("iString");
        EXPECT_TRUE(var_iString);
        ASSERT_EQ(var_iString.Shape().size(), 0);
        ASSERT_EQ(var_iString.Steps(), NSteps);

        auto var_i8 = io.InquireVariable<int8_t>("i8");
        EXPECT_TRUE(var_i8);
        ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i8.Steps(), NSteps);
        ASSERT_EQ(var_i8.Shape()[0], Ny);
        ASSERT_EQ(var_i8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i16 = io.InquireVariable<int16_t>("i16");
        EXPECT_TRUE(var_i16);
        ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i16.Steps(), NSteps);
        ASSERT_EQ(var_i16.Shape()[0], Ny);
        ASSERT_EQ(var_i16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i32 = io.InquireVariable<int32_t>("i32");
        EXPECT_TRUE(var_i32);
        ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i32.Steps(), NSteps);
        ASSERT_EQ(var_i32.Shape()[0], Ny);
        ASSERT_EQ(var_i32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i64 = io.InquireVariable<int64_t>("i64");
        EXPECT_TRUE(var_i64);
        ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i64.Steps(), NSteps);
        ASSERT_EQ(var_i64.Shape()[0], Ny);
        ASSERT_EQ(var_i64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u8 = io.InquireVariable<uint8_t>("u8");
        EXPECT_TRUE(var_u8);
        ASSERT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u8.Steps(), NSteps);
        ASSERT_EQ(var_u8.Shape()[0], Ny);
        ASSERT_EQ(var_u8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u16 = io.InquireVariable<uint16_t>("u16");
        EXPECT_TRUE(var_u16);
        ASSERT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u16.Steps(), NSteps);
        ASSERT_EQ(var_u16.Shape()[0], Ny);
        ASSERT_EQ(var_u16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u32 = io.InquireVariable<uint32_t>("u32");
        EXPECT_TRUE(var_u32);
        ASSERT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u32.Steps(), NSteps);
        ASSERT_EQ(var_u32.Shape()[0], Ny);
        ASSERT_EQ(var_u32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u64 = io.InquireVariable<uint64_t>("u64");
        EXPECT_TRUE(var_u64);
        ASSERT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u64.Steps(), NSteps);
        ASSERT_EQ(var_u64.Shape()[0], Ny);
        ASSERT_EQ(var_u64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_r32 = io.InquireVariable<float>("r32");
        EXPECT_TRUE(var_r32);
        ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r32.Steps(), NSteps);
        ASSERT_EQ(var_r32.Shape()[0], Ny);
        ASSERT_EQ(var_r32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_r64 = io.InquireVariable<double>("r64");
        EXPECT_TRUE(var_r64);
        ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r64.Steps(), NSteps);
        ASSERT_EQ(var_r64.Shape()[0], Ny);
        ASSERT_EQ(var_r64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        std::string IString;
        std::array<int8_t, Nx * Ny> I8;
        std::array<int16_t, Nx * Ny> I16;
        std::array<int32_t, Nx * Ny> I32;
        std::array<int64_t, Nx * Ny> I64;
        std::array<uint8_t, Nx * Ny> U8;
        std::array<uint16_t, Nx * Ny> U16;
        std::array<uint32_t, Nx * Ny> U32;
        std::array<uint64_t, Nx * Ny> U64;
        std::array<float, Nx * Ny> R32;
        std::array<double, Nx * Ny> R64;

        const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
        const adios2::Dims count{Ny, Nx};

        const adios2::Box<adios2::Dims> sel(start, count);

        var_i8.SetSelection(sel);
        var_i16.SetSelection(sel);
        var_i32.SetSelection(sel);
        var_i64.SetSelection(sel);

        var_u8.SetSelection(sel);
        var_u16.SetSelection(sel);
        var_u32.SetSelection(sel);
        var_u64.SetSelection(sel);

        var_r32.SetSelection(sel);
        var_r64.SetSelection(sel);

        for (size_t t = 0; t < NSteps; ++t)
        {
            var_i8.SetStepSelection({t, 1});
            var_i16.SetStepSelection({t, 1});
            var_i32.SetStepSelection({t, 1});
            var_i64.SetStepSelection({t, 1});

            var_u8.SetStepSelection({t, 1});
            var_u16.SetStepSelection({t, 1});
            var_u32.SetStepSelection({t, 1});
            var_u64.SetStepSelection({t, 1});

            var_r32.SetStepSelection({t, 1});
            var_r64.SetStepSelection({t, 1});

            hdf5Reader.Get(var_iString, IString);

            hdf5Reader.Get(var_i8, I8.data());
            hdf5Reader.Get(var_i16, I16.data());
            hdf5Reader.Get(var_i32, I32.data());
            hdf5Reader.Get(var_i64, I64.data());

            hdf5Reader.Get(var_u8, U8.data());
            hdf5Reader.Get(var_u16, U16.data());
            hdf5Reader.Get(var_u32, U32.data());
            hdf5Reader.Get(var_u64, U64.data());

            hdf5Reader.Get(var_r32, R32.data());
            hdf5Reader.Get(var_r64, R64.data());

            hdf5Reader.PerformGets();

            // Generate test data for each rank uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(t), mpiRank, mpiSize);

            EXPECT_EQ(IString, currentTestData.S1);

            for (size_t i = 0; i < Nx * Ny; ++i)
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
        }
        hdf5Reader.Close();
    }
}

// Native HDF5 write, ADIOS2 read
TEST_F(HDF5WriteReadTest, HDF5WriteADIOS2HDF5Read2D2x4)
{
    std::string fname = "HDF5WriteADIOS2HDF5Read2D2x4Test.h5";

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const std::size_t Nx = 4;
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
#ifdef _WIN32
            h5writer.CreateAndStoreVar("i64", dimSize, H5T_NATIVE_LLONG, global_dims, offset, count,
                                       currentTestData.I64.data());
#else
            h5writer.CreateAndStoreVar("i64", dimSize, H5T_NATIVE_LONG, global_dims, offset, count,
                                       currentTestData.I64.data());
#endif
            h5writer.CreateAndStoreVar("u8", dimSize, H5T_NATIVE_UCHAR, global_dims, offset, count,
                                       currentTestData.U8.data());
            h5writer.CreateAndStoreVar("u16", dimSize, H5T_NATIVE_USHORT, global_dims, offset,
                                       count, currentTestData.U16.data());
            h5writer.CreateAndStoreVar("u32", dimSize, H5T_NATIVE_UINT, global_dims, offset, count,
                                       currentTestData.U32.data());
#ifdef _WIN32
            h5writer.CreateAndStoreVar("u64", dimSize, H5T_NATIVE_ULLONG, global_dims, offset,
                                       count, currentTestData.U64.data());
#else
            h5writer.CreateAndStoreVar("u64", dimSize, H5T_NATIVE_ULONG, global_dims, offset, count,
                                       currentTestData.U64.data());
#endif
            h5writer.CreateAndStoreVar("r32", dimSize, H5T_NATIVE_FLOAT, global_dims, offset, count,
                                       currentTestData.R32.data());
            h5writer.CreateAndStoreVar("r64", dimSize, H5T_NATIVE_DOUBLE, global_dims, offset,
                                       count, currentTestData.R64.data());
            h5writer.Advance();
        }
    }

    { // read back
#ifdef TEST_HDF5_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

#ifdef TEST_HDF5_MPI
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif
        adios2::IO io = adios.DeclareIO("HDF5ReadIO");
        io.SetEngine("HDF5");

        adios2::Engine hdf5Reader = io.Open(fname, adios2::Mode::Read);

        auto var_iString = io.InquireVariable<std::string>("iString");
        EXPECT_TRUE(var_iString);
        ASSERT_EQ(var_iString.Shape().size(), 0);
        ASSERT_EQ(var_iString.Steps(), NSteps);

        auto var_i8 = io.InquireVariable<int8_t>("i8");
        EXPECT_TRUE(var_i8);
        ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i8.Steps(), NSteps);
        ASSERT_EQ(var_i8.Shape()[0], Ny);
        ASSERT_EQ(var_i8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i16 = io.InquireVariable<int16_t>("i16");
        EXPECT_TRUE(var_i16);
        ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i16.Steps(), NSteps);
        ASSERT_EQ(var_i16.Shape()[0], Ny);
        ASSERT_EQ(var_i16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i32 = io.InquireVariable<int32_t>("i32");
        EXPECT_TRUE(var_i32);
        ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i32.Steps(), NSteps);
        ASSERT_EQ(var_i32.Shape()[0], Ny);
        ASSERT_EQ(var_i32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i64 = io.InquireVariable<int64_t>("i64");
        EXPECT_TRUE(var_i64);
        ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i64.Steps(), NSteps);
        ASSERT_EQ(var_i64.Shape()[0], Ny);
        ASSERT_EQ(var_i64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u8 = io.InquireVariable<uint8_t>("u8");
        EXPECT_TRUE(var_u8);
        ASSERT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u8.Steps(), NSteps);
        ASSERT_EQ(var_u8.Shape()[0], Ny);
        ASSERT_EQ(var_u8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u16 = io.InquireVariable<uint16_t>("u16");
        EXPECT_TRUE(var_u16);
        ASSERT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u16.Steps(), NSteps);
        ASSERT_EQ(var_u16.Shape()[0], Ny);
        ASSERT_EQ(var_u16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u32 = io.InquireVariable<uint32_t>("u32");
        EXPECT_TRUE(var_u32);
        ASSERT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u32.Steps(), NSteps);
        ASSERT_EQ(var_u32.Shape()[0], Ny);
        ASSERT_EQ(var_u32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u64 = io.InquireVariable<uint64_t>("u64");
        EXPECT_TRUE(var_u64);
        ASSERT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u64.Steps(), NSteps);
        ASSERT_EQ(var_u64.Shape()[0], Ny);
        ASSERT_EQ(var_u64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_r32 = io.InquireVariable<float>("r32");
        EXPECT_TRUE(var_r32);
        ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r32.Steps(), NSteps);
        ASSERT_EQ(var_r32.Shape()[0], Ny);
        ASSERT_EQ(var_r32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_r64 = io.InquireVariable<double>("r64");
        EXPECT_TRUE(var_r64);
        ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r64.Steps(), NSteps);
        ASSERT_EQ(var_r64.Shape()[0], Ny);
        ASSERT_EQ(var_r64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        std::string IString;
        std::array<int8_t, Nx * Ny> I8;
        std::array<int16_t, Nx * Ny> I16;
        std::array<int32_t, Nx * Ny> I32;
        std::array<int64_t, Nx * Ny> I64;
        std::array<uint8_t, Nx * Ny> U8;
        std::array<uint16_t, Nx * Ny> U16;
        std::array<uint32_t, Nx * Ny> U32;
        std::array<uint64_t, Nx * Ny> U64;
        std::array<float, Nx * Ny> R32;
        std::array<double, Nx * Ny> R64;

        const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
        const adios2::Dims count{Ny, Nx};

        const adios2::Box<adios2::Dims> sel(start, count);

        var_i8.SetSelection(sel);
        var_i16.SetSelection(sel);
        var_i32.SetSelection(sel);
        var_i64.SetSelection(sel);

        var_u8.SetSelection(sel);
        var_u16.SetSelection(sel);
        var_u32.SetSelection(sel);
        var_u64.SetSelection(sel);

        var_r32.SetSelection(sel);
        var_r64.SetSelection(sel);

        for (size_t t = 0; t < NSteps; ++t)
        {
            var_i8.SetStepSelection({t, 1});
            var_i16.SetStepSelection({t, 1});
            var_i32.SetStepSelection({t, 1});
            var_i64.SetStepSelection({t, 1});

            var_u8.SetStepSelection({t, 1});
            var_u16.SetStepSelection({t, 1});
            var_u32.SetStepSelection({t, 1});
            var_u64.SetStepSelection({t, 1});

            var_r32.SetStepSelection({t, 1});
            var_r64.SetStepSelection({t, 1});

            hdf5Reader.Get(var_iString, IString);

            hdf5Reader.Get(var_i8, I8.data());
            hdf5Reader.Get(var_i16, I16.data());
            hdf5Reader.Get(var_i32, I32.data());
            hdf5Reader.Get(var_i64, I64.data());

            hdf5Reader.Get(var_u8, U8.data());
            hdf5Reader.Get(var_u16, U16.data());
            hdf5Reader.Get(var_u32, U32.data());
            hdf5Reader.Get(var_u64, U64.data());

            hdf5Reader.Get(var_r32, R32.data());
            hdf5Reader.Get(var_r64, R64.data());

            hdf5Reader.PerformGets();

            // Generate test data for each rank uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(t), mpiRank, mpiSize);

            EXPECT_EQ(IString, currentTestData.S1);

            for (size_t i = 0; i < Nx * Ny; ++i)
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
        }
        hdf5Reader.Close();
    }
}

//******************************************************************************
// 2D 4x2 test data
//******************************************************************************

// ADIOS2 write, native HDF5 read
TEST_F(HDF5WriteReadTest, ADIOS2HDF5WriteHDF5Read2D4x2)
{

    // Each process would write a 4x2 array and all processes would
    // form a 2D 4 * (NumberOfProcess * Nx) matrix where Nx is 2 here
    std::string fname = "ADIOS2HDF5WriteHDF5Read2D4x2Test.h5";

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const std::size_t Nx = 2;
    // Number of cols
    const std::size_t Ny = 4;

    // Number of steps
    const std::size_t NSteps = 3;

#ifdef TEST_HDF5_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // Write test data using ADIOS2
    {
#ifdef TEST_HDF5_MPI
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif
        adios2::IO io = adios.DeclareIO("TestIO");

        // Declare 2D variables (4 * (NumberOfProcess * Nx))
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            adios2::Dims shape{static_cast<unsigned int>(Ny),
                               static_cast<unsigned int>(mpiSize * Nx)};
            adios2::Dims start{static_cast<unsigned int>(0),
                               static_cast<unsigned int>(mpiRank * Nx)};
            adios2::Dims count{static_cast<unsigned int>(Ny), static_cast<unsigned int>(Nx)};

            auto var_iString = io.DefineVariable<std::string>("iString");
            EXPECT_TRUE(var_iString);
            auto var_i8 = io.DefineVariable<int8_t>("i8", shape, start, count);
            EXPECT_TRUE(var_i8);
            auto var_i16 = io.DefineVariable<int16_t>("i16", shape, start, count);
            EXPECT_TRUE(var_i16);
            auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count);
            EXPECT_TRUE(var_i32);
            auto var_i64 = io.DefineVariable<int64_t>("i64", shape, start, count);
            EXPECT_TRUE(var_i64);
            auto var_u8 = io.DefineVariable<uint8_t>("u8", shape, start, count);
            EXPECT_TRUE(var_u8);
            auto var_u16 = io.DefineVariable<uint16_t>("u16", shape, start, count);
            EXPECT_TRUE(var_u16);
            auto var_u32 = io.DefineVariable<uint32_t>("u32", shape, start, count);
            EXPECT_TRUE(var_u32);
            auto var_u64 = io.DefineVariable<uint64_t>("u64", shape, start, count);
            EXPECT_TRUE(var_u64);
            auto var_r32 = io.DefineVariable<float>("r32", shape, start, count);
            EXPECT_TRUE(var_r32);
            auto var_r64 = io.DefineVariable<double>("r64", shape, start, count);
            EXPECT_TRUE(var_r64);
        }

        // Create the HDF5 Engine
        io.SetEngine("HDF5");

        // HDf5 engine calls the HDF5 common object that calls the hDF5 library.
        // The IO functionality, SetParameters and AddTransports will be added
        // in the future. For now `io.AddTransport("file", {
        // "library", "MPI"}});` is omitted.
        // })
        io.AddTransport("file");

        adios2::Engine engine = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, step, mpiRank, mpiSize);

            // Retrieve the variables that previously went out of scope
            auto var_iString = io.InquireVariable<std::string>("iString");
            auto var_i8 = io.InquireVariable<int8_t>("i8");
            auto var_i16 = io.InquireVariable<int16_t>("i16");
            auto var_i32 = io.InquireVariable<int32_t>("i32");
            auto var_i64 = io.InquireVariable<int64_t>("i64");
            auto var_u8 = io.InquireVariable<uint8_t>("u8");
            auto var_u16 = io.InquireVariable<uint16_t>("u16");
            auto var_u32 = io.InquireVariable<uint32_t>("u32");
            auto var_u64 = io.InquireVariable<uint64_t>("u64");
            auto var_r32 = io.InquireVariable<float>("r32");
            auto var_r64 = io.InquireVariable<double>("r64");

            // Make a 2D selection to describe the local dimensions of the
            // variable we write and its offsets in the global spaces
            adios2::Box<adios2::Dims> sel({0, static_cast<unsigned int>(mpiRank * Nx)}, {Ny, Nx});
            var_i8.SetSelection(sel);
            var_i16.SetSelection(sel);
            var_i32.SetSelection(sel);
            var_i64.SetSelection(sel);
            var_u8.SetSelection(sel);
            var_u16.SetSelection(sel);
            var_u32.SetSelection(sel);
            var_u64.SetSelection(sel);
            var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            // Write each one
            // fill in the variable with values from starting index to
            // starting index + count
            engine.BeginStep();
            engine.Put(var_iString, currentTestData.S1);
            engine.Put(var_i8, currentTestData.I8.data());
            engine.Put(var_i16, currentTestData.I16.data());
            engine.Put(var_i32, currentTestData.I32.data());
            engine.Put(var_i64, currentTestData.I64.data());
            engine.Put(var_u8, currentTestData.U8.data());
            engine.Put(var_u16, currentTestData.U16.data());
            engine.Put(var_u32, currentTestData.U32.data());
            engine.Put(var_u64, currentTestData.U64.data());
            engine.Put(var_r32, currentTestData.R32.data());
            engine.Put(var_r64, currentTestData.R64.data());
            engine.EndStep();
        }

        engine.Close();
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
#ifdef _WIN32
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_LLONG), 1);
#else
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_LONG), 1);
#endif
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
#ifdef _WIN32
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_ULLONG), 1);
#else
            ASSERT_EQ(H5Tequal(h5Type, H5T_NATIVE_ULONG), 1);
#endif
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

// ADIOS2 write, ADIOS2 read
TEST_F(HDF5WriteReadTest, ADIOS2HDF5WriteADIOS2HDF5Read2D4x2)
{
    std::string fname = "ADIOS2HDF5WriteADIOS2HDF5Read2D4x2Test.h5";

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const std::size_t Nx = 2;
    // Number of cols
    const std::size_t Ny = 4;

    // Number of steps
    const std::size_t NSteps = 3;

#ifdef TEST_HDF5_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

#ifdef TEST_HDF5_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif

    // Write test data using ADIOS2
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        // Declare 2D variables (4 * (NumberOfProcess * Nx))
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            adios2::Dims shape{static_cast<unsigned int>(Ny),
                               static_cast<unsigned int>(mpiSize * Nx)};
            adios2::Dims start{static_cast<unsigned int>(0),
                               static_cast<unsigned int>(mpiRank * Nx)};
            adios2::Dims count{static_cast<unsigned int>(Ny), static_cast<unsigned int>(Nx)};
            auto var_iString = io.DefineVariable<std::string>("iString");
            EXPECT_TRUE(var_iString);
            auto var_i8 = io.DefineVariable<int8_t>("i8", shape, start, count);
            EXPECT_TRUE(var_i8);
            auto var_i16 = io.DefineVariable<int16_t>("i16", shape, start, count);
            EXPECT_TRUE(var_i16);
            auto var_i32 = io.DefineVariable<int32_t>("i32", shape, start, count);
            EXPECT_TRUE(var_i32);
            auto var_i64 = io.DefineVariable<int64_t>("i64", shape, start, count);
            EXPECT_TRUE(var_i64);
            auto var_u8 = io.DefineVariable<uint8_t>("u8", shape, start, count);
            EXPECT_TRUE(var_u8);
            auto var_u16 = io.DefineVariable<uint16_t>("u16", shape, start, count);
            EXPECT_TRUE(var_u16);
            auto var_u32 = io.DefineVariable<uint32_t>("u32", shape, start, count);
            EXPECT_TRUE(var_u32);
            auto var_u64 = io.DefineVariable<uint64_t>("u64", shape, start, count);
            EXPECT_TRUE(var_u64);
            auto var_r32 = io.DefineVariable<float>("r32", shape, start, count);
            EXPECT_TRUE(var_r32);
            auto var_r64 = io.DefineVariable<double>("r64", shape, start, count);
            EXPECT_TRUE(var_r64);
        }

        // Create the HDF5 Engine
        io.SetEngine("HDF5");

        // HDf5 engine calls the HDF5 common object that calls the hDF5 library.
        // The IO functionality, SetParameters and AddTransports will be added
        // in the future. For now `io.AddTransport("file", {
        // "library", "MPI"}});` is omitted.
        // })
        io.AddTransport("file");

        adios2::Engine engine = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, step, mpiRank, mpiSize);

            // Retrieve the variables that previously went out of scope
            auto var_iString = io.InquireVariable<std::string>("iString");
            auto var_i8 = io.InquireVariable<int8_t>("i8");
            auto var_i16 = io.InquireVariable<int16_t>("i16");
            auto var_i32 = io.InquireVariable<int32_t>("i32");
            auto var_i64 = io.InquireVariable<int64_t>("i64");
            auto var_u8 = io.InquireVariable<uint8_t>("u8");
            auto var_u16 = io.InquireVariable<uint16_t>("u16");
            auto var_u32 = io.InquireVariable<uint32_t>("u32");
            auto var_u64 = io.InquireVariable<uint64_t>("u64");
            auto var_r32 = io.InquireVariable<float>("r32");
            auto var_r64 = io.InquireVariable<double>("r64");

            // Make a 2D selection to describe the local dimensions of the
            // variable we write and its offsets in the global spaces
            adios2::Box<adios2::Dims> sel({0, static_cast<unsigned int>(mpiRank * Nx)}, {Ny, Nx});
            var_i8.SetSelection(sel);
            var_i16.SetSelection(sel);
            var_i32.SetSelection(sel);
            var_i64.SetSelection(sel);
            var_u8.SetSelection(sel);
            var_u16.SetSelection(sel);
            var_u32.SetSelection(sel);
            var_u64.SetSelection(sel);
            var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            // Write each one
            // fill in the variable with values from starting index to
            // starting index + count
            engine.BeginStep();
            engine.Put(var_iString, currentTestData.S1);
            engine.Put(var_i8, currentTestData.I8.data());
            engine.Put(var_i16, currentTestData.I16.data());
            engine.Put(var_i32, currentTestData.I32.data());
            engine.Put(var_i64, currentTestData.I64.data());
            engine.Put(var_u8, currentTestData.U8.data());
            engine.Put(var_u16, currentTestData.U16.data());
            engine.Put(var_u32, currentTestData.U32.data());
            engine.Put(var_u64, currentTestData.U64.data());
            engine.Put(var_r32, currentTestData.R32.data());
            engine.Put(var_r64, currentTestData.R64.data());
            engine.EndStep();
        }

        engine.Close();
    }

    { // reading back
        adios2::IO io = adios.DeclareIO("HDF5ReadIO");
        io.SetEngine("HDF5");

        adios2::Engine hdf5Reader = io.Open(fname, adios2::Mode::Read);

        auto var_iString = io.InquireVariable<std::string>("iString");
        EXPECT_TRUE(var_iString);
        ASSERT_EQ(var_iString.Shape().size(), 0);
        ASSERT_EQ(var_iString.Steps(), NSteps);

        auto var_i8 = io.InquireVariable<int8_t>("i8");
        EXPECT_TRUE(var_i8);
        ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i8.Steps(), NSteps);
        ASSERT_EQ(var_i8.Shape()[0], Ny);
        ASSERT_EQ(var_i8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i16 = io.InquireVariable<int16_t>("i16");
        EXPECT_TRUE(var_i16);
        ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i16.Steps(), NSteps);
        ASSERT_EQ(var_i16.Shape()[0], Ny);
        ASSERT_EQ(var_i16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i32 = io.InquireVariable<int32_t>("i32");
        EXPECT_TRUE(var_i32);
        ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i32.Steps(), NSteps);
        ASSERT_EQ(var_i32.Shape()[0], Ny);
        ASSERT_EQ(var_i32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i64 = io.InquireVariable<int64_t>("i64");
        EXPECT_TRUE(var_i64);
        ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i64.Steps(), NSteps);
        ASSERT_EQ(var_i64.Shape()[0], Ny);
        ASSERT_EQ(var_i64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u8 = io.InquireVariable<uint8_t>("u8");
        EXPECT_TRUE(var_u8);
        ASSERT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u8.Steps(), NSteps);
        ASSERT_EQ(var_u8.Shape()[0], Ny);
        ASSERT_EQ(var_u8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u16 = io.InquireVariable<uint16_t>("u16");
        EXPECT_TRUE(var_u16);
        ASSERT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u16.Steps(), NSteps);
        ASSERT_EQ(var_u16.Shape()[0], Ny);
        ASSERT_EQ(var_u16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u32 = io.InquireVariable<uint32_t>("u32");
        EXPECT_TRUE(var_u32);
        ASSERT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u32.Steps(), NSteps);
        ASSERT_EQ(var_u32.Shape()[0], Ny);
        ASSERT_EQ(var_u32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u64 = io.InquireVariable<uint64_t>("u64");
        EXPECT_TRUE(var_u64);
        ASSERT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u64.Steps(), NSteps);
        ASSERT_EQ(var_u64.Shape()[0], Ny);
        ASSERT_EQ(var_u64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_r32 = io.InquireVariable<float>("r32");
        EXPECT_TRUE(var_r32);
        ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r32.Steps(), NSteps);
        ASSERT_EQ(var_r32.Shape()[0], Ny);
        ASSERT_EQ(var_r32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_r64 = io.InquireVariable<double>("r64");
        EXPECT_TRUE(var_r64);
        ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r64.Steps(), NSteps);
        ASSERT_EQ(var_r64.Shape()[0], Ny);
        ASSERT_EQ(var_r64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        std::string IString;
        std::array<int8_t, Nx * Ny> I8;
        std::array<int16_t, Nx * Ny> I16;
        std::array<int32_t, Nx * Ny> I32;
        std::array<int64_t, Nx * Ny> I64;
        std::array<uint8_t, Nx * Ny> U8;
        std::array<uint16_t, Nx * Ny> U16;
        std::array<uint32_t, Nx * Ny> U32;
        std::array<uint64_t, Nx * Ny> U64;
        std::array<float, Nx * Ny> R32;
        std::array<double, Nx * Ny> R64;

        const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
        const adios2::Dims count{Ny, Nx};

        const adios2::Box<adios2::Dims> sel(start, count);

        var_i8.SetSelection(sel);
        var_i16.SetSelection(sel);
        var_i32.SetSelection(sel);
        var_i64.SetSelection(sel);

        var_u8.SetSelection(sel);
        var_u16.SetSelection(sel);
        var_u32.SetSelection(sel);
        var_u64.SetSelection(sel);

        var_r32.SetSelection(sel);
        var_r64.SetSelection(sel);

        for (size_t t = 0; t < NSteps; ++t)
        {
            var_i8.SetStepSelection({t, 1});
            var_i16.SetStepSelection({t, 1});
            var_i32.SetStepSelection({t, 1});
            var_i64.SetStepSelection({t, 1});

            var_u8.SetStepSelection({t, 1});
            var_u16.SetStepSelection({t, 1});
            var_u32.SetStepSelection({t, 1});
            var_u64.SetStepSelection({t, 1});

            var_r32.SetStepSelection({t, 1});
            var_r64.SetStepSelection({t, 1});

            hdf5Reader.Get(var_iString, IString);

            hdf5Reader.Get(var_i8, I8.data());
            hdf5Reader.Get(var_i16, I16.data());
            hdf5Reader.Get(var_i32, I32.data());
            hdf5Reader.Get(var_i64, I64.data());

            hdf5Reader.Get(var_u8, U8.data());
            hdf5Reader.Get(var_u16, U16.data());
            hdf5Reader.Get(var_u32, U32.data());
            hdf5Reader.Get(var_u64, U64.data());

            hdf5Reader.Get(var_r32, R32.data());
            hdf5Reader.Get(var_r64, R64.data());
            hdf5Reader.PerformGets();

            // Generate test data for each rank uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(t), mpiRank, mpiSize);

            EXPECT_EQ(IString, currentTestData.S1);

            for (size_t i = 0; i < Nx * Ny; ++i)
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
        }
        hdf5Reader.Close();
    }
}

// Native HDF5 write, ADIOS2 read
TEST_F(HDF5WriteReadTest, HDF5WriteADIOS2HDF5Read2D4x2)
{
    std::string fname = "HDF5WriteADIOS2HDF5Read2D4x2Test.h5";

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const std::size_t Nx = 2;
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
#ifdef _WIN32
            h5writer.CreateAndStoreVar("i64", dimSize, H5T_NATIVE_LLONG, global_dims, offset, count,
                                       currentTestData.I64.data());
#else
            h5writer.CreateAndStoreVar("i64", dimSize, H5T_NATIVE_LONG, global_dims, offset, count,
                                       currentTestData.I64.data());
#endif
            h5writer.CreateAndStoreVar("u8", dimSize, H5T_NATIVE_UCHAR, global_dims, offset, count,
                                       currentTestData.U8.data());
            h5writer.CreateAndStoreVar("u16", dimSize, H5T_NATIVE_USHORT, global_dims, offset,
                                       count, currentTestData.U16.data());
            h5writer.CreateAndStoreVar("u32", dimSize, H5T_NATIVE_UINT, global_dims, offset, count,
                                       currentTestData.U32.data());
#ifdef _WIN32
            h5writer.CreateAndStoreVar("u64", dimSize, H5T_NATIVE_ULLONG, global_dims, offset,
                                       count, currentTestData.U64.data());
#else
            h5writer.CreateAndStoreVar("u64", dimSize, H5T_NATIVE_ULONG, global_dims, offset, count,
                                       currentTestData.U64.data());
#endif
            h5writer.CreateAndStoreVar("r32", dimSize, H5T_NATIVE_FLOAT, global_dims, offset, count,
                                       currentTestData.R32.data());
            h5writer.CreateAndStoreVar("r64", dimSize, H5T_NATIVE_DOUBLE, global_dims, offset,
                                       count, currentTestData.R64.data());
            h5writer.Advance();
        }
    }

    { // read back
#ifdef TEST_HDF5_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

#ifdef TEST_HDF5_MPI
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif

        adios2::IO io = adios.DeclareIO("HDF5ReadIO");
        io.SetEngine("HDF5");

        adios2::Engine hdf5Reader = io.Open(fname, adios2::Mode::Read);

        auto var_iString = io.InquireVariable<std::string>("iString");
        EXPECT_TRUE(var_iString);
        ASSERT_EQ(var_iString.Shape().size(), 0);
        ASSERT_EQ(var_iString.Steps(), NSteps);

        auto var_i8 = io.InquireVariable<int8_t>("i8");
        EXPECT_TRUE(var_i8);
        ASSERT_EQ(var_i8.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i8.Steps(), NSteps);
        ASSERT_EQ(var_i8.Shape()[0], Ny);
        ASSERT_EQ(var_i8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i16 = io.InquireVariable<int16_t>("i16");
        EXPECT_TRUE(var_i16);
        ASSERT_EQ(var_i16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i16.Steps(), NSteps);
        ASSERT_EQ(var_i16.Shape()[0], Ny);
        ASSERT_EQ(var_i16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i32 = io.InquireVariable<int32_t>("i32");
        EXPECT_TRUE(var_i32);
        ASSERT_EQ(var_i32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i32.Steps(), NSteps);
        ASSERT_EQ(var_i32.Shape()[0], Ny);
        ASSERT_EQ(var_i32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_i64 = io.InquireVariable<int64_t>("i64");
        EXPECT_TRUE(var_i64);
        ASSERT_EQ(var_i64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_i64.Steps(), NSteps);
        ASSERT_EQ(var_i64.Shape()[0], Ny);
        ASSERT_EQ(var_i64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u8 = io.InquireVariable<uint8_t>("u8");
        EXPECT_TRUE(var_u8);
        ASSERT_EQ(var_u8.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u8.Steps(), NSteps);
        ASSERT_EQ(var_u8.Shape()[0], Ny);
        ASSERT_EQ(var_u8.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u16 = io.InquireVariable<uint16_t>("u16");
        EXPECT_TRUE(var_u16);
        ASSERT_EQ(var_u16.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u16.Steps(), NSteps);
        ASSERT_EQ(var_u16.Shape()[0], Ny);
        ASSERT_EQ(var_u16.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u32 = io.InquireVariable<uint32_t>("u32");
        EXPECT_TRUE(var_u32);
        ASSERT_EQ(var_u32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u32.Steps(), NSteps);
        ASSERT_EQ(var_u32.Shape()[0], Ny);
        ASSERT_EQ(var_u32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_u64 = io.InquireVariable<uint64_t>("u64");
        EXPECT_TRUE(var_u64);
        ASSERT_EQ(var_u64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_u64.Steps(), NSteps);
        ASSERT_EQ(var_u64.Shape()[0], Ny);
        ASSERT_EQ(var_u64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_r32 = io.InquireVariable<float>("r32");
        EXPECT_TRUE(var_r32);
        ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r32.Steps(), NSteps);
        ASSERT_EQ(var_r32.Shape()[0], Ny);
        ASSERT_EQ(var_r32.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        auto var_r64 = io.InquireVariable<double>("r64");
        EXPECT_TRUE(var_r64);
        ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var_r64.Steps(), NSteps);
        ASSERT_EQ(var_r64.Shape()[0], Ny);
        ASSERT_EQ(var_r64.Shape()[1], static_cast<size_t>(mpiSize * Nx));

        std::string IString;
        std::array<int8_t, Nx * Ny> I8;
        std::array<int16_t, Nx * Ny> I16;
        std::array<int32_t, Nx * Ny> I32;
        std::array<int64_t, Nx * Ny> I64;
        std::array<uint8_t, Nx * Ny> U8;
        std::array<uint16_t, Nx * Ny> U16;
        std::array<uint32_t, Nx * Ny> U32;
        std::array<uint64_t, Nx * Ny> U64;
        std::array<float, Nx * Ny> R32;
        std::array<double, Nx * Ny> R64;

        const adios2::Dims start{0, static_cast<size_t>(mpiRank * Nx)};
        const adios2::Dims count{Ny, Nx};

        const adios2::Box<adios2::Dims> sel(start, count);

        var_i8.SetSelection(sel);
        var_i16.SetSelection(sel);
        var_i32.SetSelection(sel);
        var_i64.SetSelection(sel);

        var_u8.SetSelection(sel);
        var_u16.SetSelection(sel);
        var_u32.SetSelection(sel);
        var_u64.SetSelection(sel);

        var_r32.SetSelection(sel);
        var_r64.SetSelection(sel);

        for (size_t t = 0; t < NSteps; ++t)
        {
            var_i8.SetStepSelection({t, 1});
            var_i16.SetStepSelection({t, 1});
            var_i32.SetStepSelection({t, 1});
            var_i64.SetStepSelection({t, 1});

            var_u8.SetStepSelection({t, 1});
            var_u16.SetStepSelection({t, 1});
            var_u32.SetStepSelection({t, 1});
            var_u64.SetStepSelection({t, 1});

            var_r32.SetStepSelection({t, 1});
            var_r64.SetStepSelection({t, 1});

            hdf5Reader.Get(var_iString, IString);

            hdf5Reader.Get(var_i8, I8.data());
            hdf5Reader.Get(var_i16, I16.data());
            hdf5Reader.Get(var_i32, I32.data());
            hdf5Reader.Get(var_i64, I64.data());

            hdf5Reader.Get(var_u8, U8.data());
            hdf5Reader.Get(var_u16, U16.data());
            hdf5Reader.Get(var_u32, U32.data());
            hdf5Reader.Get(var_u64, U64.data());

            hdf5Reader.Get(var_r32, R32.data());
            hdf5Reader.Get(var_r64, R64.data());
            hdf5Reader.PerformGets();

            // Generate test data for each rank uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(t), mpiRank, mpiSize);

            EXPECT_EQ(IString, currentTestData.S1);

            for (size_t i = 0; i < Nx * Ny; ++i)
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
        }
        hdf5Reader.Close();
    }
}

TEST_F(HDF5WriteReadTest, /*DISABLE_*/ ATTRTESTADIOS2vsHDF5)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    const std::string h5name = "ATTRTESTADIOS2vsHDF5.h5";
    const std::string bpname = "ATTRTESTADIOS2vsHDF5.bp";

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const std::size_t Nx = 8;

    // Number of steps
    const std::size_t NSteps = 3;

#ifdef TEST_HDF5_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // Write test data using ADIOS2

#ifdef TEST_HDF5_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    adios2::IO io = adios.DeclareIO("TestIO");

    std::string var1Name = "var1";
    // std::string var1Name = "/var1";
    std::string var2Name = "var2";
    std::string var3Name = "grp/var3";
    std::string var4Name = "var4";

    std::string ioAttrName = "file_notes";
    std::string ioAttrValue = "this is a comment";

    std::string v1AttrAName = "var1/meshWidth";
    int32_t v1AttrAValue = 111;
    std::string v1AttrBName = "var1/meshLength";
    int32_t v1AttrBValue = 222;

    std::string v2AttrName = "var2/meshName";
    std::string v2AttrValue = "unstructured";

    std::string v3AttrName = var3Name + "/targetDensity";
    double v3AttrValue = 55;

    std::string v4AttrName = "var4/unitName";
    std::string v4AttrValue = "kg";

    std::string fantacyPathName = "no/such/path";
    std::string longerV1AttrName1 = "var1/mesh/level";
    std::string longerV1AttrName2 = "mesh/also/level";

    int numAttr = 0;
    // Declare 1D variables (NumOfProcesses * Nx)
    // The local process' part (start, count) can be defined now or later
    // before Write().
    {
        adios2::Dims shape{static_cast<unsigned int>(Nx * mpiSize)};
        adios2::Dims start{static_cast<unsigned int>(Nx * mpiRank)};
        adios2::Dims count{static_cast<unsigned int>(Nx)};

        io.DefineVariable<int32_t>(var1Name, shape, start, count);
        io.DefineVariable<float>(var2Name, shape, start, count);
        io.DefineVariable<double>(var3Name, shape, start, count);

        io.DefineAttribute<std::string>(ioAttrName, ioAttrValue);
        numAttr++;

        io.DefineAttribute<int32_t>(v1AttrAName, v1AttrAValue);
        numAttr++;
        io.DefineAttribute<int32_t>(v1AttrBName, v1AttrBValue);
        numAttr++;

        // this expects that  var and attr can have different names
        io.DefineAttribute<double>(var1Name, v3AttrValue);
        numAttr++;

        // this allows attr can be path like (no need to exist)
        io.DefineAttribute<int32_t>(fantacyPathName, v1AttrAValue);
        numAttr++;

        // this allows variable can attach attr with long name
        io.DefineAttribute<int32_t>(longerV1AttrName1, v1AttrBValue);
        numAttr++;
        io.DefineAttribute<int32_t>(longerV1AttrName2, v1AttrBValue, var1Name);
        numAttr++;

        io.DefineAttribute<std::string>(v2AttrName, v2AttrValue);
        numAttr++;
        io.DefineAttribute<double>(v3AttrName, v3AttrValue);
        numAttr++;
    }

    // Create the HDF5 Engine
    io.SetEngine("HDF5");

    // HDf5 engine calls the HDF5 common object that calls the hDF5 library.
    // The IO functionality, SetParameters and AddTransports will be added
    // in the future. For now `io.AddTransport("file", {
    // "library", "MPI"}});` is omitted.
    // })
    // io.AddTransport("File");

    adios2::Engine engine = io.Open(h5name, adios2::Mode::Write);

    for (size_t step = 0; step < NSteps; ++step)
    {
        // Generate test data for each process uniquely
        SmallTestData currentTestData =
            generateNewSmallTestData(m_TestData, step, mpiRank, mpiSize);

        // Retrieve the variables that previously went out of scope
        auto var1 = io.InquireVariable<int32_t>(var1Name);
        auto var2 = io.InquireVariable<float>(var2Name);
        auto var3 = io.InquireVariable<double>(var3Name);

        // Make a 1D selection to describe the local dimensions of the
        // variable we write and its offsets in the global spaces

        adios2::Box<adios2::Dims> sel({mpiRank * Nx}, {Nx});

        var1.SetSelection(sel);
        var2.SetSelection(sel);
        var3.SetSelection(sel);

        // Write each one
        // fill in the variable with values from starting index to
        // starting index + count
        engine.BeginStep();
        engine.Put(var1, currentTestData.I32.data());
        engine.Put(var2, currentTestData.R32.data());
        engine.Put(var3, currentTestData.R64.data());
        if (step == 1)
        {
            adios2::Dims shape{static_cast<unsigned int>(Nx * mpiSize)};
            adios2::Dims start{static_cast<unsigned int>(Nx * mpiRank)};
            adios2::Dims count{static_cast<unsigned int>(Nx)};

            io.DefineVariable<int32_t>(var4Name, shape, start, count);
            auto var4 = io.InquireVariable<int32_t>(var4Name);
            engine.Put(var4, currentTestData.I32.data());

            io.DefineAttribute<std::string>(v4AttrName, v4AttrValue);
            numAttr++;
        }

        // Advance to the next time step
        engine.EndStep();

        if (step == 1)
        {
            io.RemoveVariable(var4Name);
            io.RemoveAttribute(v4AttrName);
        }
    }

    // Close the file
    engine.Close();

    {
        adios2::IO io = adios.DeclareIO("HDF5ReadIO");
        io.SetEngine("HDF5");

        adios2::Engine hdf5Reader = io.Open(h5name, adios2::Mode::Read);

        auto var1 = io.InquireVariable<int32_t>(var1Name);
        EXPECT_TRUE(var1);
        ASSERT_EQ(var1.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var1.Steps(), NSteps);
        ASSERT_EQ(var1.Shape()[0], mpiSize * Nx);

        auto var2 = io.InquireVariable<float>(var2Name);
        EXPECT_TRUE(var2);
        ASSERT_EQ(var2.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var2.Steps(), NSteps);
        ASSERT_EQ(var2.Shape()[0], mpiSize * Nx);

        auto var3 = io.InquireVariable<double>(var3Name);
        EXPECT_TRUE(var3);
        ASSERT_EQ(var3.ShapeID(), adios2::ShapeID::GlobalArray);
        ASSERT_EQ(var3.Steps(), NSteps);
        ASSERT_EQ(var3.Shape()[0], mpiSize * Nx);

        SmallTestData testData;

        std::array<int32_t, Nx> I32;
        std::array<float, Nx> R32;
        std::array<double, Nx> R64;

        const adios2::Dims start{mpiRank * Nx};
        const adios2::Dims count{Nx};

        const adios2::Box<adios2::Dims> sel(start, count);

        var1.SetSelection(sel);
        var2.SetSelection(sel);
        var3.SetSelection(sel);

        for (size_t t = 0; t < NSteps; ++t)
        {
            var1.SetStepSelection({t, 1});
            var2.SetStepSelection({t, 1});
            var3.SetStepSelection({t, 1});

            // Generate test data for each rank uniquely
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(t), mpiRank, mpiSize);

            hdf5Reader.BeginStep();

            hdf5Reader.Get(var1, I32.data());
            hdf5Reader.Get(var2, R32.data());
            hdf5Reader.Get(var3, R64.data());

            hdf5Reader.EndStep();

            // EXPECT_EQ(IString, currentTestData.S1);

            for (size_t i = 0; i < Nx; ++i)
            {
                std::stringstream ss;
                ss << "t=" << t << " i=" << i << " rank=" << mpiRank;
                std::string msg = ss.str();

                EXPECT_EQ(I32[i], currentTestData.I32[i]) << msg;
                EXPECT_EQ(R32[i], currentTestData.R32[i]) << msg;
                EXPECT_EQ(R64[i], currentTestData.R64[i]) << msg;
            }
        }

        const std::map<std::string, adios2::Params> &attributesInfo = io.AvailableAttributes();

        EXPECT_EQ(numAttr, attributesInfo.size());

        EXPECT_EQ(4, io.AvailableAttributes(var1Name).size());
        EXPECT_EQ(1, io.AvailableAttributes(var2Name).size());
        EXPECT_EQ(1, io.AvailableAttributes(var3Name).size());

        std::cout << " test below of attrs of var4Name is false b/c "
                     "io.m_EngineStep=2 but var4 is only in step 1, this it "
                     "thinks so such var"
                  << std::endl;
        std::cout << " need to fix semantics of io.AvailableAttributes()" << std::endl;
        // EXPECT_EQ(1, io.AvailableAttributes(var4Name).size());

        std::cout << " other tests will follow after William make changes: "
                     "e.g. GetNumAttr(var) etc +  Write a bp file and read "
                     "back in hdf5 and verify"
                  << std::endl;
        /*
        hdf5Reader.GetNumAttributes(var1); // expects 2
        hdf5Reader.GetNumAttributes(var2); // expects 1
        hdf5Reader.GetNumAttributes(var3); // expects 0

        hdf5Reader.GetNumAttributes(); // expects 4
        */
        hdf5Reader.Close();
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

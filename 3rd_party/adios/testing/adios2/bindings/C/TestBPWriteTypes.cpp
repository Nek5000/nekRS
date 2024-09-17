/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * TestBPWriteTypes.cpp : test the C bindings
 *
 *  Created on: Aug 9, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include <cstring>

#include <algorithm>

#include <adios2_c.h>

#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

#include <gtest/gtest.h>

#include "SmallTestData_c.h"

class ADIOS2_C_API : public ::testing::Test
{
public:
    ADIOS2_C_API()
    {
#if ADIOS2_USE_MPI
        adiosH = adios2_init_mpi(MPI_COMM_WORLD);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
        adiosH = adios2_init_serial();
#endif
    }

    ~ADIOS2_C_API() { adios2_finalize(adiosH); }

    adios2_adios *adiosH;
    int rank = 0;
    int size = 1;
};

TEST_F(ADIOS2_C_API, ADIOS2BPWriteTypes)
{
#if ADIOS2_USE_MPI
    const char fname[] = "ADIOS2_C_API.ADIOS2BPWriteTypes_MPI.bp";
#else
    const char fname[] = "ADIOS2_C_API.ADIOS2BPWriteTypes.bp";
#endif
    // write
    {
        // IO
        adios2_io *ioH = adios2_declare_io(adiosH, "CArrayTypes");
        // Set engine parameters
        adios2_set_engine(ioH, "BPFile");
        adios2_set_parameter(ioH, "ProfileUnits", "Microseconds");
        adios2_set_parameters(ioH, "Threads=2, CollectiveMetadata = OFF");

        size_t length;

        char *profileUnits;
        adios2_get_parameter(NULL, &length, ioH, "ProfileUnits");
        profileUnits = (char *)calloc(length + 1, sizeof(char));
        adios2_get_parameter(profileUnits, &length, ioH, "ProfileUnits");
        EXPECT_EQ(std::string(profileUnits), "Microseconds");
        free(profileUnits);

        char *threads;
        adios2_get_parameter(NULL, &length, ioH, "Threads");
        threads = (char *)calloc(length + 1, sizeof(char));
        adios2_get_parameter(threads, &length, ioH, "Threads");
        EXPECT_EQ(std::string(threads), "2");
        free(threads);

        char *collmd;
        adios2_get_parameter(NULL, &length, ioH, "CollectiveMetadata");
        collmd = (char *)calloc(length + 1, sizeof(char));
        adios2_get_parameter(collmd, &length, ioH, "CollectiveMetadata");
        EXPECT_EQ(std::string(collmd), "OFF");
        free(collmd);

        // set back the default to make sure writing/reading works
        adios2_clear_parameters(ioH);
        adios2_get_parameter(NULL, &length, ioH, "CollectiveMetadata");
        EXPECT_EQ(length, 0);

        // Set transport and parameters
        size_t transportID;
        adios2_add_transport(&transportID, ioH, "File");
        adios2_set_transport_parameter(ioH, transportID, "library", "fstream");

        // count dims are allocated in stack
        size_t shape[1];
        shape[0] = data_Nx * static_cast<size_t>(size);

        size_t start[1];
        start[0] = data_Nx * static_cast<size_t>(rank);

        size_t count[1];
        count[0] = data_Nx;

        // define attribute
        const char *strarray[] = {"first", "second", "third", "fourth"};

        adios2_define_attribute_array(ioH, "strarray", adios2_type_string, strarray, 4);

        adios2_define_attribute(ioH, "strvalue", adios2_type_string, "Hello Attribute");

        // Define variables in ioH
        {
            adios2_define_variable(ioH, "varStr", adios2_type_string, 0, NULL, NULL, NULL,
                                   adios2_constant_dims_true);
            adios2_define_variable(ioH, "varI8", adios2_type_int8_t, 1, shape, start, count,
                                   adios2_constant_dims_true);
            adios2_define_variable(ioH, "varI16", adios2_type_int16_t, 1, shape, start, count,
                                   adios2_constant_dims_true);
            adios2_define_variable(ioH, "varI32", adios2_type_int32_t, 1, shape, start, count,
                                   adios2_constant_dims_true);
            adios2_define_variable(ioH, "varI64", adios2_type_int64_t, 1, shape, start, count,
                                   adios2_constant_dims_true);

            adios2_define_variable(ioH, "varU8", adios2_type_uint8_t, 1, shape, start, count,
                                   adios2_constant_dims_true);
            adios2_define_variable(ioH, "varU16", adios2_type_uint16_t, 1, shape, start, count,
                                   adios2_constant_dims_true);
            adios2_define_variable(ioH, "varU32", adios2_type_uint32_t, 1, shape, start, count,
                                   adios2_constant_dims_true);
            adios2_define_variable(ioH, "varU64", adios2_type_uint64_t, 1, shape, start, count,
                                   adios2_constant_dims_true);

            adios2_define_variable(ioH, "varR32", adios2_type_float, 1, shape, start, count,
                                   adios2_constant_dims_true);
            adios2_define_variable(ioH, "varR64", adios2_type_double, 1, shape, start, count,
                                   adios2_constant_dims_true);
        }
        // inquire variables
        adios2_variable *varStr = adios2_inquire_variable(ioH, "varStr");
        adios2_variable *varI8 = adios2_inquire_variable(ioH, "varI8");
        adios2_variable *varI16 = adios2_inquire_variable(ioH, "varI16");
        adios2_variable *varI32 = adios2_inquire_variable(ioH, "varI32");
        adios2_variable *varI64 = adios2_inquire_variable(ioH, "varI64");

        adios2_variable *varU8 = adios2_inquire_variable(ioH, "varU8");
        adios2_variable *varU16 = adios2_inquire_variable(ioH, "varU16");
        adios2_variable *varU32 = adios2_inquire_variable(ioH, "varU32");
        adios2_variable *varU64 = adios2_inquire_variable(ioH, "varU64");

        adios2_variable *varR32 = adios2_inquire_variable(ioH, "varR32");
        adios2_variable *varR64 = adios2_inquire_variable(ioH, "varR64");

        adios2_engine *engineH = adios2_open(ioH, fname, adios2_mode_write);

        adios2_mode m;
        adios2_error e = adios2_engine_openmode(&m, engineH);
        EXPECT_EQ(e, adios2_error_none);
        EXPECT_EQ(m, adios2_mode_write);

        adios2_put(engineH, varStr, dataStr, adios2_mode_sync);
        adios2_put(engineH, varI8, data_I8, adios2_mode_deferred);
        adios2_put(engineH, varI8, data_I8, adios2_mode_deferred);
        adios2_put(engineH, varI16, data_I16, adios2_mode_deferred);
        adios2_put(engineH, varI32, data_I32, adios2_mode_deferred);
        adios2_put(engineH, varI64, data_I64, adios2_mode_deferred);

        adios2_put(engineH, varU8, data_U8, adios2_mode_deferred);
        adios2_put(engineH, varU16, data_U16, adios2_mode_deferred);
        adios2_put(engineH, varU32, data_U32, adios2_mode_deferred);
        adios2_put(engineH, varU64, data_U64, adios2_mode_deferred);

        adios2_put(engineH, varR32, data_R32, adios2_mode_deferred);
        adios2_put(engineH, varR64, data_R64, adios2_mode_deferred);

        size_t engineTypeSize;
        adios2_engine_get_type(NULL, &engineTypeSize, engineH);
        EXPECT_EQ(engineTypeSize, 9);

        char *engineType = new char[engineTypeSize + 1]();
        adios2_engine_get_type(engineType, &engineTypeSize, engineH);

        EXPECT_EQ(std::string(engineType, engineTypeSize), "BP5Writer");

        adios2_close(engineH);

        delete[] engineType;
    }
#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    {
        adios2_io *ioH = adios2_declare_io(adiosH, "Reader");
        adios2_engine *engineH = adios2_open(ioH, fname, adios2_mode_readRandomAccess);

        size_t steps;
        adios2_steps(&steps, engineH);
        EXPECT_EQ(steps, 1);

        adios2_mode m;
        adios2_error e = adios2_engine_openmode(&m, engineH);
        EXPECT_EQ(e, adios2_error_none);
        EXPECT_EQ(m, adios2_mode_readRandomAccess);

        adios2_bool result;
        char name[30];
        size_t name_size;
        adios2_type type;
        char type_str[30];
        size_t type_str_size;

        // single
        adios2_attribute *attrSingle = adios2_inquire_attribute(ioH, "strvalue");
        adios2_attribute_is_value(&result, attrSingle);
        EXPECT_EQ(result, adios2_true);

        adios2_attribute_name(name, &name_size, attrSingle);
        const std::string nameStrSingle(name, name_size);
        EXPECT_EQ(nameStrSingle, "strvalue");

        adios2_attribute_type(&type, attrSingle);
        EXPECT_EQ(type, adios2_type_string);

        adios2_attribute_type_string(type_str, &type_str_size, attrSingle);
        const std::string typeStrSingle(type_str, type_str_size);
        EXPECT_EQ(typeStrSingle, "string");

        char *dataSingle = new char[30]();
        size_t elements;
        adios2_attribute_data(dataSingle, &elements, attrSingle);
        const std::string dataSingleStr(dataSingle);
        EXPECT_EQ(dataSingleStr, "Hello Attribute");
        EXPECT_EQ(elements, 1);
        delete[] dataSingle;

        // Array
        adios2_attribute *attrArray = adios2_inquire_attribute(ioH, "strarray");
        adios2_attribute_is_value(&result, attrArray);
        EXPECT_EQ(result, adios2_false);

        adios2_attribute_name(name, &name_size, attrArray);
        const std::string nameStrArray(name, name_size);
        EXPECT_EQ(nameStrArray, "strarray");

        adios2_attribute_type(&type, attrArray);
        EXPECT_EQ(type, adios2_type_string);

        adios2_attribute_type_string(type_str, &type_str_size, attrArray);
        const std::string typeStrArray(type_str, type_str_size);
        EXPECT_EQ(typeStrArray, "string");

        char **dataArray = new char *[4];
        for (int i = 0; i < 4; ++i)
        {
            dataArray[i] = new char[30]();
        }

        adios2_attribute_data(dataArray, &elements, attrArray);
        const std::vector<std::string> dataVector(dataArray, dataArray + elements);
        EXPECT_EQ(dataVector[0], "first");
        EXPECT_EQ(dataVector[1], "second");
        EXPECT_EQ(dataVector[2], "third");
        EXPECT_EQ(dataVector[3], "fourth");
        EXPECT_EQ(elements, 4);

        // compare min and max
        auto mmI8 = std::minmax_element(&data_I8[0], &data_I8[data_Nx]);
        auto mmI16 = std::minmax_element(&data_I16[0], &data_I16[data_Nx]);
        auto mmI32 = std::minmax_element(&data_I32[0], &data_I32[data_Nx]);
        auto mmI64 = std::minmax_element(&data_I64[0], &data_I64[data_Nx]);
        auto mmU8 = std::minmax_element(&data_U8[0], &data_U8[data_Nx]);
        auto mmU16 = std::minmax_element(&data_U16[0], &data_U16[data_Nx]);
        auto mmU32 = std::minmax_element(&data_U32[0], &data_U32[data_Nx]);
        auto mmU64 = std::minmax_element(&data_U64[0], &data_U64[data_Nx]);
        auto mmR32 = std::minmax_element(&data_R32[0], &data_R32[data_Nx]);
        auto mmR64 = std::minmax_element(&data_R64[0], &data_R64[data_Nx]);

        // add min and max here
        adios2_variable *varStr = adios2_inquire_variable(ioH, "varStr");
        adios2_variable *varI8 = adios2_inquire_variable(ioH, "varI8");
        adios2_variable *varI16 = adios2_inquire_variable(ioH, "varI16");
        adios2_variable *varI32 = adios2_inquire_variable(ioH, "varI32");
        adios2_variable *varI64 = adios2_inquire_variable(ioH, "varI64");

        adios2_variable *varU8 = adios2_inquire_variable(ioH, "varU8");
        adios2_variable *varU16 = adios2_inquire_variable(ioH, "varU16");
        adios2_variable *varU32 = adios2_inquire_variable(ioH, "varU32");
        adios2_variable *varU64 = adios2_inquire_variable(ioH, "varU64");

        adios2_variable *varR32 = adios2_inquire_variable(ioH, "varR32");
        adios2_variable *varR64 = adios2_inquire_variable(ioH, "varR64");

        char inString[30] = {};
        adios2_get(engineH, varStr, inString, adios2_mode_sync);

        EXPECT_EQ(dataStr, std::string(inString));

        int8_t minI8, maxI8;
        int16_t minI16, maxI16;
        int32_t minI32, maxI32;
        int64_t minI64, maxI64;

        uint8_t minU8, maxU8;
        uint16_t minU16, maxU16;
        uint32_t minU32, maxU32;
        uint64_t minU64, maxU64;

        float minR32, maxR32;
        double minR64, maxR64;

        adios2_variable_min(&minI8, varI8);
        adios2_variable_min(&minI16, varI16);
        adios2_variable_min(&minI32, varI32);
        adios2_variable_min(&minI64, varI64);
        adios2_variable_min(&minU8, varU8);
        adios2_variable_min(&minU16, varU16);
        adios2_variable_min(&minU32, varU32);
        adios2_variable_min(&minU64, varU64);
        adios2_variable_min(&minR32, varR32);
        adios2_variable_min(&minR64, varR64);

        adios2_variable_max(&maxI8, varI8);
        adios2_variable_max(&maxI16, varI16);
        adios2_variable_max(&maxI32, varI32);
        adios2_variable_max(&maxI64, varI64);
        adios2_variable_max(&maxU8, varU8);
        adios2_variable_max(&maxU16, varU16);
        adios2_variable_max(&maxU32, varU32);
        adios2_variable_max(&maxU64, varU64);
        adios2_variable_max(&maxR32, varR32);
        adios2_variable_max(&maxR64, varR64);

        EXPECT_EQ(minI8, *mmI8.first);
        EXPECT_EQ(minI16, *mmI16.first);
        EXPECT_EQ(minI32, *mmI32.first);
        EXPECT_EQ(minI64, *mmI64.first);
        EXPECT_EQ(minU8, *mmU8.first);
        EXPECT_EQ(minU16, *mmU16.first);
        EXPECT_EQ(minU32, *mmU32.first);
        EXPECT_EQ(minU64, *mmU64.first);
        EXPECT_EQ(minR32, *mmR32.first);
        EXPECT_EQ(minR64, *mmR64.first);

        EXPECT_EQ(maxI8, *mmI8.second);
        EXPECT_EQ(maxI16, *mmI16.second);
        EXPECT_EQ(maxI32, *mmI32.second);
        EXPECT_EQ(maxI64, *mmI64.second);
        EXPECT_EQ(maxU8, *mmU8.second);
        EXPECT_EQ(maxU16, *mmU16.second);
        EXPECT_EQ(maxU32, *mmU32.second);
        EXPECT_EQ(maxU64, *mmU64.second);
        EXPECT_EQ(maxR32, *mmR32.second);
        EXPECT_EQ(maxR64, *mmR64.second);

        adios2_close(engineH);

        for (size_t i = 0; i < 4; ++i)
        {
            delete[] dataArray[i];
        }
        delete[] dataArray;
    }
}

class ADIOS2_C_API_IO : public ADIOS2_C_API
{
public:
    ADIOS2_C_API_IO() { ioH = adios2_declare_io(adiosH, "C_API_TestIO"); }

    adios2_io *ioH;
};

namespace testing // helpers for testing
{
/**
 * helper that tests the C API for adios2_engine_type,
 * and then converts the result to a std::string for easier
 * handling
 */
std::string adios2_engine_type_as_string(adios2_io *ioH)
{
    int ierr;
    size_t engine_type_size;
    ierr = adios2_engine_type(NULL, &engine_type_size, ioH);
    EXPECT_EQ(ierr, 0);
    char *engine_type = (char *)malloc(engine_type_size + 1);
    ierr = adios2_engine_type(engine_type, &engine_type_size, ioH);
    EXPECT_EQ(ierr, 0);
    engine_type[engine_type_size] = '\0';
    std::string type(engine_type);
    free(engine_type);

    return type;
}

std::string adios2_engine_name_as_string(adios2_engine *engineH)
{
    int ierr;
    size_t engine_name_size;
    ierr = adios2_engine_name(NULL, &engine_name_size, engineH);
    EXPECT_EQ(ierr, 0);
    char *engine_name = (char *)malloc(engine_name_size + 1);
    ierr = adios2_engine_name(engine_name, &engine_name_size, engineH);
    EXPECT_EQ(ierr, 0);
    engine_name[engine_name_size] = '\0';
    std::string name(engine_name);
    free(engine_name);

    return name;
}
}

TEST_F(ADIOS2_C_API_IO, Engine)
{
#if ADIOS2_USE_MPI
    const char fname[] = "ADIOS2_C_API_IO.engine_MPI.bp";
#else
    const char fname[] = "ADIOS2_C_API_IO.engine.bp";
#endif
    int ierr;

    ierr = adios2_set_engine(ioH, "bpfile");
    EXPECT_EQ(ierr, 0);

    std::string engine_type = testing::adios2_engine_type_as_string(ioH);
    EXPECT_EQ(engine_type, "bpfile");

    adios2_engine *engineH = adios2_open(ioH, fname, adios2_mode_write);

    // FIXME, I'd like to check that the engine type itself is correct, but
    // there's no API to get it
    std::string engine_name = testing::adios2_engine_name_as_string(engineH);
    EXPECT_EQ(engine_name, fname);
}

TEST_F(ADIOS2_C_API_IO, EngineDefault)
{
    const char fname[] = "ADIOS2_C_API_IO.EngineDefault.bp";
    int ierr;

    ierr = adios2_set_engine(ioH, "");
    EXPECT_EQ(ierr, 0);

    std::string engine_type = testing::adios2_engine_type_as_string(ioH);
    EXPECT_EQ(engine_type, "");

    adios2_engine *engineH = adios2_open(ioH, fname, adios2_mode_write);

    // FIXME, I'd like to check that the engine type itself is correct, but
    // there's no API to get it
    std::string engine_name = testing::adios2_engine_name_as_string(engineH);
    EXPECT_EQ(engine_name, fname);
}

TEST_F(ADIOS2_C_API_IO, ReturnedStrings)
{
    // create some objects
    adios2_set_engine(ioH, "BP3");

    size_t shape[1] = {2};
    size_t start[1] = {0};
    size_t count[1] = {2};
    adios2_variable *var = adios2_define_variable(ioH, "varI8", adios2_type_int8_t, 1, shape, start,
                                                  count, adios2_constant_dims_true);

    int32_t value = 99;
    adios2_attribute *attr = adios2_define_attribute(ioH, "intAttr", adios2_type_int32_t, &value);

#ifdef ADIOS2_HAVE_BZIP2
    adios2_operator *op = adios2_define_operator(adiosH, "testOp", "bzip2");
#endif

    // now test the APIs that return strings
    std::string engine_type = testing::adios2_engine_type_as_string(ioH);

    size_t var_name_size;
    adios2_variable_name(NULL, &var_name_size, var);
    char *var_name = (char *)malloc(var_name_size + 1);
    adios2_variable_name(var_name, &var_name_size, var);
    var_name[var_name_size] = '\0';

    size_t attr_name_size;
    adios2_attribute_name(NULL, &attr_name_size, attr);
    char *attr_name = (char *)malloc(attr_name_size + 1);
    adios2_attribute_name(attr_name, &attr_name_size, attr);
    attr_name[attr_name_size] = '\0';

#ifdef ADIOS2_HAVE_BZIP2
    size_t op_type_size;
    adios2_operator_type(NULL, &op_type_size, op);
    char *op_type = (char *)malloc(op_type_size + 1);
    adios2_operator_type(op_type, &op_type_size, op);
    op_type[op_type_size] = '\0';
#endif

    // remove all IOs here to specifically check that the strings
    // are still accessible (which is kinda obvious)
    adios2_remove_all_ios(adiosH);

    EXPECT_EQ(engine_type, "BP3");
    EXPECT_EQ(std::string(var_name), "varI8");
    EXPECT_EQ(std::string(attr_name), "intAttr");
#ifdef ADIOS2_HAVE_BZIP2
    EXPECT_EQ(std::string(op_type), "bzip2");
#endif

    free(var_name);
    free(attr_name);
#ifdef ADIOS2_HAVE_BZIP2
    free(op_type);
#endif
}

//******************************************************************************
// main
//******************************************************************************

int main(int argc, char **argv)
{
#if ADIOS2_USE_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(nullptr, nullptr, MPI_THREAD_MULTIPLE, &provided);
#endif

    int result;
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return result;
}

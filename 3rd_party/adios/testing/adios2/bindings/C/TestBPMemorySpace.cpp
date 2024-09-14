/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * TestBPMemorySpace.cpp : test the C bindings related to memory space
 */

#include <adios2_c.h>
#include <gtest/gtest.h>

class ADIOS2_C_API : public ::testing::Test
{
public:
    ADIOS2_C_API() { adiosH = adios2_init_serial(); }

    ~ADIOS2_C_API() { adios2_finalize(adiosH); }

    adios2_adios *adiosH;
};

#ifdef ADIOS2_HAVE_GPU_SUPPORT
TEST_F(ADIOS2_C_API, ADIOS2BPMemorySpaceGPU)
{
    adios2_io *ioH = adios2_declare_io(adiosH, "CMemSpace");
    adios2_set_engine(ioH, "BPFile");

    size_t shape[1];
    shape[0] = 10;

    size_t start[1];
    start[0] = 0;

    size_t count[1];
    count[0] = 10;

    adios2_define_variable(ioH, "varI32", adios2_type_int32_t, 1, shape, start, count,
                           adios2_constant_dims_true);
    adios2_variable *varI32 = adios2_inquire_variable(ioH, "varI32");

    // test that the default memory space is Detect
    adios2_memory_space mem;
    adios2_error e = adios2_get_memory_space(&mem, varI32);
    EXPECT_EQ(e, adios2_error_none);
    EXPECT_EQ(mem, adios2_memory_space_detect);
    // test that the set memory space is GPU
    adios2_set_memory_space(varI32, adios2_memory_space_gpu);
    e = adios2_get_memory_space(&mem, varI32);
    EXPECT_EQ(e, adios2_error_none);
    EXPECT_EQ(mem, adios2_memory_space_gpu);
}
#endif

TEST_F(ADIOS2_C_API, ADIOS2BPMemorySpaceShape)
{
#if ADIOS2_USE_MPI
    const char fname[] = "ADIOS2_C_API.ADIOS2BPMemorySpace_MPI.bp";
#else
    const char fname[] = "ADIOS2_C_API.ADIOS2BPMemorySpace.bp";
#endif
    // write
    {
        adios2_io *ioH = adios2_declare_io(adiosH, "CMemSpace");
        adios2_set_engine(ioH, "BPFile");

        size_t shape[2];
        shape[0] = 5;
        shape[1] = 2;

        size_t start[2];
        start[0] = 0;
        start[1] = 0;

        size_t count[2];
        count[0] = 5;
        count[1] = 2;

        int32_t data_I32[10] = {131072, 131073,  -131070, 131075,  -131068,
                                131077, -131066, 131079,  -131064, 131081};
        adios2_define_variable(ioH, "varI32", adios2_type_int32_t, 2, shape, start, count,
                               adios2_constant_dims_true);
        adios2_variable *varI32 = adios2_inquire_variable(ioH, "varI32");

        // test that the set memory space is Host
        adios2_memory_space mem;
        adios2_set_memory_space(varI32, adios2_memory_space_host);
        adios2_error e = adios2_get_memory_space(&mem, varI32);
        EXPECT_EQ(e, adios2_error_none);
        EXPECT_EQ(mem, adios2_memory_space_host);

        adios2_engine *engineH = adios2_open(ioH, fname, adios2_mode_write);
        adios2_put(engineH, varI32, data_I32, adios2_mode_deferred);
        adios2_close(engineH);
    }
    // read shape
    {
        adios2_io *ioH = adios2_declare_io(adiosH, "Reader");
        adios2_engine *engineH = adios2_open(ioH, fname, adios2_mode_readRandomAccess);
        size_t steps;
        adios2_steps(&steps, engineH);
        EXPECT_EQ(steps, 1);

        adios2_variable *varI32 = adios2_inquire_variable(ioH, "varI32");

        // test that the shape function returns the correct dimensions
        size_t cpu_shape[2];
        adios2_error e = adios2_variable_shape(cpu_shape, varI32);
        EXPECT_EQ(e, adios2_error_none);
        EXPECT_EQ(cpu_shape[0], 5);
        EXPECT_EQ(cpu_shape[1], 2);
        e = adios2_variable_shape_with_memory_space(cpu_shape, varI32, adios2_memory_space_host);
        EXPECT_EQ(e, adios2_error_none);
        EXPECT_EQ(cpu_shape[0], 5);
        EXPECT_EQ(cpu_shape[1], 2);
#ifdef ADIOS2_HAVE_GPU_SUPPORT
        size_t gpu_shape[2];
        e = adios2_variable_shape_with_memory_space(gpu_shape, varI32, adios2_memory_space_gpu);
        EXPECT_EQ(e, adios2_error_none);
        EXPECT_EQ(gpu_shape[0], 2);
        EXPECT_EQ(gpu_shape[1], 5);
#endif
        adios2_close(engineH);
    }
}

int main(int argc, char **argv)
{
    int result;
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();

    return result;
}

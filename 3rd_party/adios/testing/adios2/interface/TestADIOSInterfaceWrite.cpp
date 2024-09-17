#include <cstdint>

#include <iostream>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

class ADIOSInterfaceWriteTest : public ::testing::Test
{
public:
#if ADIOS2_USE_MPI
    ADIOSInterfaceWriteTest() : adios(MPI_COMM_SELF), io(adios.DeclareIO("TestIO")) {}
#else
    ADIOSInterfaceWriteTest() : adios(), io(adios.DeclareIO("TestIO")) {}
#endif

protected:
    // virtual void SetUp() { }

    // virtual void TearDown() { }

    adios2::ADIOS adios;
    adios2::IO io;
};

// 1x10
TEST_F(ADIOSInterfaceWriteTest, DefineVar_int8_t_1x10)
{
    int mpiRank = 0;
#if ADIOS2_USE_MPI
    int mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    std::string name = std::string("i8_") + std::to_string(mpiRank);

    // Define ADIOS variables for each type

    auto var_int8_t = io.DefineVariable<int8_t>(name, {}, {}, adios2::Dims{10});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(var_int8_t), adios2::Variable<int8_t>>();

    // Verify exceptions are thrown upon duplicate variable names
    EXPECT_THROW(auto foo = io.DefineVariable<int8_t>(name, {}, {}, adios2::Dims{10});
                 (void)foo, std::invalid_argument);

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(var_int8_t.Shape().size(), 0);
    EXPECT_EQ(var_int8_t.Start().size(), 0);
    EXPECT_EQ(var_int8_t.Count().size(), 1);
    EXPECT_EQ(var_int8_t.Count()[0], 10);
    EXPECT_EQ(var_int8_t.Name(), name);
    EXPECT_EQ(var_int8_t.Type(), "int8_t");
}

TEST_F(ADIOSInterfaceWriteTest, DefineVar_int16_t_1x10)
{
    int mpiRank = 0;
#if ADIOS2_USE_MPI
    int mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    std::string name = std::string("i16_") + std::to_string(mpiRank);

    // Define ADIOS variables for each type

    auto var_int16_t = io.DefineVariable<int16_t>(name, {}, {}, adios2::Dims{10});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(var_int16_t), adios2::Variable<int16_t>>();

    // Verify exceptions are thrown upon duplicate variable names
    EXPECT_THROW(auto foo = io.DefineVariable<int16_t>(name, {}, {}, adios2::Dims{10});
                 (void)foo, std::invalid_argument);

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(var_int16_t.Shape().size(), 0);
    EXPECT_EQ(var_int16_t.Start().size(), 0);
    EXPECT_EQ(var_int16_t.Count().size(), 1);
    EXPECT_EQ(var_int16_t.Count()[0], 10);
    EXPECT_EQ(var_int16_t.Name(), name);
    EXPECT_EQ(var_int16_t.Type(), "int16_t");
}

TEST_F(ADIOSInterfaceWriteTest, DefineVar_int32_t_1x10)
{
    int mpiRank = 0;
#if ADIOS2_USE_MPI
    int mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    std::string name = std::string("i32_") + std::to_string(mpiRank);

    // Define ADIOS variables for each type

    auto var_int32_t = io.DefineVariable<int32_t>(name, {}, {}, adios2::Dims{10});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(var_int32_t), adios2::Variable<int32_t>>();

    // Verify exceptions are thrown upon duplicate variable names
    EXPECT_THROW(auto foo = io.DefineVariable<int32_t>(name, {}, {}, adios2::Dims{10});
                 (void)foo, std::invalid_argument);

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(var_int32_t.Shape().size(), 0);
    EXPECT_EQ(var_int32_t.Start().size(), 0);
    EXPECT_EQ(var_int32_t.Count().size(), 1);
    EXPECT_EQ(var_int32_t.Count()[0], 10);
    EXPECT_EQ(var_int32_t.Name(), name);
    EXPECT_EQ(var_int32_t.Type(), "int32_t");
}

TEST_F(ADIOSInterfaceWriteTest, DefineVar_int64_t_1x10)
{
    int mpiRank = 0;

#if ADIOS2_USE_MPI
    int mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    std::string name = std::string("i64_") + std::to_string(mpiRank);

    // Define ADIOS variables for each type

    auto var_int64_t = io.DefineVariable<int64_t>(name, {}, {}, adios2::Dims{10});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(var_int64_t), adios2::Variable<int64_t>>();

    // Verify exceptions are thrown upon duplicate variable names
    EXPECT_THROW(auto foo = io.DefineVariable<int64_t>(name, {}, {}, adios2::Dims{10});
                 (void)foo, std::invalid_argument);

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(var_int64_t.Shape().size(), 0);
    EXPECT_EQ(var_int64_t.Start().size(), 0);
    EXPECT_EQ(var_int64_t.Count().size(), 1);
    EXPECT_EQ(var_int64_t.Count()[0], 10);
    EXPECT_EQ(var_int64_t.Name(), name);
    EXPECT_EQ(var_int64_t.Sizeof(), 8);
    EXPECT_EQ(var_int64_t.Type(), "int64_t");
}

TEST_F(ADIOSInterfaceWriteTest, DefineVar_uint8_t_1x10)
{
    int mpiRank = 0;

#if ADIOS2_USE_MPI
    int mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    std::string name = std::string("u8_") + std::to_string(mpiRank);

    // Define ADIOS variables for each type

    auto var_uint8_t = io.DefineVariable<uint8_t>(name, {}, {}, adios2::Dims{10});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(var_uint8_t), adios2::Variable<uint8_t>>();

    // Verify exceptions are thrown upon duplicate variable names
    EXPECT_THROW(auto foo = io.DefineVariable<uint8_t>(name, {}, {}, adios2::Dims{10});
                 (void)foo, std::invalid_argument);

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(var_uint8_t.Shape().size(), 0);
    EXPECT_EQ(var_uint8_t.Start().size(), 0);
    EXPECT_EQ(var_uint8_t.Count().size(), 1);
    EXPECT_EQ(var_uint8_t.Count()[0], 10);
    EXPECT_EQ(var_uint8_t.Name(), name);
    EXPECT_EQ(var_uint8_t.Type(), "uint8_t");
}

TEST_F(ADIOSInterfaceWriteTest, DefineVar_uint16_t_1x10)
{
    int mpiRank = 0;

#if ADIOS2_USE_MPI
    int mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    std::string name = std::string("u16_") + std::to_string(mpiRank);

    // Define ADIOS variables for each type

    auto var_uint16_t = io.DefineVariable<uint16_t>(name, {}, {}, adios2::Dims{10});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(var_uint16_t), adios2::Variable<uint16_t>>();

    // Verify exceptions are thrown upon duplicate variable names
    EXPECT_THROW(auto foo = io.DefineVariable<uint16_t>(name, {}, {}, adios2::Dims{10});
                 (void)foo, std::invalid_argument);

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(var_uint16_t.Shape().size(), 0);
    EXPECT_EQ(var_uint16_t.Start().size(), 0);
    EXPECT_EQ(var_uint16_t.Count().size(), 1);
    EXPECT_EQ(var_uint16_t.Count()[0], 10);
    EXPECT_EQ(var_uint16_t.Name(), name);
    EXPECT_EQ(var_uint16_t.Type(), "uint16_t");
}

TEST_F(ADIOSInterfaceWriteTest, DefineVar_uint32_t_1x10)
{
    int mpiRank = 0;

#if ADIOS2_USE_MPI
    int mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    std::string name = std::string("u32_") + std::to_string(mpiRank);

    // Define ADIOS variables for each type

    auto var_uint32_t = io.DefineVariable<uint32_t>(name, {}, {}, adios2::Dims{10});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(var_uint32_t), adios2::Variable<uint32_t>>();

    // Verify exceptions are thrown upon duplicate variable names
    EXPECT_THROW(auto foo = io.DefineVariable<uint32_t>(name, {}, {}, adios2::Dims{10});
                 (void)foo, std::invalid_argument);

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(var_uint32_t.Shape().size(), 0);
    EXPECT_EQ(var_uint32_t.Start().size(), 0);
    EXPECT_EQ(var_uint32_t.Count().size(), 1);
    EXPECT_EQ(var_uint32_t.Count()[0], 10);
    EXPECT_EQ(var_uint32_t.Name(), name);
    EXPECT_EQ(var_uint32_t.Type(), "uint32_t");
}

TEST_F(ADIOSInterfaceWriteTest, DefineVar_uint64_t_1x10)
{
    int mpiRank = 0;

#if ADIOS2_USE_MPI
    int mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    std::string name = std::string("u64_") + std::to_string(mpiRank);

    // Define ADIOS variables for each type

    auto var_uint64_t = io.DefineVariable<uint64_t>(name, {}, {}, adios2::Dims{10});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(var_uint64_t), adios2::Variable<uint64_t>>();

    // Verify exceptions are thrown upon duplicate variable names
    EXPECT_THROW(auto foo = io.DefineVariable<uint64_t>(name, {}, {}, adios2::Dims{10});
                 (void)foo, std::invalid_argument);

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(var_uint64_t.Shape().size(), 0);
    EXPECT_EQ(var_uint64_t.Start().size(), 0);
    EXPECT_EQ(var_uint64_t.Count().size(), 1);
    EXPECT_EQ(var_uint64_t.Count()[0], 10);
    EXPECT_EQ(var_uint64_t.Name(), name);
    EXPECT_EQ(var_uint64_t.Sizeof(), 8);
    EXPECT_EQ(var_uint64_t.Type(), "uint64_t");
}

// Rinse  and repeat for remaining types

// 2x5
TEST_F(ADIOSInterfaceWriteTest, DefineVar_int8_t_2x5)
{
    int mpiRank = 0;

#if ADIOS2_USE_MPI
    int mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    std::string name = std::string("i8_") + std::to_string(mpiRank);

    // Define ADIOS variables for each type

    auto var_int8_t = io.DefineVariable<int8_t>(name, {}, {}, adios2::Dims{2, 5});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(var_int8_t), adios2::Variable<int8_t>>();

    // Verify exceptions are thrown upon duplicate variable names
    EXPECT_THROW(auto foo = io.DefineVariable<int8_t>(name, {}, {}, adios2::Dims{2, 5});
                 (void)foo, std::invalid_argument);

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(var_int8_t.Shape().size(), 0);
    EXPECT_EQ(var_int8_t.Start().size(), 0);
    EXPECT_EQ(var_int8_t.Count().size(), 2);
    EXPECT_EQ(var_int8_t.Count()[0], 2);
    EXPECT_EQ(var_int8_t.Count()[1], 5);
    EXPECT_EQ(var_int8_t.Name(), name);
    EXPECT_EQ(var_int8_t.Type(), "int8_t");
}

TEST_F(ADIOSInterfaceWriteTest, DefineVar_int16_t_2x5)
{
    int mpiRank = 0;

#if ADIOS2_USE_MPI
    int mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    std::string name = std::string("i16_") + std::to_string(mpiRank);

    // Define ADIOS variables for each type

    auto var_int16_t = io.DefineVariable<int16_t>(name, {}, {}, adios2::Dims{2, 5});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(var_int16_t), adios2::Variable<int16_t>>();

    // Verify exceptions are thrown upon duplicate variable names
    EXPECT_THROW(auto foo = io.DefineVariable<int16_t>(name, {}, {}, adios2::Dims{2, 5});
                 (void)foo, std::invalid_argument);

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(var_int16_t.Shape().size(), 0);
    EXPECT_EQ(var_int16_t.Start().size(), 0);
    EXPECT_EQ(var_int16_t.Count().size(), 2);
    EXPECT_EQ(var_int16_t.Count()[0], 2);
    EXPECT_EQ(var_int16_t.Count()[1], 5);
    EXPECT_EQ(var_int16_t.Name(), name);
    EXPECT_EQ(var_int16_t.Type(), "int16_t");
}

TEST_F(ADIOSInterfaceWriteTest, DefineVar_int32_t_2x5)
{
    int mpiRank = 0;

#if ADIOS2_USE_MPI
    int mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    std::string name = std::string("i32_") + std::to_string(mpiRank);

    // Define ADIOS variables for each type

    auto var_int32_t = io.DefineVariable<int32_t>(name, {}, {}, adios2::Dims{2, 5});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(var_int32_t), adios2::Variable<int32_t>>();

    // Verify exceptions are thrown upon duplicate variable names
    EXPECT_THROW(auto foo = io.DefineVariable<int32_t>(name, {}, {}, adios2::Dims{2, 5});
                 (void)foo, std::invalid_argument);

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(var_int32_t.Shape().size(), 0);
    EXPECT_EQ(var_int32_t.Start().size(), 0);
    EXPECT_EQ(var_int32_t.Count().size(), 2);
    EXPECT_EQ(var_int32_t.Count()[0], 2);
    EXPECT_EQ(var_int32_t.Count()[1], 5);
    EXPECT_EQ(var_int32_t.Name(), name);
    EXPECT_EQ(var_int32_t.Type(), "int32_t");
}

TEST_F(ADIOSInterfaceWriteTest, DefineVar_int64_t_2x5)
{
    int mpiRank = 0;

#if ADIOS2_USE_MPI
    int mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    std::string name = std::string("i64_") + std::to_string(mpiRank);

    // Define ADIOS variables for each type

    auto var_int64_t = io.DefineVariable<int64_t>(name, {}, {}, adios2::Dims{2, 5});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(var_int64_t), adios2::Variable<int64_t>>();

    // Verify exceptions are thrown upon duplicate variable names
    EXPECT_THROW(auto foo = io.DefineVariable<int64_t>(name, {}, {}, adios2::Dims{2, 5});
                 (void)foo, std::invalid_argument);

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(var_int64_t.Shape().size(), 0);
    EXPECT_EQ(var_int64_t.Start().size(), 0);
    EXPECT_EQ(var_int64_t.Count().size(), 2);
    EXPECT_EQ(var_int64_t.Count()[0], 2);
    EXPECT_EQ(var_int64_t.Count()[1], 5);
    EXPECT_EQ(var_int64_t.Name(), name);
    EXPECT_EQ(var_int64_t.Sizeof(), 8);
    EXPECT_EQ(var_int64_t.Type(), "int64_t");
}

TEST_F(ADIOSInterfaceWriteTest, DefineVar_uint8_t_2x5)
{
    int mpiRank = 0;

#if ADIOS2_USE_MPI
    int mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    std::string name = std::string("u8_") + std::to_string(mpiRank);

    // Define ADIOS variables for each type

    auto var_uint8_t = io.DefineVariable<uint8_t>(name, {}, {}, adios2::Dims{2, 5});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(var_uint8_t), adios2::Variable<uint8_t>>();

    // Verify exceptions are thrown upon duplicate variable names
    EXPECT_THROW(auto foo = io.DefineVariable<uint8_t>(name, {}, {}, adios2::Dims{2, 5});
                 (void)foo, std::invalid_argument);

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(var_uint8_t.Shape().size(), 0);
    EXPECT_EQ(var_uint8_t.Start().size(), 0);
    EXPECT_EQ(var_uint8_t.Count().size(), 2);
    EXPECT_EQ(var_uint8_t.Count()[0], 2);
    EXPECT_EQ(var_uint8_t.Count()[1], 5);
    EXPECT_EQ(var_uint8_t.Name(), name);
    EXPECT_EQ(var_uint8_t.Type(), "uint8_t");
}

TEST_F(ADIOSInterfaceWriteTest, DefineVar_uint16_t_2x5)
{
    int mpiRank = 0;

#if ADIOS2_USE_MPI
    int mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    std::string name = std::string("u16_") + std::to_string(mpiRank);

    // Define ADIOS variables for each type

    auto var_uint16_t = io.DefineVariable<uint16_t>(name, {}, {}, adios2::Dims{2, 5});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(var_uint16_t), adios2::Variable<uint16_t>>();

    // Verify exceptions are thrown upon duplicate variable names
    EXPECT_THROW(auto foo = io.DefineVariable<uint16_t>(name, {}, {}, adios2::Dims{2, 5});
                 (void)foo, std::invalid_argument);

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(var_uint16_t.Shape().size(), 0);
    EXPECT_EQ(var_uint16_t.Start().size(), 0);
    EXPECT_EQ(var_uint16_t.Count().size(), 2);
    EXPECT_EQ(var_uint16_t.Count()[0], 2);
    EXPECT_EQ(var_uint16_t.Count()[1], 5);
    EXPECT_EQ(var_uint16_t.Name(), name);
    EXPECT_EQ(var_uint16_t.Type(), "uint16_t");
}

TEST_F(ADIOSInterfaceWriteTest, DefineVar_uint32_t_2x5)
{
    int mpiRank = 0;

#if ADIOS2_USE_MPI
    int mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    std::string name = std::string("u32_") + std::to_string(mpiRank);

    // Define ADIOS variables for each type

    auto var_uint32_t = io.DefineVariable<uint32_t>(name, {}, {}, adios2::Dims{2, 5});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(var_uint32_t), adios2::Variable<uint32_t>>();

    // Verify exceptions are thrown upon duplicate variable names
    EXPECT_THROW(auto foo = io.DefineVariable<uint32_t>(name, {}, {}, adios2::Dims{2, 5});
                 (void)foo, std::invalid_argument);

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(var_uint32_t.Shape().size(), 0);
    EXPECT_EQ(var_uint32_t.Start().size(), 0);
    EXPECT_EQ(var_uint32_t.Count().size(), 2);
    EXPECT_EQ(var_uint32_t.Count()[0], 2);
    EXPECT_EQ(var_uint32_t.Count()[1], 5);
    EXPECT_EQ(var_uint32_t.Name(), name);
    EXPECT_EQ(var_uint32_t.Type(), "uint32_t");
}

TEST_F(ADIOSInterfaceWriteTest, DefineVar_uint64_t_2x5)
{
    int mpiRank = 0;

#if ADIOS2_USE_MPI
    int mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    std::string name = std::string("u64_") + std::to_string(mpiRank);

    // Define ADIOS variables for each type

    auto var_uint64_t = io.DefineVariable<uint64_t>(name, {}, {}, adios2::Dims{2, 5});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(var_uint64_t), adios2::Variable<uint64_t>>();

    // Verify exceptions are thrown upon duplicate variable names
    EXPECT_THROW(auto foo = io.DefineVariable<uint64_t>(name, {}, {}, adios2::Dims{2, 5});
                 (void)foo, std::invalid_argument);

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(var_uint64_t.Shape().size(), 0);
    EXPECT_EQ(var_uint64_t.Start().size(), 0);
    EXPECT_EQ(var_uint64_t.Count().size(), 2);
    EXPECT_EQ(var_uint64_t.Count()[0], 2);
    EXPECT_EQ(var_uint64_t.Count()[1], 5);
    EXPECT_EQ(var_uint64_t.Name(), name);
    EXPECT_EQ(var_uint64_t.Sizeof(), 8);
    EXPECT_EQ(var_uint64_t.Type(), "uint64_t");
}

TEST_F(ADIOSInterfaceWriteTest, Exceptions)
{
    EXPECT_EQ(io.EngineType(), "File");

    adios2::Operator invalidOp = adios.InquireOperator("InvalidOp");
    EXPECT_FALSE(invalidOp);
    EXPECT_THROW(adios.AtIO("IOnull"), std::invalid_argument);
    EXPECT_THROW(adios.DefineOperator("WrongOp", "UnsupportedType"), std::invalid_argument);

#ifdef ADIOS2_HAVE_BZIP2
    EXPECT_NO_THROW(adios.DefineOperator("bzip2Op", "bzip2"));
    EXPECT_THROW(adios.DefineOperator("bzip2Op", "bzip2"), std::invalid_argument);
#else
    EXPECT_THROW(adios.DefineOperator("bzip2Op", "bzip2"), std::invalid_argument);
#endif

#ifdef ADIOS2_HAVE_ZFP
    EXPECT_NO_THROW(adios.DefineOperator("zfpOp", "zfp"));
    EXPECT_THROW(adios.DefineOperator("zfpOp", "zfp"), std::invalid_argument);
#else
    EXPECT_THROW(adios.DefineOperator("zfpOp", "zfp"), std::invalid_argument);
#endif

#ifdef ADIOS2_HAVE_SZ
    EXPECT_NO_THROW(adios.DefineOperator("szOp", "sz"));
    EXPECT_THROW(adios.DefineOperator("szOp", "sz"), std::invalid_argument);
#else
    EXPECT_THROW(adios.DefineOperator("szOp", "sz"), std::invalid_argument);
#endif
}

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

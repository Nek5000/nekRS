#include <cstdint>

#include <array>
#include <iostream>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

class ADIOSDefineVariableTest : public ::testing::Test
{
public:
#if ADIOS2_USE_MPI
    ADIOSDefineVariableTest() : adios(MPI_COMM_WORLD), io(adios.DeclareIO("TestIO"))
#else
    ADIOSDefineVariableTest() : adios(), io(adios.DeclareIO("TestIO"))
#endif
    {
    }

protected:
    adios2::ADIOS adios;
    adios2::IO io;
};

TEST_F(ADIOSDefineVariableTest, DefineGlobalValue)
{
    std::string name = std::string("globalValue");

    // Define ADIOS global value
    auto globalvalue = io.DefineVariable<int>(name);

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(globalvalue), adios2::Variable<int>>();

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(globalvalue.Shape().size(), 0);
    EXPECT_EQ(globalvalue.Start().size(), 0);
    EXPECT_EQ(globalvalue.Count().size(), 0);
    EXPECT_EQ(globalvalue.Name(), name);
    EXPECT_EQ(globalvalue.Type(), "int32_t");
}

TEST_F(ADIOSDefineVariableTest, DefineLocalValue)
{
    // Define ADIOS local value (a value changing across processes)
    auto localvalue = io.DefineVariable<int>("localvalue", {adios2::LocalValueDim});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(localvalue), adios2::Variable<int>>();

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(localvalue.Shape().size(), 1);
    EXPECT_EQ(localvalue.Shape()[0], adios2::LocalValueDim);
    EXPECT_EQ(localvalue.Start().size(), 1);
    EXPECT_EQ(localvalue.Count().size(), 1);
    EXPECT_EQ(localvalue.Name(), "localvalue");
    EXPECT_EQ(localvalue.Type(), "int32_t");
}

TEST_F(ADIOSDefineVariableTest, DefineGlobalArray)
{
    int mpiRank = 0, mpiSize = 1;
#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    const std::size_t Nx(10), Ny(20), Nz(30);

    adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), static_cast<size_t>(Ny * mpiSize),
                       static_cast<size_t>(Nz * mpiSize)};
    adios2::Dims start{static_cast<size_t>(Nx * mpiRank), static_cast<size_t>(Ny * mpiRank),
                       static_cast<size_t>(Nz * mpiRank)};
    adios2::Dims count{static_cast<size_t>(Nx), static_cast<size_t>(Ny), static_cast<size_t>(Nz)};
    // Define ADIOS global array
    auto globalarray = io.DefineVariable<int>("globalarray", shape, start, count);

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(globalarray), adios2::Variable<int>>();

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(globalarray.Shape().size(), 3);
    EXPECT_EQ(globalarray.Shape()[0], Nx * mpiSize);
    EXPECT_EQ(globalarray.Shape()[1], Ny * mpiSize);
    EXPECT_EQ(globalarray.Shape()[2], Nz * mpiSize);
    EXPECT_EQ(globalarray.Start().size(), 3);
    EXPECT_EQ(globalarray.Start()[0], Nx * mpiRank);
    EXPECT_EQ(globalarray.Start()[1], Ny * mpiRank);
    EXPECT_EQ(globalarray.Start()[2], Nz * mpiRank);
    EXPECT_EQ(globalarray.Count().size(), 3);
    EXPECT_EQ(globalarray.Count()[0], Nx);
    EXPECT_EQ(globalarray.Count()[1], Ny);
    EXPECT_EQ(globalarray.Count()[2], Nz);
    EXPECT_EQ(globalarray.Name(), "globalarray");
    EXPECT_EQ(globalarray.Type(), "int32_t");
}

TEST_F(ADIOSDefineVariableTest, DefineGlobalArrayWithSelections)
{
    // Define ADIOS global array with postponed size definition in SetSelection

    int mpiRank = 0, mpiSize = 1;
#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    const std::size_t Nx(10), Ny(20), Nz(30);

    adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), static_cast<size_t>(Ny * mpiSize),
                       static_cast<size_t>(Nz * mpiSize)};
    adios2::Dims start{static_cast<size_t>(Nx * mpiRank), static_cast<size_t>(Ny * mpiRank),
                       static_cast<size_t>(Nz * mpiRank)};
    adios2::Dims count{static_cast<size_t>(Nx), static_cast<size_t>(Ny), static_cast<size_t>(Nz)};
    // Define ADIOS global array
    auto globalarray = io.DefineVariable<int>("globalarray", shape);

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(globalarray), adios2::Variable<int>>();

    // Make a 3D selection to describe the local dimensions of the
    // variable we write and its offsets in the global spaces
    adios2::Box<adios2::Dims> sel(start, count);
    globalarray.SetSelection(sel);

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(globalarray.Shape().size(), 3);
    EXPECT_EQ(globalarray.Shape()[0], Nx * mpiSize);
    EXPECT_EQ(globalarray.Shape()[1], Ny * mpiSize);
    EXPECT_EQ(globalarray.Shape()[2], Nz * mpiSize);
    EXPECT_EQ(globalarray.Start().size(), 3);
    EXPECT_EQ(globalarray.Start()[0], Nx * mpiRank);
    EXPECT_EQ(globalarray.Start()[1], Ny * mpiRank);
    EXPECT_EQ(globalarray.Start()[2], Nz * mpiRank);
    EXPECT_EQ(globalarray.Count().size(), 3);
    EXPECT_EQ(globalarray.Count()[0], Nx);
    EXPECT_EQ(globalarray.Count()[1], Ny);
    EXPECT_EQ(globalarray.Count()[2], Nz);
    EXPECT_EQ(globalarray.Name(), "globalarray");
    EXPECT_EQ(globalarray.Type(), "int32_t");
}

TEST_F(ADIOSDefineVariableTest, DefineGlobalArrayConstantDims)
{
    // Define ADIOS global array with locked-down dimensions
    int mpiRank = 0, mpiSize = 1;
#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    const std::size_t Nx(10), Ny(20), Nz(30);

    adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), static_cast<size_t>(Ny * mpiSize),
                       static_cast<size_t>(Nz * mpiSize)};
    adios2::Dims start{static_cast<size_t>(Nx * mpiRank), static_cast<size_t>(Ny * mpiRank),
                       static_cast<size_t>(Nz * mpiRank)};
    adios2::Dims count{static_cast<size_t>(Nx), static_cast<size_t>(Ny), static_cast<size_t>(Nz)};
    // Define ADIOS global array
    auto globalarray = io.DefineVariable<int>("globalarray", shape, start, count, true);

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(globalarray), adios2::Variable<int>>();

    // Make a 3D selection to describe the local dimensions of the
    // variable we write and its offsets in the global spaces
    adios2::Box<adios2::Dims> sel(start, count);
    EXPECT_THROW(globalarray.SetSelection(sel), std::invalid_argument);

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(globalarray.Shape().size(), 3);
    EXPECT_EQ(globalarray.Shape()[0], Nx * mpiSize);
    EXPECT_EQ(globalarray.Shape()[1], Ny * mpiSize);
    EXPECT_EQ(globalarray.Shape()[2], Nz * mpiSize);
    EXPECT_EQ(globalarray.Start().size(), 3);
    EXPECT_EQ(globalarray.Start()[0], Nx * mpiRank);
    EXPECT_EQ(globalarray.Start()[1], Ny * mpiRank);
    EXPECT_EQ(globalarray.Start()[2], Nz * mpiRank);
    EXPECT_EQ(globalarray.Count().size(), 3);
    EXPECT_EQ(globalarray.Count()[0], Nx);
    EXPECT_EQ(globalarray.Count()[1], Ny);
    EXPECT_EQ(globalarray.Count()[2], Nz);
    EXPECT_EQ(globalarray.Name(), "globalarray");
    EXPECT_EQ(globalarray.Type(), "int32_t");
}

TEST_F(ADIOSDefineVariableTest, DefineGlobalArrayInvalidLocalValueDim)
{
    // Define ADIOS global array
    std::size_t n = 50;
    EXPECT_THROW(io.DefineVariable<int>("globalarray", {100, adios2::LocalValueDim, 30},
                                        {50, n / 2, 0}, {10, n / 2, 30}),
                 std::invalid_argument);
}

TEST_F(ADIOSDefineVariableTest, DefineLocalArray)
{
    // Define ADIOS local array (no global dimensions, no offsets)
    std::size_t n = 50;
    auto localarray = io.DefineVariable<int>("localarray", {}, {}, {10, n / 2, 30});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(localarray), adios2::Variable<int>>();

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(localarray.Shape().size(), 0);
    EXPECT_EQ(localarray.Start().size(), 0);
    EXPECT_EQ(localarray.Count().size(), 3);
    EXPECT_EQ(localarray.Count()[0], 10);
    EXPECT_EQ(localarray.Count()[1], n / 2);
    EXPECT_EQ(localarray.Count()[2], 30);
    EXPECT_EQ(localarray.Name(), "localarray");
    EXPECT_EQ(localarray.Type(), "int32_t");
    EXPECT_EQ(localarray.ShapeID(), adios2::ShapeID::LocalArray);
}

TEST_F(ADIOSDefineVariableTest, DefineLocalArrayWithSelection)
{
    int mpiRank = 0, mpiSize = 1;
#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    const std::size_t Nx(10), Ny(20), Nz(30);

    adios2::Dims shape{static_cast<size_t>(Nx * mpiSize), static_cast<size_t>(Ny * mpiSize),
                       static_cast<size_t>(Nz * mpiSize)};
    adios2::Dims start{static_cast<size_t>(Nx * mpiRank), static_cast<size_t>(Ny * mpiRank),
                       static_cast<size_t>(Nz * mpiRank)};
    adios2::Dims count{static_cast<size_t>(Nx), static_cast<size_t>(Ny), static_cast<size_t>(Nz)};
    // Define ADIOS global array
    auto localArray = io.DefineVariable<int>(
        "localArray", {}, {}, {adios2::UnknownDim, adios2::UnknownDim, adios2::UnknownDim});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(localArray), adios2::Variable<int>>();
    ASSERT_EQ(localArray.Shape().size(), 0);
    EXPECT_EQ(localArray.Start().size(), 0);
    EXPECT_EQ(localArray.Count().size(), 3);
    EXPECT_EQ(localArray.Count()[0], 0);
    EXPECT_EQ(localArray.Count()[1], 0);
    EXPECT_EQ(localArray.Count()[2], 0);
    EXPECT_EQ(localArray.Name(), "localArray");
    EXPECT_EQ(localArray.Type(), "int32_t");
    EXPECT_EQ(localArray.ShapeID(), adios2::ShapeID::LocalArray);

    // Make a 3D selection to describe the local dimensions of the
    // variable we write
    adios2::Box<adios2::Dims> sel({}, {Nx, Ny, Nz});
    localArray.SetSelection(sel);

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(localArray.Shape().size(), 0);
    EXPECT_EQ(localArray.Start().size(), 0);
    EXPECT_EQ(localArray.Count().size(), 3);
    EXPECT_EQ(localArray.Count()[0], Nx);
    EXPECT_EQ(localArray.Count()[1], Ny);
    EXPECT_EQ(localArray.Count()[2], Nz);
    EXPECT_EQ(localArray.Name(), "localArray");
    EXPECT_EQ(localArray.Type(), "int32_t");
    EXPECT_EQ(localArray.ShapeID(), adios2::ShapeID::LocalArray);

    // TODO: Put must not allow start in LocalArrays
    // adios2::Box<adios2::Dims> selbad(start, count);
    // EXPECT_THROW(localArray.SetSelection(selbad), std::invalid_argument);
}

TEST_F(ADIOSDefineVariableTest, DefineLocalArrayConstantDims)
{
#if ADIOS2_USE_MPI
    int mpiRank = 0, mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    const std::size_t Nx(10), Ny(20), Nz(30);

    adios2::Dims count{Nx, Ny, Nz};

    // Define ADIOS global array
    auto localArray = io.DefineVariable<int>("localArray", {}, {}, count, true);

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(localArray), adios2::Variable<int>>();

    adios2::Box<adios2::Dims> sel({}, count);
    EXPECT_THROW(localArray.SetSelection(sel), std::invalid_argument);

    ASSERT_EQ(localArray.Shape().size(), 0);
    EXPECT_EQ(localArray.Start().size(), 0);
    EXPECT_EQ(localArray.Count().size(), 3);
    EXPECT_EQ(localArray.Count()[0], Nx);
    EXPECT_EQ(localArray.Count()[1], Ny);
    EXPECT_EQ(localArray.Count()[2], Nz);
    EXPECT_EQ(localArray.Name(), "localArray");
    EXPECT_EQ(localArray.Type(), "int32_t");
    EXPECT_EQ(localArray.ShapeID(), adios2::ShapeID::LocalArray);
}

TEST_F(ADIOSDefineVariableTest, DefineLocalArrayInvalidOffsets)
{
    // Define ADIOS local array but try to add offsets
    std::size_t n = 50;

    EXPECT_THROW(io.DefineVariable<int>("localarray", {}, {50, n / 2, 0}, {10, n / 2, 30}),
                 std::invalid_argument);
}

TEST_F(ADIOSDefineVariableTest, DefineJoinedArrayFirstDim)
{
    // Define ADIOS joined array
    std::size_t n = 50;
    auto joinedarray =
        io.DefineVariable<int>("joinedarray", {adios2::JoinedDim, n, 30}, {}, {10, n, 30});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(joinedarray), adios2::Variable<int>>();

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(joinedarray.Shape().size(), 3);
    EXPECT_EQ(joinedarray.Shape()[0], adios2::JoinedDim);
    EXPECT_EQ(joinedarray.Shape()[1], n);
    EXPECT_EQ(joinedarray.Shape()[2], 30);
    EXPECT_EQ(joinedarray.Start().size(), 0);
    EXPECT_EQ(joinedarray.Count().size(), 3);
    EXPECT_EQ(joinedarray.Count()[0], 10);
    EXPECT_EQ(joinedarray.Count()[1], n);
    EXPECT_EQ(joinedarray.Count()[2], 30);
    EXPECT_EQ(joinedarray.Name(), "joinedarray");
    EXPECT_EQ(joinedarray.Type(), "int32_t");
}

TEST_F(ADIOSDefineVariableTest, DefineJoinedArraySecondDim)
{
    // Define ADIOS joined array
    std::size_t n = 50;
    auto joinedarray =
        io.DefineVariable<int>("joinedarray", {n, adios2::JoinedDim, 30}, {0, 0, 0}, {n, 10, 30});

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(joinedarray), adios2::Variable<int>>();

    // Verify the dimensions, name, and type are correct
    ASSERT_EQ(joinedarray.Shape().size(), 3);
    EXPECT_EQ(joinedarray.Shape()[0], n);
    EXPECT_EQ(joinedarray.Shape()[1], adios2::JoinedDim);
    EXPECT_EQ(joinedarray.Shape()[2], 30);
    EXPECT_EQ(joinedarray.Start().size(), 3);
    EXPECT_EQ(joinedarray.Start()[0], 0);
    EXPECT_EQ(joinedarray.Start()[1], 0);
    EXPECT_EQ(joinedarray.Start()[2], 0);
    EXPECT_EQ(joinedarray.Count().size(), 3);
    EXPECT_EQ(joinedarray.Count()[0], n);
    EXPECT_EQ(joinedarray.Count()[1], 10);
    EXPECT_EQ(joinedarray.Count()[2], 30);
    EXPECT_EQ(joinedarray.Name(), "joinedarray");
    EXPECT_EQ(joinedarray.Type(), "int32_t");
}

TEST_F(ADIOSDefineVariableTest, DefineJoinedArrayTooManyJoinedDims)
{
    // Define ADIOS joined array
    std::size_t n = 50;

    EXPECT_THROW(io.DefineVariable<int>("joinedarray", {n, adios2::JoinedDim, adios2::JoinedDim},
                                        {}, {n, 50, 30}),
                 std::invalid_argument);
}

TEST_F(ADIOSDefineVariableTest, DefineJoinedArrayInvalidStart)
{
    // Define ADIOS joined array
    std::size_t n = 10;
    std::size_t WrongValue = 1;
    // Start must be empty or full zero array
    EXPECT_THROW(
        io.DefineVariable<int>("joinedarray", {adios2::JoinedDim, 50}, {0, WrongValue}, {n, 50}),
        std::invalid_argument);
}

TEST_F(ADIOSDefineVariableTest, DefineString)
{
    // Define ADIOS local array but try to add offsets
    const std::size_t n = 50;
    EXPECT_THROW(
        io.DefineVariable<std::string>("invalidString1", {}, {50, n / 2, 0}, {10, n / 2, 30}),
        std::invalid_argument);
    EXPECT_THROW(io.DefineVariable<std::string>("invalidString2", {}, {}, {1}),
                 std::invalid_argument);

    EXPECT_NO_THROW(io.DefineVariable<std::string>("validString1"));
    EXPECT_NO_THROW(io.DefineVariable<std::string>("validString2", {}, {}, {}));
}

TEST_F(ADIOSDefineVariableTest, DefineAndRemove)
{
    auto lf_CheckRemove = [&](const std::string variableName) {
        const bool isRemoved = io.RemoveVariable(variableName);
        EXPECT_EQ(isRemoved, true);
    };

    const adios2::Dims shape = {10};
    const adios2::Dims start = {0};
    const adios2::Dims count = {10};

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
    io.DefineVariable<std::complex<float>>("c32", shape, start, count);
    io.DefineVariable<std::complex<double>>("c64", shape, start, count);
    io.DefineVariable<char>("char", shape, start, count);
    io.DefineVariable<long>("l", shape, start, count);
    io.DefineVariable<unsigned long>("ul", shape, start, count);

    lf_CheckRemove("iString");
    lf_CheckRemove("i8");
    lf_CheckRemove("i16");
    lf_CheckRemove("i32");
    lf_CheckRemove("i64");

    lf_CheckRemove("u8");
    lf_CheckRemove("u16");
    lf_CheckRemove("u32");
    lf_CheckRemove("u64");

    lf_CheckRemove("r32");
    lf_CheckRemove("r64");

    lf_CheckRemove("c32");
    lf_CheckRemove("c64");

    lf_CheckRemove("char");
    lf_CheckRemove("l");
    lf_CheckRemove("ul");

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
    auto var_c32 = io.InquireVariable<std::complex<float>>("c32");
    auto var_c64 = io.InquireVariable<std::complex<double>>("c64");
    auto var_char = io.InquireVariable<char>("char");
    auto var_l = io.InquireVariable<char>("l");
    auto var_ul = io.InquireVariable<char>("ul");

    EXPECT_FALSE(var_iString);
    EXPECT_FALSE(var_i8);
    EXPECT_FALSE(var_i16);
    EXPECT_FALSE(var_i32);
    EXPECT_FALSE(var_i64);

    EXPECT_FALSE(var_u8);
    EXPECT_FALSE(var_u16);
    EXPECT_FALSE(var_u32);
    EXPECT_FALSE(var_u64);

    EXPECT_FALSE(var_r32);
    EXPECT_FALSE(var_r64);

    EXPECT_FALSE(var_c32);
    EXPECT_FALSE(var_c64);

    EXPECT_FALSE(var_char);
    EXPECT_FALSE(var_l);
    EXPECT_FALSE(var_ul);
}

TEST_F(ADIOSDefineVariableTest, DefineRemoveDefine)
{
    auto lf_CheckRemove = [&](const std::string variableName) {
        const bool isRemoved = io.RemoveVariable(variableName);
        EXPECT_EQ(isRemoved, true);
    };

    const adios2::Dims shape = {10};
    const adios2::Dims start = {0};
    const adios2::Dims count = {10};

    io.DefineVariable<int8_t>("i8_0", shape, start, count);
    io.DefineVariable<int8_t>("i8_1", shape, start, count);

    lf_CheckRemove("i8_0");

    std::array<adios2::Variable<int8_t>, 2> i8_vars;

    i8_vars[0] = io.InquireVariable<int8_t>("i8_0");
    EXPECT_FALSE(i8_vars[0]);

    i8_vars[1] = io.InquireVariable<int8_t>("i8_1");
    EXPECT_TRUE(i8_vars[1]);
    EXPECT_EQ(i8_vars[1].Name(), "i8_1");

    io.DefineVariable<int8_t>("i8_0", shape, start, count);

    // check again after defining variable
    i8_vars[1] = io.InquireVariable<int8_t>("i8_1");
    EXPECT_TRUE(i8_vars[1]);
    EXPECT_EQ(i8_vars[1].Name(), "i8_1");

    i8_vars[0] = io.InquireVariable<int8_t>("i8_0");
    EXPECT_TRUE(i8_vars[0]);
    EXPECT_EQ(i8_vars[0].Name(), "i8_0");

    auto i8_var2 = io.DefineVariable<int8_t>("i8_2", shape, start, count);
    EXPECT_TRUE(i8_var2);
    EXPECT_EQ(i8_var2.Name(), "i8_2");
}

TEST_F(ADIOSDefineVariableTest, DefineAndRemoveAll)
{
    const adios2::Dims shape = {10};
    const adios2::Dims start = {0};
    const adios2::Dims count = {10};

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
    io.DefineVariable<std::complex<float>>("c32", shape, start, count);
    io.DefineVariable<std::complex<double>>("c64", shape, start, count);
    io.DefineVariable<char>("char", shape, start, count);
    io.DefineVariable<long>("l", shape, start, count);
    io.DefineVariable<unsigned long>("ul", shape, start, count);

    io.RemoveAllVariables();

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
    auto var_c32 = io.InquireVariable<std::complex<float>>("c32");
    auto var_c64 = io.InquireVariable<std::complex<double>>("c64");
    auto var_char = io.InquireVariable<char>("char");
    auto var_l = io.InquireVariable<char>("l");
    auto var_ul = io.InquireVariable<char>("ul");

    EXPECT_FALSE(var_iString);
    EXPECT_FALSE(var_i8);
    EXPECT_FALSE(var_i16);
    EXPECT_FALSE(var_i32);
    EXPECT_FALSE(var_i64);

    EXPECT_FALSE(var_u8);
    EXPECT_FALSE(var_u16);
    EXPECT_FALSE(var_u32);
    EXPECT_FALSE(var_u64);

    EXPECT_FALSE(var_r32);
    EXPECT_FALSE(var_r64);

    EXPECT_FALSE(var_c32);
    EXPECT_FALSE(var_c64);

    EXPECT_FALSE(var_char);
    EXPECT_FALSE(var_l);
    EXPECT_FALSE(var_ul);
}

TEST_F(ADIOSDefineVariableTest, DefineCheckType)
{
    const adios2::Dims shape = {10};
    const adios2::Dims start = {0};
    const adios2::Dims count = {10};

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
    io.DefineVariable<std::complex<float>>("c32", shape, start, count);
    io.DefineVariable<std::complex<double>>("c64", shape, start, count);
    io.DefineVariable<char>("char", shape, start, count);
    io.DefineVariable<long>("l", shape, start, count);
    io.DefineVariable<unsigned long>("ul", shape, start, count);

    EXPECT_EQ(io.VariableType("iString"), adios2::GetType<std::string>());
    EXPECT_EQ(io.VariableType("i8"), adios2::GetType<int8_t>());
    EXPECT_EQ(io.VariableType("i16"), adios2::GetType<int16_t>());
    EXPECT_EQ(io.VariableType("i32"), adios2::GetType<int32_t>());
    EXPECT_EQ(io.VariableType("i64"), adios2::GetType<int64_t>());
    EXPECT_EQ(io.VariableType("u8"), adios2::GetType<uint8_t>());
    EXPECT_EQ(io.VariableType("u16"), adios2::GetType<uint16_t>());
    EXPECT_EQ(io.VariableType("u32"), adios2::GetType<uint32_t>());
    EXPECT_EQ(io.VariableType("u64"), adios2::GetType<uint64_t>());
    EXPECT_EQ(io.VariableType("r32"), adios2::GetType<float>());
    EXPECT_EQ(io.VariableType("r64"), adios2::GetType<double>());
    EXPECT_EQ(io.VariableType("c32"), adios2::GetType<std::complex<float>>());
    EXPECT_EQ(io.VariableType("c64"), adios2::GetType<std::complex<double>>());
    EXPECT_EQ(io.VariableType("char"), adios2::GetType<char>());
    EXPECT_EQ(io.VariableType("l"), adios2::GetType<long>());
    EXPECT_EQ(io.VariableType("ul"), adios2::GetType<unsigned long>());
}

TEST_F(ADIOSDefineVariableTest, DefineStructVariable)
{
    const adios2::Dims shape = {10};
    const adios2::Dims start = {0};
    const adios2::Dims count = {10};

    typedef struct def1
    {
        int8_t a;
        int32_t b[5];
    } def1;
    auto struct1 = io.DefineStruct("def1", sizeof(def1));
    struct1.AddField("a", offsetof(def1, a), adios2::DataType::Int8);
    struct1.AddField("b", offsetof(def1, b), adios2::DataType::Int32, 5);
    struct1.Freeze();
    EXPECT_THROW(struct1.AddField("c", 0, adios2::DataType::Int32), std::runtime_error);

    typedef struct def2
    {
        int8_t a;
        int32_t b[5];
        int32_t c;
    } def2;
    auto struct2 = io.DefineStruct("def2", sizeof(def2));
    struct2.AddField("a", offsetof(def2, a), adios2::DataType::Int8);
    struct2.AddField("b", offsetof(def2, b), adios2::DataType::Int32, 5);
    struct2.AddField("c", 24, adios2::DataType::Int32);
    EXPECT_THROW(struct2.AddField("c", 27, adios2::DataType::Int32), std::runtime_error);

    auto structVar = io.DefineStructVariable("particle", struct1, shape, start, count);

    EXPECT_EQ(structVar.Shape().size(), 1);
    EXPECT_EQ(structVar.Start().size(), 1);
    EXPECT_EQ(structVar.Count().size(), 1);
    EXPECT_EQ(structVar.Name(), "particle");
    EXPECT_EQ(structVar.Type(), "struct");
    EXPECT_EQ(structVar.Sizeof(), 24);
    EXPECT_EQ(structVar.StructFields(), 2);
    EXPECT_EQ(structVar.StructFieldName(0), "a");
    EXPECT_EQ(structVar.StructFieldName(1), "b");
    EXPECT_EQ(structVar.StructFieldOffset(0), 0);
    EXPECT_EQ(structVar.StructFieldOffset(1), offsetof(def1, b));
    EXPECT_EQ(structVar.StructFieldType(0), adios2::DataType::Int8);
    EXPECT_EQ(structVar.StructFieldType(1), adios2::DataType::Int32);
    EXPECT_EQ(structVar.StructFieldElementCount(0), 1);
    EXPECT_EQ(structVar.StructFieldElementCount(1), 5);

    EXPECT_THROW(structVar.StructFieldName(2), std::invalid_argument);
    EXPECT_THROW(structVar.StructFieldOffset(2), std::invalid_argument);
    EXPECT_THROW(structVar.StructFieldType(2), std::invalid_argument);
    EXPECT_THROW(structVar.StructFieldElementCount(2), std::invalid_argument);

    auto inquire1 = io.InquireVariable("particle");
    EXPECT_TRUE(inquire1);

    auto inquire2 = io.InquireStructVariable("particle");
    EXPECT_TRUE(inquire2);

    auto inquire3 = io.InquireStructVariable("particle", struct1);
    EXPECT_TRUE(inquire3);

    auto inquire4 = io.InquireStructVariable("particle", struct2);
    EXPECT_FALSE(inquire4);
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

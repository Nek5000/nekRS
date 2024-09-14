#include <cstdint>

#include <array>
#include <iostream>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

#include "../engine/SmallTestData.h"

class ADIOSDefineAttributeTest : public ::testing::Test
{
public:
#if ADIOS2_USE_MPI
    ADIOSDefineAttributeTest() : adios(MPI_COMM_WORLD), io(adios.DeclareIO("TestIO"))
#else
    ADIOSDefineAttributeTest() : adios(), io(adios.DeclareIO("TestIO"))
#endif
    {
    }

    SmallTestData m_TestData;

protected:
    adios2::ADIOS adios;
    adios2::IO io;
};

TEST_F(ADIOSDefineAttributeTest, DefineAttributeNameException)
{
    int mpiRank = 0;

#if ADIOS2_USE_MPI
    int mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    std::string name = std::string("attributeString") + std::to_string(mpiRank);

    // Attribute should be unique per process
    io.DefineAttribute<std::string>(name, "-1");
    auto availableAttributes = io.AvailableAttributes();
    EXPECT_EQ(availableAttributes.size(), 1);

    // Redefinition is not allowed (non-modifiable attribute)
    EXPECT_THROW(io.DefineAttribute<std::string>(name, "0"), std::invalid_argument);

    auto attributeString1 = io.InquireAttribute<std::string>("NonExistingAttribute");
    EXPECT_FALSE(attributeString1);

    auto attributeString2 = io.InquireAttribute<std::string>(name);
    EXPECT_TRUE(attributeString2);

    /* Modifiable attribute can change its value(s) ... */
    io.DefineAttribute<std::string>("modifiable", "initial", "", "", true);
    io.DefineAttribute<std::string>("modifiable", "modified", "", "", true);

    auto attributeString3 = io.InquireAttribute<std::string>("modifiable");
    EXPECT_TRUE(attributeString3);
    auto attributeString3Value = attributeString3.Data();
    ASSERT_EQ(attributeString3Value.size() == 1, true);
    EXPECT_EQ(attributeString3Value[0], "modified");

    /* ... but not its type */
    EXPECT_THROW(io.DefineAttribute<double>("modifiable", 1.0, "", "", true),
                 std::invalid_argument);
}

TEST_F(ADIOSDefineAttributeTest, DefineAttributeTypeByValue)
{
    int mpiRank = 0;
    int mpiSize = 1;
#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // Define unique data for each process
    SmallTestData currentTestData = generateNewSmallTestData(m_TestData, 0, mpiRank, mpiSize);

    const std::string mpiRankString(std::to_string(mpiRank));
    const std::string s1_Single("s1_Single_" + mpiRankString);
    const std::string i8_Single("i8_Single_" + mpiRankString);
    const std::string i16_Single("i16_Single_" + mpiRankString);
    const std::string i32_Single("i32_Single_" + mpiRankString);
    const std::string i64_Single("i64_Single_" + mpiRankString);
    const std::string u8_Single("u8_Single_" + mpiRankString);
    const std::string u16_Single("u16_Single_" + mpiRankString);
    const std::string u32_Single("u32_Single_" + mpiRankString);
    const std::string u64_Single("u64_Single_" + mpiRankString);
    const std::string float_Single("float_Single_" + mpiRankString);
    const std::string double_Single("double_Single_" + mpiRankString);

    // Define ADIOS global value
    auto attributeS1 = io.DefineAttribute<std::string>(s1_Single, currentTestData.S1);
    auto attributeI8 = io.DefineAttribute<int8_t>(i8_Single, currentTestData.I8.front());
    auto attributeI16 = io.DefineAttribute<int16_t>(i16_Single, currentTestData.I16.front());
    auto attributeI32 = io.DefineAttribute<int32_t>(i32_Single, currentTestData.I32.front());
    auto attributeI64 = io.DefineAttribute<int64_t>(i64_Single, currentTestData.I64.front());

    auto attributeU8 = io.DefineAttribute<uint8_t>(u8_Single, currentTestData.U8.front());
    auto attributeU16 = io.DefineAttribute<uint16_t>(u16_Single, currentTestData.U16.front());
    auto attributeU32 = io.DefineAttribute<uint32_t>(u32_Single, currentTestData.U32.front());
    auto attributeU64 = io.DefineAttribute<uint64_t>(u64_Single, currentTestData.U64.front());

    auto attributeFloat = io.DefineAttribute<float>(float_Single, currentTestData.R32.front());
    auto attributeDouble = io.DefineAttribute<double>(double_Single, currentTestData.R64.front());

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(attributeS1), adios2::Attribute<std::string>>();
    ::testing::StaticAssertTypeEq<decltype(attributeI8), adios2::Attribute<int8_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeI16), adios2::Attribute<int16_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeI32), adios2::Attribute<int32_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeI64), adios2::Attribute<int64_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeU8), adios2::Attribute<uint8_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeU16), adios2::Attribute<uint16_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeU32), adios2::Attribute<uint32_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeU64), adios2::Attribute<uint64_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeFloat), adios2::Attribute<float>>();
    ::testing::StaticAssertTypeEq<decltype(attributeDouble), adios2::Attribute<double>>();

    // Verify the members are correct
    std::cout << "STRING name: " << attributeS1.Name() << "\n";

    ASSERT_EQ(attributeS1.Data().size() == 1, true);
    ASSERT_EQ(attributeS1.Data().empty(), false);
    EXPECT_EQ(attributeS1.Name(), s1_Single);
    EXPECT_EQ(attributeS1.Data()[0], currentTestData.S1);
    EXPECT_EQ(attributeS1.Data().size(), 1);
    EXPECT_EQ(attributeS1.Type(), "string");
    EXPECT_TRUE(attributeS1.IsValue());

    ASSERT_EQ(attributeI8.Data().size() == 1, true);
    ASSERT_EQ(attributeI8.Data().empty(), false);
    EXPECT_EQ(attributeI8.Name(), i8_Single);
    EXPECT_EQ(attributeI8.Data()[0], currentTestData.I8.front());
    EXPECT_EQ(attributeI8.Data().size(), 1);
    EXPECT_EQ(attributeI8.Type(), "int8_t");
    EXPECT_TRUE(attributeI8.IsValue());

    ASSERT_EQ(attributeI16.Data().size() == 1, true);
    ASSERT_EQ(attributeI16.Data().empty(), false);
    EXPECT_EQ(attributeI16.Name(), i16_Single);
    EXPECT_EQ(attributeI16.Data()[0], currentTestData.I16.front());
    EXPECT_EQ(attributeI16.Data().size(), 1);
    EXPECT_EQ(attributeI16.Type(), "int16_t");
    EXPECT_TRUE(attributeI16.IsValue());

    ASSERT_EQ(attributeI32.Data().size() == 1, true);
    ASSERT_EQ(attributeI32.Data().empty(), false);
    EXPECT_EQ(attributeI32.Name(), i32_Single);
    EXPECT_EQ(attributeI32.Data()[0], currentTestData.I32.front());
    EXPECT_EQ(attributeI32.Data().size(), 1);
    EXPECT_EQ(attributeI32.Type(), "int32_t");
    EXPECT_TRUE(attributeI32.IsValue());

    ASSERT_EQ(attributeI64.Data().size() == 1, true);
    ASSERT_EQ(attributeI64.Data().empty(), false);
    EXPECT_EQ(attributeI64.Name(), i64_Single);
    EXPECT_EQ(attributeI64.Data()[0], currentTestData.I64.front());
    EXPECT_EQ(attributeI64.Data().size(), 1);
    EXPECT_EQ(attributeI64.Type(), "int64_t");
    EXPECT_TRUE(attributeI64.IsValue());
    EXPECT_EQ(sizeof(attributeI64.Data()[0]), 8);

    ASSERT_EQ(attributeU8.Data().size() == 1, true);
    ASSERT_EQ(attributeU8.Data().empty(), false);
    EXPECT_EQ(attributeU8.Name(), u8_Single);
    EXPECT_EQ(attributeU8.Data()[0], currentTestData.U8.front());
    EXPECT_EQ(attributeU8.Data().size(), 1);
    EXPECT_EQ(attributeU8.Type(), "uint8_t");
    EXPECT_TRUE(attributeU8.IsValue());

    ASSERT_EQ(attributeU16.Data().size() == 1, true);
    ASSERT_EQ(attributeU16.Data().empty(), false);
    EXPECT_EQ(attributeU16.Name(), u16_Single);
    EXPECT_EQ(attributeU16.Data()[0], currentTestData.U16.front());
    EXPECT_EQ(attributeU16.Data().size(), 1);
    EXPECT_EQ(attributeU16.Type(), "uint16_t");
    EXPECT_TRUE(attributeU16.IsValue());

    ASSERT_EQ(attributeU32.Data().size() == 1, true);
    ASSERT_EQ(attributeU32.Data().empty(), false);
    EXPECT_EQ(attributeU32.Name(), u32_Single);
    EXPECT_EQ(attributeU32.Data()[0], currentTestData.U32.front());
    EXPECT_EQ(attributeU32.Data().size(), 1);
    EXPECT_EQ(attributeU32.Type(), "uint32_t");
    EXPECT_TRUE(attributeU32.IsValue());

    ASSERT_EQ(attributeU64.Data().size() == 1, true);
    ASSERT_EQ(attributeU64.Data().empty(), false);
    EXPECT_EQ(attributeU64.Name(), u64_Single);
    EXPECT_EQ(attributeU64.Data()[0], currentTestData.U64.front());
    EXPECT_EQ(attributeU64.Data().size(), 1);
    EXPECT_EQ(attributeU64.Type(), "uint64_t");
    EXPECT_TRUE(attributeU64.IsValue());
    EXPECT_EQ(sizeof(attributeU64.Data()[0]), 8);

    ASSERT_EQ(attributeFloat.Data().size() == 1, true);
    ASSERT_EQ(attributeFloat.Data().empty(), false);
    EXPECT_EQ(attributeFloat.Name(), float_Single);
    EXPECT_EQ(attributeFloat.Data()[0], currentTestData.R32.front());
    EXPECT_EQ(attributeFloat.Data().size(), 1);
    EXPECT_EQ(attributeFloat.Type(), "float");
    EXPECT_TRUE(attributeFloat.IsValue());

    ASSERT_EQ(attributeDouble.Data().size() == 1, true);
    ASSERT_EQ(attributeDouble.Data().empty(), false);
    EXPECT_EQ(attributeDouble.Name(), double_Single);
    EXPECT_EQ(attributeDouble.Data()[0], currentTestData.R64.front());
    EXPECT_EQ(attributeDouble.Data().size(), 1);
    EXPECT_EQ(attributeDouble.Type(), "double");
    EXPECT_TRUE(attributeDouble.IsValue());
}

TEST_F(ADIOSDefineAttributeTest, DefineAttributeTypeByReference)
{
    int mpiRank = 0, mpiSize = 1;
    size_t numberOfElements = 10;
#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // Define unique data for each process
    SmallTestData currentTestData = generateNewSmallTestData(m_TestData, 0, mpiRank, mpiSize);

    std::string mpiRankString = std::to_string(mpiRank);
    std::string s3_Single = std::string("s3_Single_") + mpiRankString;
    std::string i8_Single = std::string("i8_Single_") + mpiRankString;
    std::string i16_Single = std::string("i16_Single_") + mpiRankString;
    std::string i32_Single = std::string("i32_Single_") + mpiRankString;
    std::string i64_Single = std::string("i64_Single_") + mpiRankString;
    std::string u8_Single = std::string("u8_Single_") + mpiRankString;
    std::string u16_Single = std::string("u16_Single_") + mpiRankString;
    std::string u32_Single = std::string("u32_Single_") + mpiRankString;
    std::string u64_Single = std::string("u64_Single_") + mpiRankString;
    std::string float_Single = std::string("float_Single_") + mpiRankString;
    std::string double_Single = std::string("double_Single_") + mpiRankString;

    // Define ADIOS global value
    auto attributeS3 = io.DefineAttribute<std::string>(s3_Single, currentTestData.S3.data(), 3);
    auto attributeI8 =
        io.DefineAttribute<int8_t>(i8_Single, currentTestData.I8.data(), numberOfElements);
    auto attributeI16 =
        io.DefineAttribute<int16_t>(i16_Single, currentTestData.I16.data(), numberOfElements);
    auto attributeI32 =
        io.DefineAttribute<int32_t>(i32_Single, currentTestData.I32.data(), numberOfElements);
    auto attributeI64 =
        io.DefineAttribute<int64_t>(i64_Single, currentTestData.I64.data(), numberOfElements);

    auto attributeU8 =
        io.DefineAttribute<uint8_t>(u8_Single, currentTestData.U8.data(), numberOfElements);
    auto attributeU16 =
        io.DefineAttribute<uint16_t>(u16_Single, currentTestData.U16.data(), numberOfElements);
    auto attributeU32 =
        io.DefineAttribute<uint32_t>(u32_Single, currentTestData.U32.data(), numberOfElements);
    auto attributeU64 =
        io.DefineAttribute<uint64_t>(u64_Single, currentTestData.U64.data(), numberOfElements);

    auto attributeFloat =
        io.DefineAttribute<float>(float_Single, currentTestData.R32.data(), numberOfElements);
    auto attributeDouble =
        io.DefineAttribute<double>(double_Single, currentTestData.R64.data(), numberOfElements);

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(attributeS3), adios2::Attribute<std::string>>();
    ::testing::StaticAssertTypeEq<decltype(attributeI8), adios2::Attribute<int8_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeI16), adios2::Attribute<int16_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeI32), adios2::Attribute<int32_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeI64), adios2::Attribute<int64_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeU8), adios2::Attribute<uint8_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeU16), adios2::Attribute<uint16_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeU32), adios2::Attribute<uint32_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeU64), adios2::Attribute<uint64_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeFloat), adios2::Attribute<float>>();
    ::testing::StaticAssertTypeEq<decltype(attributeDouble), adios2::Attribute<double>>();

    // Verify the members are correct
    ASSERT_EQ(attributeS3.Data().size() == 1, false);
    ASSERT_EQ(attributeS3.Data().empty(), false);
    EXPECT_EQ(attributeS3.Name(), s3_Single);
    EXPECT_EQ(attributeS3.Data().size(), 3);
    EXPECT_EQ(attributeS3.Type(), "string");

    ASSERT_EQ(attributeI8.Data().size() == 1, false);
    ASSERT_EQ(attributeI8.Data().empty(), false);
    EXPECT_EQ(attributeI8.Name(), i8_Single);
    EXPECT_EQ(attributeI8.Data().size(), numberOfElements);
    EXPECT_EQ(attributeI8.Type(), "int8_t");

    ASSERT_EQ(attributeI16.Data().size() == 1, false);
    ASSERT_EQ(attributeI16.Data().empty(), false);
    EXPECT_EQ(attributeI16.Name(), i16_Single);
    EXPECT_EQ(attributeI16.Data().size(), numberOfElements);
    EXPECT_EQ(attributeI16.Type(), "int16_t");

    ASSERT_EQ(attributeI32.Data().size() == 1, false);
    ASSERT_EQ(attributeI32.Data().empty(), false);
    EXPECT_EQ(attributeI32.Name(), i32_Single);
    EXPECT_EQ(attributeI32.Data().size(), numberOfElements);
    EXPECT_EQ(attributeI32.Type(), "int32_t");

    ASSERT_EQ(attributeI64.Data().size() == 1, false);
    ASSERT_EQ(attributeI64.Data().empty(), false);
    EXPECT_EQ(attributeI64.Name(), i64_Single);
    EXPECT_EQ(attributeI64.Data().size(), numberOfElements);
    EXPECT_EQ(attributeI64.Type(), "int64_t");
    EXPECT_EQ(sizeof(attributeI64.Data()[0]), 8);

    ASSERT_EQ(attributeU8.Data().size() == 1, false);
    ASSERT_EQ(attributeU8.Data().empty(), false);
    EXPECT_EQ(attributeU8.Name(), u8_Single);
    EXPECT_EQ(attributeU8.Data().size(), numberOfElements);
    EXPECT_EQ(attributeU8.Type(), "uint8_t");

    ASSERT_EQ(attributeU16.Data().size() == 1, false);
    ASSERT_EQ(attributeU16.Data().empty(), false);
    EXPECT_EQ(attributeU16.Name(), u16_Single);
    EXPECT_EQ(attributeU16.Data().size(), numberOfElements);
    EXPECT_EQ(attributeU16.Type(), "uint16_t");

    ASSERT_EQ(attributeU32.Data().size() == 1, false);
    ASSERT_EQ(attributeU32.Data().empty(), false);
    EXPECT_EQ(attributeU32.Name(), u32_Single);
    EXPECT_EQ(attributeU32.Data().size(), numberOfElements);
    EXPECT_EQ(attributeU32.Type(), "uint32_t");

    ASSERT_EQ(attributeU64.Data().size() == 1, false);
    ASSERT_EQ(attributeU64.Data().empty(), false);
    EXPECT_EQ(attributeU64.Name(), u64_Single);
    EXPECT_EQ(attributeU64.Data().size(), numberOfElements);
    EXPECT_EQ(attributeU64.Type(), "uint64_t");
    EXPECT_EQ(sizeof(attributeU64.Data()[0]), 8);

    ASSERT_EQ(attributeFloat.Data().size() == 1, false);
    ASSERT_EQ(attributeFloat.Data().empty(), false);
    EXPECT_EQ(attributeFloat.Name(), float_Single);
    EXPECT_EQ(attributeFloat.Data().size(), numberOfElements);
    EXPECT_EQ(attributeFloat.Type(), "float");

    ASSERT_EQ(attributeDouble.Data().size() == 1, false);
    ASSERT_EQ(attributeDouble.Data().empty(), false);
    EXPECT_EQ(attributeDouble.Name(), double_Single);
    EXPECT_EQ(attributeDouble.Data().size(), numberOfElements);
    EXPECT_EQ(attributeDouble.Type(), "double");

    // Verify data
    for (size_t index = 0; index < numberOfElements; index++)
    {
        EXPECT_EQ(attributeI8.Data()[index], currentTestData.I8.at(index));
        EXPECT_EQ(attributeI16.Data()[index], currentTestData.I16.at(index));
        EXPECT_EQ(attributeI32.Data()[index], currentTestData.I32.at(index));
        EXPECT_EQ(attributeU8.Data()[index], currentTestData.U8.at(index));
        EXPECT_EQ(attributeU16.Data()[index], currentTestData.U16.at(index));
        EXPECT_EQ(attributeU32.Data()[index], currentTestData.U32.at(index));
        EXPECT_EQ(attributeFloat.Data()[index], currentTestData.R32.at(index));
        EXPECT_EQ(attributeDouble.Data()[index], currentTestData.R64.at(index));
    }
}

TEST_F(ADIOSDefineAttributeTest, GetAttribute)
{
    int mpiRank = 0, mpiSize = 1;
    size_t numberOfElements = 10;
#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // Define unique data for each process
    SmallTestData currentTestData = generateNewSmallTestData(m_TestData, 0, mpiRank, mpiSize);

    std::string mpiRankString = std::to_string(mpiRank);
    std::string s3_Single = std::string("s3_Array_") + mpiRankString;
    std::string i8_Single = std::string("i8_Array_") + mpiRankString;
    std::string i16_Single = std::string("i16_Array_") + mpiRankString;
    std::string i32_Single = std::string("i32_Array_") + mpiRankString;
    std::string i64_Single = std::string("i64_Array_") + mpiRankString;
    std::string u8_Single = std::string("u8_Array_") + mpiRankString;
    std::string u16_Single = std::string("u16_Array_") + mpiRankString;
    std::string u32_Single = std::string("u32_Array_") + mpiRankString;
    std::string u64_Single = std::string("u64_Array_") + mpiRankString;
    std::string float_Single = std::string("float_Array_") + mpiRankString;
    std::string double_Single = std::string("double_Array_") + mpiRankString;

    // Define ADIOS global value
    {
        io.DefineAttribute<std::string>(s3_Single, currentTestData.S3.data(), 3);
        io.DefineAttribute<int8_t>(i8_Single, currentTestData.I8.data(), numberOfElements);
        io.DefineAttribute<int16_t>(i16_Single, currentTestData.I16.data(), numberOfElements);
        io.DefineAttribute<int32_t>(i32_Single, currentTestData.I32.data(), numberOfElements);
        io.DefineAttribute<int64_t>(i64_Single, currentTestData.I64.data(), numberOfElements);
        io.DefineAttribute<uint8_t>(u8_Single, currentTestData.U8.data(), numberOfElements);
        io.DefineAttribute<uint16_t>(u16_Single, currentTestData.U16.data(), numberOfElements);
        io.DefineAttribute<uint32_t>(u32_Single, currentTestData.U32.data(), numberOfElements);
        io.DefineAttribute<uint64_t>(u64_Single, currentTestData.U64.data(), numberOfElements);
        io.DefineAttribute<float>(float_Single, currentTestData.R32.data(), numberOfElements);
        io.DefineAttribute<double>(double_Single, currentTestData.R64.data(), numberOfElements);
    }

    auto attributeS3 = io.InquireAttribute<std::string>(s3_Single);
    auto attributeI8 = io.InquireAttribute<int8_t>(i8_Single);
    auto attributeI16 = io.InquireAttribute<int16_t>(i16_Single);
    auto attributeI32 = io.InquireAttribute<int32_t>(i32_Single);
    auto attributeI64 = io.InquireAttribute<int64_t>(i64_Single);
    auto attributeU8 = io.InquireAttribute<uint8_t>(u8_Single);
    auto attributeU16 = io.InquireAttribute<uint16_t>(u16_Single);
    auto attributeU32 = io.InquireAttribute<uint32_t>(u32_Single);
    auto attributeU64 = io.InquireAttribute<uint64_t>(u64_Single);
    auto attributeFloat = io.InquireAttribute<float>(float_Single);
    auto attributeDouble = io.InquireAttribute<double>(double_Single);

    // Verify the return type is as expected
    ::testing::StaticAssertTypeEq<decltype(attributeS3), adios2::Attribute<std::string>>();
    ::testing::StaticAssertTypeEq<decltype(attributeI8), adios2::Attribute<int8_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeI16), adios2::Attribute<int16_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeI32), adios2::Attribute<int32_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeI64), adios2::Attribute<int64_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeU8), adios2::Attribute<uint8_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeU16), adios2::Attribute<uint16_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeU32), adios2::Attribute<uint32_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeU64), adios2::Attribute<uint64_t>>();
    ::testing::StaticAssertTypeEq<decltype(attributeFloat), adios2::Attribute<float>>();
    ::testing::StaticAssertTypeEq<decltype(attributeDouble), adios2::Attribute<double>>();

    // Verify the members are correct
    ASSERT_EQ(attributeS3.Data().size() == 1, false);
    ASSERT_EQ(attributeS3.Data().empty(), false);
    EXPECT_EQ(attributeS3.Name(), s3_Single);
    EXPECT_EQ(attributeS3.Data().size(), 3);
    EXPECT_EQ(attributeS3.Type(), "string");
    EXPECT_FALSE(attributeS3.IsValue());

    ASSERT_EQ(attributeI8.Data().size() == 1, false);
    ASSERT_EQ(attributeI8.Data().empty(), false);
    EXPECT_EQ(attributeI8.Name(), i8_Single);
    EXPECT_EQ(attributeI8.Data().size(), numberOfElements);
    EXPECT_EQ(attributeI8.Type(), "int8_t");
    EXPECT_FALSE(attributeI8.IsValue());

    ASSERT_EQ(attributeI16.Data().size() == 1, false);
    ASSERT_EQ(attributeI16.Data().empty(), false);
    EXPECT_EQ(attributeI16.Name(), i16_Single);
    EXPECT_EQ(attributeI16.Data().size(), numberOfElements);
    EXPECT_EQ(attributeI16.Type(), "int16_t");
    EXPECT_FALSE(attributeI16.IsValue());

    ASSERT_EQ(attributeI32.Data().size() == 1, false);
    ASSERT_EQ(attributeI32.Data().empty(), false);
    EXPECT_EQ(attributeI32.Name(), i32_Single);
    EXPECT_EQ(attributeI32.Data().size(), numberOfElements);
    EXPECT_EQ(attributeI32.Type(), "int32_t");
    EXPECT_FALSE(attributeI32.IsValue());

    ASSERT_EQ(attributeI64.Data().size() == 1, false);
    ASSERT_EQ(attributeI64.Data().empty(), false);
    EXPECT_EQ(attributeI64.Name(), i64_Single);
    EXPECT_EQ(attributeI64.Data().size(), numberOfElements);
    EXPECT_EQ(attributeI64.Type(), "int64_t");
    EXPECT_FALSE(attributeI64.IsValue());
    EXPECT_EQ(sizeof(attributeI64.Data()[0]), 8);

    ASSERT_EQ(attributeU8.Data().size() == 1, false);
    ASSERT_EQ(attributeU8.Data().empty(), false);
    EXPECT_EQ(attributeU8.Name(), u8_Single);
    EXPECT_EQ(attributeU8.Data().size(), numberOfElements);
    EXPECT_EQ(attributeU8.Type(), "uint8_t");
    EXPECT_FALSE(attributeU8.IsValue());

    ASSERT_EQ(attributeU16.Data().size() == 1, false);
    ASSERT_EQ(attributeU16.Data().empty(), false);
    EXPECT_EQ(attributeU16.Name(), u16_Single);
    EXPECT_EQ(attributeU16.Data().size(), numberOfElements);
    EXPECT_EQ(attributeU16.Type(), "uint16_t");
    EXPECT_FALSE(attributeU16.IsValue());

    ASSERT_EQ(attributeU32.Data().size() == 1, false);
    ASSERT_EQ(attributeU32.Data().empty(), false);
    EXPECT_EQ(attributeU32.Name(), u32_Single);
    EXPECT_EQ(attributeU32.Data().size(), numberOfElements);
    EXPECT_EQ(attributeU32.Type(), "uint32_t");
    EXPECT_FALSE(attributeU32.IsValue());

    ASSERT_EQ(attributeU64.Data().size() == 1, false);
    ASSERT_EQ(attributeU64.Data().empty(), false);
    EXPECT_EQ(attributeU64.Name(), u64_Single);
    EXPECT_EQ(attributeU64.Data().size(), numberOfElements);
    EXPECT_EQ(attributeU64.Type(), "uint64_t");
    EXPECT_FALSE(attributeU64.IsValue());
    EXPECT_EQ(sizeof(attributeU64.Data()[0]), 8);

    ASSERT_EQ(attributeFloat.Data().size() == 1, false);
    ASSERT_EQ(attributeFloat.Data().empty(), false);
    EXPECT_EQ(attributeFloat.Name(), float_Single);
    EXPECT_EQ(attributeFloat.Data().size(), numberOfElements);
    EXPECT_EQ(attributeFloat.Type(), "float");
    EXPECT_FALSE(attributeFloat.IsValue());

    ASSERT_EQ(attributeDouble.Data().size() == 1, false);
    ASSERT_EQ(attributeDouble.Data().empty(), false);
    EXPECT_EQ(attributeDouble.Name(), double_Single);
    EXPECT_EQ(attributeDouble.Data().size(), numberOfElements);
    EXPECT_EQ(attributeDouble.Type(), "double");
    EXPECT_FALSE(attributeDouble.IsValue());

    // Verify data
    for (size_t index = 0; index < numberOfElements; index++)
    {
        EXPECT_EQ(attributeI8.Data()[index], currentTestData.I8.at(index));
        EXPECT_EQ(attributeI16.Data()[index], currentTestData.I16.at(index));
        EXPECT_EQ(attributeI32.Data()[index], currentTestData.I32.at(index));
        EXPECT_EQ(attributeU8.Data()[index], currentTestData.U8.at(index));
        EXPECT_EQ(attributeU16.Data()[index], currentTestData.U16.at(index));
        EXPECT_EQ(attributeU32.Data()[index], currentTestData.U32.at(index));
        EXPECT_EQ(attributeFloat.Data()[index], currentTestData.R32.at(index));
        EXPECT_EQ(attributeDouble.Data()[index], currentTestData.R64.at(index));
    }
}

TEST_F(ADIOSDefineAttributeTest, DefineAndRemove)
{
    auto lf_CheckRemove = [&](const std::string attributeName) {
        const bool isRemoved = io.RemoveAttribute(attributeName);
        EXPECT_EQ(isRemoved, true);
    };

    const adios2::Dims shape = {10};
    const adios2::Dims start = {0};
    const adios2::Dims count = {10};

    io.DefineAttribute<std::string>("iString", "String Attribute");
    io.DefineAttribute<int8_t>("i8", -8);
    io.DefineAttribute<int16_t>("i16", -16);
    io.DefineAttribute<int32_t>("i32", -32);
    io.DefineAttribute<int64_t>("i64", -64);
    io.DefineAttribute<uint8_t>("u8", 8);
    io.DefineAttribute<uint16_t>("u16", 16);
    io.DefineAttribute<uint32_t>("u32", 32);
    io.DefineAttribute<uint64_t>("u64", 64);
    io.DefineAttribute<float>("r32", 32);
    io.DefineAttribute<double>("r64", 64);

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

    auto attr_iString = io.InquireAttribute<std::string>("iString");
    auto attr_i8 = io.InquireAttribute<int8_t>("i8");
    auto attr_i16 = io.InquireAttribute<int16_t>("i16");
    auto attr_i32 = io.InquireAttribute<int32_t>("i32");
    auto attr_i64 = io.InquireAttribute<int64_t>("i64");
    auto attr_u8 = io.InquireAttribute<uint8_t>("u8");
    auto attr_u16 = io.InquireAttribute<uint16_t>("u16");
    auto attr_u32 = io.InquireAttribute<uint32_t>("u32");
    auto attr_u64 = io.InquireAttribute<uint64_t>("u64");
    auto attr_r32 = io.InquireAttribute<float>("r32");
    auto attr_r64 = io.InquireAttribute<double>("r64");

    EXPECT_FALSE(attr_iString);
    EXPECT_FALSE(attr_i8);
    EXPECT_FALSE(attr_i16);
    EXPECT_FALSE(attr_i32);
    EXPECT_FALSE(attr_i64);
    EXPECT_FALSE(attr_u8);
    EXPECT_FALSE(attr_u16);
    EXPECT_FALSE(attr_u32);
    EXPECT_FALSE(attr_u64);
    EXPECT_FALSE(attr_r32);
    EXPECT_FALSE(attr_r64);
}

TEST_F(ADIOSDefineAttributeTest, DefineRemoveDefine)
{
    auto lf_CheckRemove = [&](const std::string attributeName) {
        const bool isRemoved = io.RemoveAttribute(attributeName);
        EXPECT_EQ(isRemoved, true);
    };

    io.DefineAttribute<std::string>("string_0", "attribute_0");
    io.DefineAttribute<std::string>("string_1", "attribute_1");

    lf_CheckRemove("string_0");

    std::array<adios2::Attribute<std::string>, 2> attributes;

    attributes[0] = io.InquireAttribute<std::string>("string_0");
    EXPECT_FALSE(attributes[0]);

    attributes[1] = io.InquireAttribute<std::string>("string_1");
    EXPECT_TRUE(attributes[1]);
    EXPECT_EQ(attributes[1].Name(), "string_1");
    EXPECT_EQ(attributes[1].Data().front(), "attribute_1");

    io.DefineAttribute<std::string>("string_0", "attribute_0_new");

    // check again after defining variable
    attributes[1] = io.InquireAttribute<std::string>("string_1");
    EXPECT_TRUE(attributes[1]);
    EXPECT_EQ(attributes[1].Name(), "string_1");
    EXPECT_EQ(attributes[1].Data().front(), "attribute_1");

    attributes[0] = io.InquireAttribute<std::string>("string_0");
    EXPECT_TRUE(attributes[0]);
    EXPECT_EQ(attributes[0].Name(), "string_0");
    EXPECT_EQ(attributes[0].Data().front(), "attribute_0_new");

    auto attribute2 = io.DefineAttribute<std::string>("string_2", "attribute_2");
    EXPECT_TRUE(attribute2);
    EXPECT_EQ(attribute2.Name(), "string_2");
    EXPECT_EQ(attribute2.Data().front(), "attribute_2");
}

TEST_F(ADIOSDefineAttributeTest, DefineAndRemoveAll)
{
    const adios2::Dims shape = {10};
    const adios2::Dims start = {0};
    const adios2::Dims count = {10};

    io.DefineAttribute<std::string>("iString", "String Attribute");
    io.DefineAttribute<int8_t>("i8", -8);
    io.DefineAttribute<int16_t>("i16", -16);
    io.DefineAttribute<int32_t>("i32", -32);
    io.DefineAttribute<int64_t>("i64", -64);
    io.DefineAttribute<uint8_t>("u8", 8);
    io.DefineAttribute<uint16_t>("u16", 16);
    io.DefineAttribute<uint32_t>("u32", 32);
    io.DefineAttribute<uint64_t>("u64", 64);
    io.DefineAttribute<float>("r32", 32);
    io.DefineAttribute<double>("r64", 64);

    io.RemoveAllAttributes();

    auto attr_iString = io.InquireAttribute<std::string>("iString");
    auto attr_i8 = io.InquireAttribute<int8_t>("i8");
    auto attr_i16 = io.InquireAttribute<int16_t>("i16");
    auto attr_i32 = io.InquireAttribute<int32_t>("i32");
    auto attr_i64 = io.InquireAttribute<int64_t>("i64");
    auto attr_u8 = io.InquireAttribute<uint8_t>("u8");
    auto attr_u16 = io.InquireAttribute<uint16_t>("u16");
    auto attr_u32 = io.InquireAttribute<uint32_t>("u32");
    auto attr_u64 = io.InquireAttribute<uint64_t>("u64");
    auto attr_r32 = io.InquireAttribute<float>("r32");
    auto attr_r64 = io.InquireAttribute<double>("r64");

    EXPECT_FALSE(attr_iString);
    EXPECT_FALSE(attr_i8);
    EXPECT_FALSE(attr_i16);
    EXPECT_FALSE(attr_i32);
    EXPECT_FALSE(attr_i64);
    EXPECT_FALSE(attr_u8);
    EXPECT_FALSE(attr_u16);
    EXPECT_FALSE(attr_u32);
    EXPECT_FALSE(attr_u64);
    EXPECT_FALSE(attr_r32);
    EXPECT_FALSE(attr_r64);
}

TEST_F(ADIOSDefineAttributeTest, DefineCheckType)
{
    const adios2::Dims shape = {10};
    const adios2::Dims start = {0};
    const adios2::Dims count = {10};

    io.DefineAttribute<std::string>("iString", "String Attribute");
    io.DefineAttribute<int8_t>("i8", -8);
    io.DefineAttribute<int16_t>("i16", -16);
    io.DefineAttribute<int32_t>("i32", -32);
    io.DefineAttribute<int64_t>("i64", -64);
    io.DefineAttribute<uint8_t>("u8", 8);
    io.DefineAttribute<uint16_t>("u16", 16);
    io.DefineAttribute<uint32_t>("u32", 32);
    io.DefineAttribute<uint64_t>("u64", 64);
    io.DefineAttribute<float>("r32", 32);
    io.DefineAttribute<double>("r64", 64);

    EXPECT_EQ(io.AttributeType("iString"), adios2::GetType<std::string>());
    EXPECT_EQ(io.AttributeType("i8"), adios2::GetType<int8_t>());
    EXPECT_EQ(io.AttributeType("i16"), adios2::GetType<int16_t>());
    EXPECT_EQ(io.AttributeType("i32"), adios2::GetType<int32_t>());
    EXPECT_EQ(io.AttributeType("i64"), adios2::GetType<int64_t>());
    EXPECT_EQ(io.AttributeType("u8"), adios2::GetType<uint8_t>());
    EXPECT_EQ(io.AttributeType("u16"), adios2::GetType<uint16_t>());
    EXPECT_EQ(io.AttributeType("u32"), adios2::GetType<uint32_t>());
    EXPECT_EQ(io.AttributeType("u64"), adios2::GetType<uint64_t>());
    EXPECT_EQ(io.AttributeType("r32"), adios2::GetType<float>());
    EXPECT_EQ(io.AttributeType("r64"), adios2::GetType<double>());
}

TEST_F(ADIOSDefineAttributeTest, VariableException)
{
    const std::string separator = "/";

    std::vector<std::string> numbers = {"one", "two", "three"};

// Write test data using BP
#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        EXPECT_THROW(io.DefineAttribute<std::string>("Hello Value", "Value", "myVar1"),
                     std::invalid_argument);

        EXPECT_THROW(io.DefineAttribute<std::string>("Hello Array", numbers.data(), numbers.size(),
                                                     "myVar1", separator),
                     std::invalid_argument);

        io.DefineVariable<int>("myVar1");

        EXPECT_NO_THROW(io.DefineAttribute<std::string>("Hello Value", "Value", "myVar1"));

        EXPECT_NO_THROW(io.DefineAttribute<std::string>("Hello Array", numbers.data(),
                                                        numbers.size(), "myVar1", separator));
    }
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

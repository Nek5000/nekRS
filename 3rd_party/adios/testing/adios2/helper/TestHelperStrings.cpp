/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cmath>
#include <cstdint>
#include <cstring>

#include <iostream>
#include <limits>
#include <stdexcept>

#include <adios2.h>
#include <adios2/common/ADIOSTypes.h>
#include <adios2/helper/adiosString.h>

#include <gtest/gtest.h>

TEST(ADIOS2HelperString, ADIOS2HelperStringFNF)
{

    const std::string fname("nosuchfile.txt");
    const std::string hint("");

    ASSERT_THROW(adios2::helper::FileToString(fname, hint), std::ios_base::failure);
}

TEST(ADIOS2HelperString, ADIOS2HelperStringParameterMapFromVector)
{

    const std::vector<std::string> badparam_in = {"badparam"};
    const std::vector<std::string> emptyparam_in = {"emptyparam="};
    const std::vector<std::string> dupparam_in = {"dupparam=1", "dupparam=2"};
    const std::vector<std::string> param_in = {"param1=1", "param2=2", "param3=3"};

    adios2::Params parameters = adios2::helper::BuildParametersMap(param_in, '=');

    ASSERT_EQ(parameters.find("param1")->second, "1");
    ASSERT_EQ(parameters.find("param2")->second, "2");
    ASSERT_EQ(parameters.find("param3")->second, "3");

    ASSERT_THROW(adios2::helper::BuildParametersMap(badparam_in, '='), std::invalid_argument);
    ASSERT_THROW(adios2::helper::BuildParametersMap(emptyparam_in, '='), std::invalid_argument);
    ASSERT_THROW(adios2::helper::BuildParametersMap(dupparam_in, '='), std::invalid_argument);
}

TEST(ADIOS2HelperString, ADIOS2HelperStringParameterMapFromString)
{

    const std::string badparam_in = "badparam";
    const std::string emptyparam_in = "emptyparam=";
    const std::string dupparam_in = "dupparam = 1 , dupparam=2";
    const std::string param_in = "param1=1, param2=2,                                       "
                                 "                    param3=3";

    adios2::Params parameters = adios2::helper::BuildParametersMap(param_in, '=', ',');

    ASSERT_EQ(parameters.find("param1")->second, "1");
    ASSERT_EQ(parameters.find("param2")->second, "2");
    ASSERT_EQ(parameters.find("param3")->second, "3");

    ASSERT_THROW(adios2::helper::BuildParametersMap(badparam_in, '=', ','), std::invalid_argument);
    ASSERT_THROW(adios2::helper::BuildParametersMap(emptyparam_in, '=', ','),
                 std::invalid_argument);
    ASSERT_THROW(adios2::helper::BuildParametersMap(dupparam_in, '=', ','), std::invalid_argument);
}

TEST(ADIOS2HelperString, ADIOS2HelperStringAddExtension)
{

    const std::string abc("abc");
    const std::string abcbp("abc.bp");
    const std::string ext(".bp");

    ASSERT_EQ(adios2::helper::AddExtension(abc, ext), abcbp);
    ASSERT_EQ(adios2::helper::AddExtension(abcbp, ext), abcbp);
}

TEST(ADIOS2HelperString, ADIOS2HelperStringEndsWithCaseSensitive)
{

    const std::string abcdefgh("abcd.efgh");
    const std::string theend(".efgh");
    const std::string notanend(".ephs");
    const std::string shortstr("abc");

    ASSERT_TRUE(adios2::helper::EndsWith(abcdefgh, theend));
    ASSERT_FALSE(adios2::helper::EndsWith(abcdefgh, notanend));
    ASSERT_FALSE(adios2::helper::EndsWith(shortstr, theend));
}

TEST(ADIOS2HelperString, ADIOS2HelperStringEndsWithCaseInsensitive)
{

    const std::string abcdefgh("abCd.eFgH");
    const std::string end1(".efGh");
    const std::string end2(".efgh");
    const std::string noend1(".ephs");
    const std::string noend2(".efdHs");
    const std::string shortstr("ABC");

    ASSERT_TRUE(adios2::helper::EndsWith(abcdefgh, end1, false));
    ASSERT_TRUE(adios2::helper::EndsWith(abcdefgh, end2, false));
    ASSERT_FALSE(adios2::helper::EndsWith(abcdefgh, noend1, false));
    ASSERT_FALSE(adios2::helper::EndsWith(abcdefgh, noend2, false));
    ASSERT_FALSE(adios2::helper::EndsWith(shortstr, end1, false));
}

TEST(ADIOS2HelperString, ADIOS2HelperStringConversion)
{

    const std::string dbl("123.1230");
    const std::string uint("123");
    const std::string notnum("notnum");
    const std::string hint("");

    const double diff = std::abs(adios2::helper::StringTo<double>(dbl, hint) - 123.123);
    ASSERT_LT(diff, 1E-4);

    ASSERT_THROW(adios2::helper::StringTo<double>(notnum, hint), std::invalid_argument);
    ASSERT_EQ(adios2::helper::StringTo<uint32_t>(uint, hint), 123);
    ASSERT_THROW(adios2::helper::StringTo<uint32_t>(notnum, hint), std::invalid_argument);
}

TEST(ADIOS2HelperString, ADIOS2HelperDimString)
{

    const adios2::Dims dimensions = {1, 2, 3};
    const std::string dimstr("Dims(3):[1, 2, 3]");

    ASSERT_EQ(adios2::helper::DimsToString(dimensions), dimstr);
}

TEST(ADIOS2HelperString, ADIOS2HelperGlobalName)
{

    const std::string localName("myfile.bp");
    const std::string prefix("mydir");
    const std::string separator("/");
    const std::string global("mydir/myfile.bp");

    ASSERT_EQ(adios2::helper::GlobalName(localName, prefix, separator), global);
}

int main(int argc, char **argv)
{

    int result;
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();

    return result;
}

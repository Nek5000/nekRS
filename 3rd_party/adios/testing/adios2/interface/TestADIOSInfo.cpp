#include <adios2/core/Info.h>
#include <gtest/gtest.h>

TEST(ADIOSInterface, info_available_features)
{
    size_t nfeatures = 0;
    const char *const *list_features = nullptr;
    adios2_available_features(&nfeatures, &list_features);

    EXPECT_GE(nfeatures, 0);
    EXPECT_NE(list_features, nullptr);
    EXPECT_EQ(list_features[nfeatures], nullptr);
}

TEST(ADIOSInterface, info_available_engines)
{
    size_t nengines = 0;
    const char *const *list_engines = nullptr;
    adios2_available_engines(&nengines, &list_engines);

    EXPECT_GE(nengines, 1);
    EXPECT_NE(list_engines, nullptr);
    EXPECT_EQ(list_engines[nengines], nullptr);
}

TEST(ADIOSInterface, info_available_operators)
{
    size_t noperators = 0;
    const char *const *list_operators = nullptr;
    adios2_available_operators(&noperators, &list_operators);

    EXPECT_GE(noperators, 0);
    EXPECT_NE(list_operators, nullptr);
    EXPECT_EQ(list_operators[noperators], nullptr);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

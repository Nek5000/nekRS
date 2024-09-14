#include <adios2/helper/adiosRangeFilter.h>
#include <fstream>
#include <gtest/gtest.h>

TEST(ADIOS2RangeFilter, Elements)
{
    adios2::helper::RangeFilter r;

    r.ParseSelection("2");
    EXPECT_FALSE(r.IsSelected(0));
    EXPECT_FALSE(r.IsSelected(1));
    EXPECT_TRUE(r.IsSelected(2));
    EXPECT_FALSE(r.IsSelected(3));

    r.ParseSelection("2 3 0");
    EXPECT_TRUE(r.IsSelected(0));
    EXPECT_FALSE(r.IsSelected(1));
    EXPECT_TRUE(r.IsSelected(2));
    EXPECT_TRUE(r.IsSelected(3));
}

TEST(ADIOS2RangeFilter, Ranges)
{
    adios2::helper::RangeFilter r;

    r.ParseSelection("1:2");
    EXPECT_FALSE(r.IsSelected(0));
    EXPECT_TRUE(r.IsSelected(1));
    EXPECT_TRUE(r.IsSelected(2));
    EXPECT_FALSE(r.IsSelected(3));

    r.ParseSelection("1:2:1");
    EXPECT_FALSE(r.IsSelected(0));
    EXPECT_TRUE(r.IsSelected(1));
    EXPECT_TRUE(r.IsSelected(2));
    EXPECT_FALSE(r.IsSelected(3));

    r.ParseSelection("1:2:2");
    EXPECT_FALSE(r.IsSelected(0));
    EXPECT_TRUE(r.IsSelected(1));
    EXPECT_FALSE(r.IsSelected(2));
    EXPECT_FALSE(r.IsSelected(3));

    r.ParseSelection("1:6:3");
    EXPECT_FALSE(r.IsSelected(0));
    EXPECT_TRUE(r.IsSelected(1));
    EXPECT_FALSE(r.IsSelected(2));
    EXPECT_FALSE(r.IsSelected(3));
    EXPECT_TRUE(r.IsSelected(4));
    EXPECT_FALSE(r.IsSelected(5));
    EXPECT_FALSE(r.IsSelected(6));

    r.ParseSelection("1:3:2 0:4:2");
    EXPECT_TRUE(r.IsSelected(0));
    EXPECT_TRUE(r.IsSelected(1));
    EXPECT_TRUE(r.IsSelected(2));
    EXPECT_TRUE(r.IsSelected(3));
    EXPECT_TRUE(r.IsSelected(4));
    EXPECT_FALSE(r.IsSelected(5));
    EXPECT_FALSE(r.IsSelected(6));
}

TEST(ADIOS2RangeFilter, UnlimitedRanges)
{
    adios2::helper::RangeFilter r;

    r.ParseSelection("2:n");
    EXPECT_FALSE(r.IsSelected(0));
    EXPECT_FALSE(r.IsSelected(1));
    EXPECT_TRUE(r.IsSelected(2));
    EXPECT_TRUE(r.IsSelected(3000));

    r.ParseSelection("2:N:1");
    EXPECT_FALSE(r.IsSelected(0));
    EXPECT_FALSE(r.IsSelected(1));
    EXPECT_TRUE(r.IsSelected(2));
    EXPECT_TRUE(r.IsSelected(3000));

    r.ParseSelection("1:n:2");
    EXPECT_FALSE(r.IsSelected(0));
    EXPECT_TRUE(r.IsSelected(1));
    EXPECT_FALSE(r.IsSelected(2));
    EXPECT_TRUE(r.IsSelected(3));
    EXPECT_FALSE(r.IsSelected(3000));
    EXPECT_TRUE(r.IsSelected(3001));

    r.ParseSelection("1:n:2   0:N:2");
    EXPECT_TRUE(r.IsSelected(0));
    EXPECT_TRUE(r.IsSelected(1));
    EXPECT_TRUE(r.IsSelected(2));
    EXPECT_TRUE(r.IsSelected(3));
    EXPECT_TRUE(r.IsSelected(4));
    EXPECT_TRUE(r.IsSelected(5));
    EXPECT_TRUE(r.IsSelected(6));
    EXPECT_TRUE(r.IsSelected(3000));
    EXPECT_TRUE(r.IsSelected(3001));
}

TEST(ADIOS2RangeFilter, InvalidSelection)
{
    adios2::helper::RangeFilter r;

    EXPECT_THROW(r.ParseSelection("b"), std::invalid_argument);
    EXPECT_THROW(r.ParseSelection("n"), std::invalid_argument);
    EXPECT_THROW(r.ParseSelection("1:n:2,3"), std::invalid_argument);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

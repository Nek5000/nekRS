/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <array>
#include <limits>
#include <stdexcept>
#include <vector>

#include <adios2.h>
#include <adios2/common/ADIOSTypes.h>

#include <gtest/gtest.h>

std::string engineName; // comes from command line

class ADIOSSelectionColumnMajorTest : public ::testing::Test
{
public:
    ADIOSSelectionColumnMajorTest() = default;
};

void copySelection2D_CM(const double *a, const adios2::Dims &shape, const adios2::Dims &start,
                        const adios2::Dims &count, double *b)
{
    std::cout << " copy Block: shape = " << adios2::ToString(shape)
              << " start = " << adios2::ToString(start) << " count = " << adios2::ToString(count)
              << std::endl;
    double *bp = b;
    for (size_t y = 0; y < count[1]; ++y)
    {
        size_t aidx = (start[1] + y) * shape[0] + start[0];
        for (size_t x = 0; x < count[0]; ++x)
        {
            /*if (x == 0 && y == 0)
            {
                std::cout << " idx of pos = " << adios2::ToString({x, y})
                          << " => " << aidx << ", a = " << a[aidx]
                          << std::endl;
            }*/
            *bp = a[aidx];
            ++bp;
            ++aidx;
        }
    }
}

/*bool compareSelection2D(const double *a, const adios2::Dims &shape,
                        const adios2::Dims &start, const adios2::Dims &count,
                        double *b, adios2::Dims &firstNonEqPoint)
{
    std::cout << " compare Block: shape = " << adios2::ToString(shape)
              << " start = " << adios2::ToString(start)
              << " count = " << adios2::ToString(count) << std::endl;
    double *bp = b;
    for (size_t x = 0; x < count[0]; ++x)
    {
        for (size_t y = 0; y < count[1]; ++y)
        {
            const size_t aidx = (start[0] + x) * shape[1] + start[1] + y;

            if (*bp != a[aidx])
            {
                firstNonEqPoint = {x, y};
                std::cout << "   Non-match at pos = "
                          << adios2::ToString(firstNonEqPoint)
                          << " : a = " << a[aidx] << ", b = " << *bp
                          << std::endl;
                return false;
            }
            ++bp;
        }
    }
    return true;
}*/

bool compareSelection2D_RM(const double *a, const adios2::Dims &shape, const adios2::Dims &start,
                           const adios2::Dims &count, double *b, adios2::Dims &firstNonEqPoint)
{
    std::cout << " compare Block: shape = " << adios2::ToString(shape)
              << " start = " << adios2::ToString(start) << " count = " << adios2::ToString(count)
              << std::endl;
    bool match = true;
    double *bp = b;
    for (size_t x = 0; x < count[0]; ++x)
    {
        size_t aidx = (start[0] + x) * shape[1] + start[1];
        for (size_t y = 0; y < count[1]; ++y)
        {
            if (*bp != a[aidx])
            {
                firstNonEqPoint = {x, y};
                std::cout << "   Non-match at pos = " << adios2::ToString(firstNonEqPoint)
                          << " : a = " << a[aidx] << ", b = " << *bp << std::endl;
                match = false;
            }
            ++bp;
            ++aidx;
        }
    }
    return match;
}

bool compareSelection2D_CM(const double *a, const adios2::Dims &shape, const adios2::Dims &start,
                           const adios2::Dims &count, double *b, adios2::Dims &firstNonEqPoint)
{
    std::cout << " compare Block: shape = " << adios2::ToString(shape)
              << " start = " << adios2::ToString(start) << " count = " << adios2::ToString(count)
              << std::endl;
    bool match = true;
    double *bp = b;
    for (size_t y = 0; y < count[1]; ++y)
    {
        size_t aidx = (start[1] + y) * shape[0] + start[0];
        for (size_t x = 0; x < count[0]; ++x)
        {
            if (*bp != a[aidx])
            {
                firstNonEqPoint = {y, x};
                std::cout << "   Non-match at pos = " << adios2::ToString(firstNonEqPoint)
                          << " : a = " << a[aidx] << ", b = " << *bp << std::endl;
                match = false;
            }
            ++bp;
            ++aidx;
        }
    }
    return match;
}

TEST_F(ADIOSSelectionColumnMajorTest, 2D)
{
    /* !!! In this test DIM1 (C1) is the "fastest" dimension as in Fortran */

    constexpr size_t C1 = 5;
    constexpr size_t C2 = 4;
    constexpr size_t DIM1 = 3 * C1;
    constexpr size_t DIM2 = 3 * C2;
    const std::string filename = "ADIOSSelectionColumnMajor2D.bp";
    double a[DIM1 * DIM2];

    double *ap = a;
    for (size_t y = 1; y <= DIM2; ++y)
    {
        for (size_t x = 1; x <= DIM1; ++x)
        {
            *ap = y * 1.0 + x / 100.0;
            ++ap;
        }
    }
    /* ./bin/bpls testing/adios2/engine/bp/bp4/ADIOSSelectionColumnMajor2D.bp
        -d  a -f "%5.2f " -n 15
    */

    // Write test data as if we were in Fortran (Column-major language)
    {

        adios2::ADIOS adiosF("", "Fortran");
        adios2::IO ioWrite = adiosF.DeclareIO("TestIOWrite");
        if (!engineName.empty())
        {
            ioWrite.SetEngine(engineName);
        }
        else
        {
            // Create the BP Engine
            ioWrite.SetEngine("BPFile");
        }

        std::cout << "Write data as Column-major..." << std::endl;

        adios2::Engine engine = ioWrite.Open(filename, adios2::Mode::Write);
        const adios2::Dims shape = {DIM1, DIM2};
        const adios2::Dims count = {C1, C2};
        adios2::Dims start{0, 0};

        double b[C1 * C2];

        adios2::Variable<double> var = ioWrite.DefineVariable<double>("a", shape, start, count);

        // adios2::Variable<double> vara1 =
        //    ioWrite.DefineVariable<double>("a1", shape, start, shape);

        engine.BeginStep();
        // engine.Put(vara1, a);

        for (size_t y = 0; y < DIM2; y += count[1])
        {
            for (size_t x = 0; x < DIM1; x += count[0])
            {
                start = {x, y};
                copySelection2D_CM(a, shape, start, count, b);
                var.SetSelection({start, count});
                engine.Put(var, b, adios2::Mode::Sync);
            }
        }

        engine.EndStep();
        engine.Close();
    }

    // Read data pretending to be a Column-major language
    {
        adios2::ADIOS adiosF("", "Fortran");
        adios2::IO ioRead = adiosF.DeclareIO("TestIOReadF");
        ioRead.SetEngine("File");
        adios2::Engine engine = ioRead.Open(filename, adios2::Mode::Read);
        EXPECT_TRUE(engine);
        engine.BeginStep();
        adios2::Variable<double> var = ioRead.InquireVariable<double>("a");
        EXPECT_TRUE(var);
        const adios2::Dims shapeF = {DIM1, DIM2};
        const adios2::Dims countF = {C1, C2};
        adios2::Dims s{0, 0};
        adios2::Dims c = shapeF;
        adios2::Dims firstNonMatch{0, 0, 0};
        std::vector<double> res;

        std::cout << "Column-major read selections with entire blocks..." << std::endl;

        /* Entire array */
        s = {0, 0};
        c = shapeF;
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), DIM1 * DIM2);
        EXPECT_TRUE(compareSelection2D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* First block */
        s = {0, 0};
        c = countF;
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* Single block in the center */
        s = {5, 4};
        c = countF;
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* Two blocks in X (fast) direction */
        s = {5, 4};
        c = {2 * countF[0], countF[1]};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* Three blocks in Y (slow) direction */
        s = {5, 0};
        c = {countF[0], 3 * countF[1]};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* Four blocks in X-Y direction */
        s = {5, 4};
        c = {2 * countF[0], 2 * countF[1]};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        std::cout << "Column-major read selections with partial blocks..." << std::endl;

        /* center part of single block in center */
        s = {6, 5};
        c = {countF[0] - 2, countF[1] - 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into two blocks in X direction */
        s = {6, 5};
        c = {2 * countF[0] - 2, countF[1] - 1};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into two blocks in Y direction */
        s = {6, 5};
        c = {countF[0] - 1, 2 * countF[1] - 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into three blocks in Y direction */
        s = {6, 1};
        c = {countF[0] - 1, 3 * countF[1] - 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into four blocks in X-Y direction */
        s = {6, 1};
        c = {2 * countF[0] - 2, 2 * countF[1] - 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* Center block plus 1 in each direction, cutting into all blocks */
        /* ./bin/bpls -la testing/adios2/engine/bp/bp4/ADIOSSelection2D.bp
            -d a  -f "%6.3f " -s "4,3" -c "7,6" -n 6
        */
        s = {4, 3};
        c = {countF[0] + 2, countF[1] + 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        engine.EndStep();
        engine.Close();
    }

    // Read data with Row Major (as in C++)
    {
        adios2::ADIOS adios;
        adios2::IO ioRead = adios.DeclareIO("TestIORead");
        ioRead.SetEngine("File");
        adios2::Engine engine = ioRead.Open(filename, adios2::Mode::Read);
        EXPECT_TRUE(engine);
        engine.BeginStep();
        adios2::Variable<double> var = ioRead.InquireVariable<double>("a");
        EXPECT_TRUE(var);
        const adios2::Dims shape = {DIM2, DIM1};
        const adios2::Dims count = {C2, C1};
        adios2::Dims s{0, 0};
        adios2::Dims c = shape;
        adios2::Dims firstNonMatch{0, 0};
        std::vector<double> res;

        std::cout << "Row-major read selections with entire blocks..." << std::endl;

        /* Entire array */
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), DIM1 * DIM2);
        EXPECT_TRUE(compareSelection2D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* Single block in the center */
        s = {4, 5};
        c = count;
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* Two blocks in X direction */
        s = {4, 5};
        c = {2 * count[0], count[1]};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* Three blocks in Y direction */
        s = {4, 0};
        c = {count[0], 3 * count[1]};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* Four blocks in X-Y direction */
        s = {4, 5};
        c = {2 * count[0], 2 * count[1]};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /*
         *   Partial blocks
         */

        std::cout << "Row-major read selections with partial blocks..." << std::endl;

        /* center part of single block in center */
        s = {5, 6};
        c = {count[0] - 2, count[1] - 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into two blocks in X direction */
        s = {5, 6};
        c = {2 * count[0] - 2, count[1] - 1};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into two blocks in Y direction */
        s = {5, 6};
        c = {count[0] - 1, 2 * count[1] - 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into three blocks in Y direction */
        s = {5, 1};
        c = {count[0] - 1, 3 * count[1] - 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into four blocks in X-Y direction */
        s = {1, 6};
        c = {2 * count[0] - 2, 2 * count[1] - 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* Center block plus 1 in each direction, cutting into all blocks */
        /* ./bin/bpls -la testing/adios2/engine/bp/bp4/ADIOSSelection2D.bp
            -d a  -f "%6.3f " -s "4,3" -c "7,6" -n 6
        */
        s = {3, 4};
        c = {count[0] + 2, count[1] + 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1]);
        EXPECT_TRUE(compareSelection2D_RM(a, shape, s, c, res.data(), firstNonMatch));

        engine.EndStep();
        engine.Close();
    }
}

void copySelection3D_CM(const double *a, const adios2::Dims &shape, const adios2::Dims &start,
                        const adios2::Dims &count, double *b)
{
    /*std::cout << " copy Block: shape = " << adios2::ToString(shape)
              << " start = " << adios2::ToString(start)
              << " count = " << adios2::ToString(count) << std::endl;*/
    double *bp = b;
    for (size_t z = 0; z < count[2]; ++z)
    {
        for (size_t y = 0; y < count[1]; ++y)
        {
            size_t aidx =
                (start[2] + z) * shape[0] * shape[1] + (start[1] + y) * shape[0] + start[0];
            for (size_t x = 0; x < count[0]; ++x)
            {
                *bp = a[aidx];
                ++bp;
                ++aidx;
            }
        }
    }
}

bool compareSelection3D_RM(const double *a, const adios2::Dims &shape, const adios2::Dims &start,
                           const adios2::Dims &count, double *b, adios2::Dims &firstNonEqPoint)
{
    std::cout << " compare Block: shape = " << adios2::ToString(shape)
              << " start = " << adios2::ToString(start) << " count = " << adios2::ToString(count)
              << std::endl;
    double *bp = b;
    for (size_t x = 0; x < count[0]; ++x)
    {
        for (size_t y = 0; y < count[1]; ++y)
        {
            size_t aidx =
                (start[0] + x) * shape[1] * shape[2] + (start[1] + y) * shape[2] + start[2];
            for (size_t z = 0; z < count[2]; ++z)
            {
                if (*bp != a[aidx])
                {
                    firstNonEqPoint = {x, y, z};
                    std::cout << "   Non-match at pos = " << adios2::ToString(firstNonEqPoint)
                              << " : a = " << a[aidx] << ", b = " << *bp << std::endl;
                    return false;
                }
                ++bp;
                ++aidx;
            }
        }
    }
    return true;
}

bool compareSelection3D_CM(const double *a, const adios2::Dims &shape, const adios2::Dims &start,
                           const adios2::Dims &count, double *b, adios2::Dims &firstNonEqPoint)
{
    std::cout << " compare Block: shape = " << adios2::ToString(shape)
              << " start = " << adios2::ToString(start) << " count = " << adios2::ToString(count)
              << std::endl;
    double *bp = b;
    for (size_t z = 0; z < count[2]; ++z)
    {
        for (size_t y = 0; y < count[1]; ++y)
        {
            size_t aidx =
                (start[2] + z) * shape[0] * shape[1] + (start[1] + y) * shape[0] + start[0];
            for (size_t x = 0; x < count[0]; ++x)
            {
                if (*bp != a[aidx])
                {
                    firstNonEqPoint = {x, y, z};
                    std::cout << "   Non-match at pos = " << adios2::ToString(firstNonEqPoint)
                              << " : a = " << a[aidx] << ", b = " << *bp << std::endl;
                    return false;
                }
                ++bp;
                ++aidx;
            }
        }
    }
    return true;
}

TEST_F(ADIOSSelectionColumnMajorTest, 3D)
{
    /* !!! In this test DIM1 (C1) is the "fastest" dimension as in Fortran */

    constexpr size_t C1 = 5;
    constexpr size_t C2 = 4;
    constexpr size_t C3 = 3;
    constexpr size_t DIM1 = 3 * C1;
    constexpr size_t DIM2 = 3 * C2;
    constexpr size_t DIM3 = 3 * C3;
    const std::string filename = "ADIOSSelectionColumnMajor3D.bp";
    double a[DIM1 * DIM2 * DIM3];

    double *ap = a;
    for (size_t z = 1; z <= DIM3; ++z)
    {
        for (size_t y = 1; y <= DIM2; ++y)
        {
            for (size_t x = 1; x <= DIM1; ++x)
            {

                *ap = z * 1.0 + y / 100.0 + x / 10000.0;
                ++ap;
            }
        }
    }
    /* ./bin/bpls testing/adios2/engine/bp/bp4/ADIOSSelectionColumnMajor3D.bp
        -d a -n 15 -f "%7.4f"
    */

    // Write test data using BP
    {
        adios2::ADIOS adiosF("", "Fortran");
        adios2::IO ioWrite = adiosF.DeclareIO("TestIOWrite");
        if (!engineName.empty())
        {
            ioWrite.SetEngine(engineName);
        }
        else
        {
            // Create the BP Engine
            ioWrite.SetEngine("BPFile");
        }

        std::cout << "Write data as Column-major..." << std::endl;

        adios2::Engine engine = ioWrite.Open(filename, adios2::Mode::Write);
        const adios2::Dims shape = {DIM1, DIM2, DIM3};
        const adios2::Dims count = {C1, C2, C3};
        adios2::Dims start{0, 0, 0};
        double b[C1 * C2 * C3];

        adios2::Variable<double> var = ioWrite.DefineVariable<double>("a", shape, start, count);

        /*adios2::Variable<double> vara1 =
            ioWrite.DefineVariable<double>("a1", shape, start, shape);*/

        engine.BeginStep();
        /*engine.Put(vara1, a);*/

        for (size_t z = 0; z < DIM3; z += count[2])
        {
            for (size_t y = 0; y < DIM2; y += count[1])
            {
                for (size_t x = 0; x < DIM1; x += count[0])
                {
                    start = {x, y, z};
                    copySelection3D_CM(a, shape, start, count, b);
                    var.SetSelection({start, count});
                    engine.Put(var, b, adios2::Mode::Sync);
                }
            }
        }

        engine.EndStep();
        engine.Close();
    }

    // Read data pretending to be a Column-major language
    {
        adios2::ADIOS adiosF("", "Fortran");
        adios2::IO ioRead = adiosF.DeclareIO("TestIORead");
        ioRead.SetEngine("File");
        adios2::Engine engine = ioRead.Open(filename, adios2::Mode::Read);
        EXPECT_TRUE(engine);
        engine.BeginStep();
        adios2::Variable<double> var = ioRead.InquireVariable<double>("a");
        EXPECT_TRUE(var);
        const adios2::Dims shapeF = {DIM1, DIM2, DIM3};
        const adios2::Dims countF = {C1, C2, C3};
        adios2::Dims s{0, 0, 0};
        adios2::Dims c = shapeF;
        adios2::Dims firstNonMatch{0, 0, 0};
        std::vector<double> res;

        std::cout << "Row-major read selections with entire blocks..." << std::endl;

        /* Entire array */
        s = {0, 0, 0};
        c = shapeF;
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), DIM1 * DIM2 * DIM3);
        EXPECT_TRUE(compareSelection3D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* Single block in the center */
        s = {5, 4, 3};
        c = countF;
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* Two blocks in X direction */
        s = {5, 4, 3};
        c = {2 * countF[0], countF[1], countF[2]};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* Three blocks in Y direction */
        s = {5, 0, 3};
        c = {countF[0], 3 * countF[1], countF[2]};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* Two blocks in Z direction */
        s = {5, 4, 0};
        c = {countF[0], countF[1], 2 * countF[2]};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* Four blocks in X-Z direction */
        s = {5, 4, 3};
        c = {2 * countF[0], countF[1], 2 * countF[2]};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* 8 blocks in X-Y-Z direction */
        s = {5, 4, 3};
        c = {2 * countF[0], 2 * countF[1], 2 * countF[2]};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /*
         *   Partial blocks
         */

        std::cout << "Row-major read selections with partial blocks..." << std::endl;

        /* center part of single block in center */
        s = {6, 5, 4};
        c = {countF[0] - 2, countF[1] - 2, countF[2] - 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into two blocks in X direction */
        s = {6, 5, 4};
        c = {2 * countF[0] - 2, countF[1] - 1, countF[2] - 1};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into two blocks in Y direction */
        s = {6, 5, 4};
        c = {countF[0] - 1, 2 * countF[1] - 2, countF[2] - 1};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into three blocks in Y direction
           while contiguous in Z direction, making the read of the
           center of the three blocks as a single contiguous 4xYxZ copy */
        s = {6, 1, 3};
        c = {countF[0] - 1, 3 * countF[1] - 2, countF[2]};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into two blocks in Z direction */
        s = {6, 5, 4};
        c = {countF[0] - 1, countF[1] - 1, 2 * countF[2] - 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into four blocks in X-Y direction */
        s = {1, 5, 4};
        c = {2 * countF[0] - 2, 2 * countF[1] - 2, countF[2] - 1};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into four blocks in X-Z direction */
        s = {1, 5, 4};
        c = {2 * countF[0] - 2, countF[1] - 2, 2 * countF[2] - 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into four blocks in Y-Z direction */
        s = {6, 5, 4};
        c = {countF[0] - 2, 2 * countF[1] - 2, 2 * countF[2] - 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into six blocks in X-Z direction */
        s = {1, 5, 4};
        c = {3 * countF[0] - 2, countF[1] - 2, 2 * countF[2] - 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        /* Center block plus 1 in each direction, cutting into all blocks */
        /* ./bin/bpls -la testing/adios2/engine/bp/bp4/ADIOSSelection3D.bp
            -d a  -f "%6.3f " -s "4,3,2" -c "7,6,5" -n 5
        */
        s = {4, 3, 2};
        c = {countF[0] + 2, countF[1] + 2, countF[2] + 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_CM(a, shapeF, s, c, res.data(), firstNonMatch));

        engine.EndStep();
        engine.Close();
    }

    // Read data from Row-major language (C++)
    {
        adios2::ADIOS adios;
        adios2::IO ioRead = adios.DeclareIO("TestIOReadF");
        ioRead.SetEngine("File");
        adios2::Engine engine = ioRead.Open(filename, adios2::Mode::Read);
        EXPECT_TRUE(engine);
        engine.BeginStep();
        adios2::Variable<double> var = ioRead.InquireVariable<double>("a");
        EXPECT_TRUE(var);
        const adios2::Dims shape = {DIM3, DIM2, DIM1};
        const adios2::Dims count = {C3, C2, C1};
        adios2::Dims s{0, 0, 0};
        adios2::Dims c = shape;
        adios2::Dims firstNonMatch{0, 0, 0};
        std::vector<double> res;

        std::cout << "Column-major read selections with entire blocks..." << std::endl;

        /* Entire array */
        s = {0, 0, 0};
        c = shape;
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), DIM1 * DIM2 * DIM3);
        EXPECT_TRUE(compareSelection3D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* Single block in the center */
        s = {3, 4, 5};
        c = count;
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* Two blocks in X direction */
        s = {3, 4, 5};
        c = {2 * count[0], count[1], count[2]};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* Three blocks in Y direction */
        s = {3, 0, 5};
        c = {count[0], 3 * count[1], count[2]};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* Two blocks in Z direction */
        s = {0, 4, 5};
        c = {count[0], count[1], 2 * count[2]};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* Four blocks in X-Z direction */
        s = {3, 4, 5};
        c = {2 * count[0], count[1], 2 * count[2]};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* 8 blocks in X-Y-Z direction */
        s = {3, 4, 5};
        c = {2 * count[0], 2 * count[1], 2 * count[2]};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_RM(a, shape, s, c, res.data(), firstNonMatch));

        std::cout << "Column-major read selections with partial blocks..." << std::endl;

        /* center part of single block in center */
        s = {4, 5, 6};
        c = {count[0] - 2, count[1] - 2, count[2] - 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into two blocks in X direction */
        s = {4, 5, 6};
        c = {2 * count[0] - 2, count[1] - 1, count[2] - 1};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into two blocks in Y direction */
        s = {4, 5, 6};
        c = {count[0] - 1, 2 * count[1] - 2, count[2] - 1};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into three blocks in Y direction
           while contiguous in Z direction, making the read of the
           center of the three blocks as a single contiguous 4xYxZ copy */
        s = {3, 1, 6};
        c = {count[0] - 1, 3 * count[1] - 2, count[2]};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into two blocks in Z direction */
        s = {4, 5, 6};
        c = {count[0] - 1, count[1] - 1, 2 * count[2] - 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into four blocks in X-Y direction */
        s = {1, 5, 6};
        c = {2 * count[0] - 2, 2 * count[1] - 2, count[2] - 1};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into four blocks in X-Z direction */
        s = {1, 5, 6};
        c = {2 * count[0] - 2, count[1] - 2, 2 * count[2] - 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into four blocks in Y-Z direction */
        s = {4, 5, 6};
        c = {count[0] - 2, 2 * count[1] - 2, 2 * count[2] - 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* partial selection cutting into six blocks in X-Z direction */
        s = {1, 5, 6};
        c = {3 * count[0] - 2, count[1] - 2, 2 * count[2] - 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_RM(a, shape, s, c, res.data(), firstNonMatch));

        /* Center block plus 1 in each direction, cutting into all blocks */
        /* ./bin/bpls -la testing/adios2/engine/bp/bp4/ADIOSSelection3D.bp
            -d a  -f "%6.3f " -s "4,3,2" -c "7,6,5" -n 5
        */
        s = {2, 3, 4};
        c = {count[0] + 2, count[1] + 2, count[2] + 2};
        var.SetSelection({s, c});
        engine.Get<double>(var, res, adios2::Mode::Sync);
        EXPECT_EQ(res.size(), c[0] * c[1] * c[2]);
        EXPECT_TRUE(compareSelection3D_RM(a, shape, s, c, res.data(), firstNonMatch));

        engine.EndStep();
        engine.Close();
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc > 1)
    {
        engineName = std::string(argv[1]);
    }
    int result = RUN_ALL_TESTS();
    return result;
}

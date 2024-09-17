/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#include "operations/CudaRoutines.h"

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

class ADIOSSelectionCUDATest : public ::testing::Test
{
public:
    ADIOSSelectionCUDATest() = default;
};

void copySelection2D(const double *a, const adios2::Dims &shape, const adios2::Dims &start,
                     const adios2::Dims &count, double *b)
{
    double *bp = b;
    for (size_t x = 0; x < count[0]; ++x)
    {
        for (size_t y = 0; y < count[1]; ++y)
        {
            const size_t aidx = (start[0] + x) * shape[1] + start[1] + y;
            *bp = a[aidx];
            ++bp;
        }
    }
}

bool compareSelection2D(const double *a, const adios2::Dims &shape, const adios2::Dims &start,
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

bool compareSelection2D_F(const double *a, const adios2::Dims &shape, const adios2::Dims &start,
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

TEST_F(ADIOSSelectionCUDATest, 2D)
{
    constexpr size_t C1 = 5;
    constexpr size_t C2 = 4;
    constexpr size_t DIM1 = 3 * C1;
    constexpr size_t DIM2 = 3 * C2;
    const std::string filename = "ADIOSSelectionCuda2D.bp";
    double a[DIM1 * DIM2];

    double *ap = a;
    for (size_t x = 1; x <= DIM1; ++x)
    {
        for (size_t y = 1; y <= DIM2; ++y)
        {
            *ap = x * 1.0 + y / 100.0;
            ++ap;
        }
    }
    /* ./bin/bpls -la testing/adios2/engine/bp/bp4/ADIOSSelection2D.bp
        -d a -n 12 -f "%6.3f"
    */

    // Write test data using BP
    {
        adios2::ADIOS adios;
        adios2::IO ioWrite = adios.DeclareIO("TestIOWrite");
        if (!engineName.empty())
        {
            ioWrite.SetEngine(engineName);
        }
        else
        {
            // Create the BP Engine
            ioWrite.SetEngine("BPFile");
        }

        std::cout << "Write data as CUDA buffer..." << std::endl;

        adios2::Engine engine = ioWrite.Open(filename, adios2::Mode::Write);
        const adios2::Dims shape = {DIM1, DIM2};
        const adios2::Dims count = {C1, C2};
        adios2::Dims start{0, 0};

        double b[C1 * C2];

        adios2::Variable<double> var = ioWrite.DefineVariable<double>("a", shape, start, count);

        /*adios2::Variable<double> vara1 =
            ioWrite.DefineVariable<double>("a1", shape, start, shape);*/

        engine.BeginStep();

        double *gpuSimData;
        cudaMalloc(&gpuSimData, C1 * C2 * sizeof(double));
        for (size_t x = 0; x < DIM1; x += count[0])
        {
            for (size_t y = 0; y < DIM2; y += count[1])
            {
                start = {x, y};
                copySelection2D(a, shape, start, count, b);
                var.SetSelection({start, count});
                cudaMemcpy(gpuSimData, b, C1 * C2 * sizeof(double), cudaMemcpyHostToDevice);
                var.SetMemorySpace(adios2::MemorySpace::GPU);
                engine.Put(var, gpuSimData, adios2::Mode::Sync);
            }
        }

        engine.EndStep();
        engine.Close();
    }

    // Read data
    {
        adios2::ADIOS adios;
        adios2::IO ioRead = adios.DeclareIO("TestIORead");
        ioRead.SetEngine("File");
        adios2::Engine engine = ioRead.Open(filename, adios2::Mode::Read);
        EXPECT_TRUE(engine);
        engine.BeginStep();
        adios2::Variable<double> var = ioRead.InquireVariable<double>("a");
        EXPECT_TRUE(var);
        const adios2::Dims shape = {DIM1, DIM2};
        const adios2::Dims count = {C1, C2};
        adios2::Dims s{0, 0};
        adios2::Dims c = shape;
        adios2::Dims firstNonMatch{0, 0};

        std::cout << "CUDA read selections with entire blocks..." << std::endl;

        /* Entire array */
        {
            std::vector<double> res(DIM1 * DIM2);
            double *gpuRet;
            cudaMalloc(&gpuRet, DIM1 * DIM2 * sizeof(double));
            var.SetSelection({s, c});
            var.SetMemorySpace(adios2::MemorySpace::GPU);
            engine.Get<double>(var, gpuRet, adios2::Mode::Sync);
            cudaMemcpy(res.data(), gpuRet, DIM1 * DIM2 * sizeof(double), cudaMemcpyDeviceToHost);
            EXPECT_EQ(res.size(), DIM1 * DIM2);
            EXPECT_TRUE(compareSelection2D(a, shape, s, c, res.data(), firstNonMatch));
        }

        /* Single block in the center */
        {
            s = {5, 4};
            c = count;
            std::vector<double> res(c[0] * c[1]);
            double *gpuRet;
            cudaMalloc(&gpuRet, c[0] * c[1] * sizeof(double));
            var.SetSelection({s, c});
            var.SetMemorySpace(adios2::MemorySpace::GPU);
            engine.Get<double>(var, gpuRet, adios2::Mode::Sync);
            cudaMemcpy(res.data(), gpuRet, c[0] * c[1] * sizeof(double), cudaMemcpyDeviceToHost);
            EXPECT_EQ(res.size(), c[0] * c[1]);
            EXPECT_TRUE(compareSelection2D(a, shape, s, c, res.data(), firstNonMatch));
        }

        /* Four blocks in X-Y direction */
        {
            s = {5, 4};
            c = {2 * count[0], 2 * count[1]};
            std::vector<double> res(c[0] * c[1]);
            double *gpuRet;
            cudaMalloc(&gpuRet, c[0] * c[1] * sizeof(double));
            var.SetSelection({s, c});
            var.SetMemorySpace(adios2::MemorySpace::GPU);
            engine.Get<double>(var, gpuRet, adios2::Mode::Sync);
            cudaMemcpy(res.data(), gpuRet, c[0] * c[1] * sizeof(double), cudaMemcpyDeviceToHost);
            EXPECT_EQ(res.size(), c[0] * c[1]);
            EXPECT_TRUE(compareSelection2D(a, shape, s, c, res.data(), firstNonMatch));
        }

        /*
         *   Partial blocks
         */

        std::cout << "CUDA read selections with partial blocks..." << std::endl;

        /* center part of single block in center */
        {
            s = {6, 5};
            c = {count[0] - 2, count[1] - 2};
            std::vector<double> res(c[0] * c[1]);
            double *gpuRet;
            cudaMalloc(&gpuRet, c[0] * c[1] * sizeof(double));
            var.SetSelection({s, c});
            var.SetMemorySpace(adios2::MemorySpace::GPU);
            engine.Get<double>(var, gpuRet, adios2::Mode::Sync);
            cudaMemcpy(res.data(), gpuRet, c[0] * c[1] * sizeof(double), cudaMemcpyDeviceToHost);
            EXPECT_EQ(res.size(), c[0] * c[1]);
            EXPECT_TRUE(compareSelection2D(a, shape, s, c, res.data(), firstNonMatch));
        }

        /* Center block plus 1 in each direction, cutting into all blocks */
        /* ./bin/bpls -la testing/adios2/engine/bp/bp4/ADIOSSelection2D.bp
            -d a  -f "%6.3f " -s "4,3" -c "7,6" -n 6
        */
        {
            s = {4, 3};
            c = {count[0] + 2, count[1] + 2};
            std::vector<double> res(c[0] * c[1]);
            double *gpuRet;
            cudaMalloc(&gpuRet, c[0] * c[1] * sizeof(double));
            var.SetSelection({s, c});
            var.SetMemorySpace(adios2::MemorySpace::GPU);
            engine.Get<double>(var, gpuRet, adios2::Mode::Sync);
            cudaMemcpy(res.data(), gpuRet, c[0] * c[1] * sizeof(double), cudaMemcpyDeviceToHost);
            EXPECT_EQ(res.size(), c[0] * c[1]);
            EXPECT_TRUE(compareSelection2D(a, shape, s, c, res.data(), firstNonMatch));
        }

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

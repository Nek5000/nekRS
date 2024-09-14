/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#ifndef TESTING_ADIOS2_ENGINE_SMALLTESTDATA_H_
#define TESTING_ADIOS2_ENGINE_SMALLTESTDATA_H_

#include <cstdint>

#include <algorithm>
#include <array>
#include <limits>
#include <string>
#include <vector>

#ifdef WIN32
#define NOMINMAX
#endif

// Test data for each type.  Make sure our values exceed the range of the
// previous size to make sure we all bytes for each element
struct SmallTestData
{
    std::string S1 = "Testing ADIOS2 String type";

    // These shoudl be able to use std::array like the rest of the pieces
    // but the XL compiler seems to have some bad code generation surounding
    // it that results in a double-free corruption.  Switching to std::vector
    // bypasses the problem
    std::vector<std::string> S1array = {"one"};
    std::vector<std::string> S3 = {"one", "two", "three"};

    std::array<int8_t, 10> I8 = {{0, 1, -2, 3, -4, 5, -6, 7, -8, 9}};
    std::array<int16_t, 10> I16 = {{512, 513, -510, 515, -508, 517, -506, 519, -504, 521}};
    std::array<int32_t, 10> I32 = {
        {131072, 131073, -131070, 131075, -131068, 131077, -131066, 131079, -131064, 131081}};
    std::array<int64_t, 10> I64 = {{8589934592, 8589934593, -8589934590, 8589934595, -8589934588,
                                    8589934597, -8589934586, 8589934599, -8589934584, 8589934601}};
    std::array<uint8_t, 10> U8 = {{128, 129, 130, 131, 132, 133, 134, 135, 136, 137}};
    std::array<uint16_t, 10> U16 = {
        {32768, 32769, 32770, 32771, 32772, 32773, 32774, 32775, 32776, 32777}};
    std::array<uint32_t, 10> U32 = {{2147483648, 2147483649, 2147483650, 2147483651, 2147483652,
                                     2147483653, 2147483654, 2147483655, 2147483656, 2147483657}};
    std::array<uint64_t, 10> U64 = {
        {9223372036854775808UL, 9223372036854775809UL, 9223372036854775810UL, 9223372036854775811UL,
         9223372036854775812UL, 9223372036854775813UL, 9223372036854775814UL, 9223372036854775815UL,
         9223372036854775816UL, 9223372036854775817UL}};
    std::array<float, 10> R32 = {{0.1f, 1.1f, 2.1f, 3.1f, 4.1f, 5.1f, 6.1f, 7.1f, 8.1f, 9.1f}};
    std::array<double, 10> R64 = {{10.2, 11.2, 12.2, 13.2, 14.2, 15.2, 16.2, 17.2, 18.2, 19.2}};
    std::array<long double, 10> R128 = {
        {410.2, 411.2, 412.2, 413.2, 414.2, 415.2, 416.2, 417.2, 418.2, 419.2}};

    std::array<std::complex<float>, 10> CR32 = {
        {std::complex<float>(0.1f, 1.1f), std::complex<float>(1.1f, 2.1f),
         std::complex<float>(2.1f, 3.1f), std::complex<float>(3.1f, 4.1f),
         std::complex<float>(4.1f, 5.1f), std::complex<float>(5.1f, 6.1f),
         std::complex<float>(6.1f, 7.1f), std::complex<float>(7.1f, 8.1f),
         std::complex<float>(8.1f, 9.1f), std::complex<float>(9.1f, 10.1f)}};

    std::array<std::complex<double>, 10> CR64 = {
        {std::complex<double>(10.2, 11.2), std::complex<double>(11.2, 12.2),
         std::complex<double>(12.2, 13.2), std::complex<double>(13.2, 14.2),
         std::complex<double>(14.2, 15.2), std::complex<double>(15.2, 16.2),
         std::complex<double>(16.2, 17.2), std::complex<double>(17.2, 18.2),
         std::complex<double>(18.2, 19.2), std::complex<double>(19.2, 20.2)}};

    std::array<char, 10> CHAR = {{'a', 'b', 'c', 'y', 'z', 'A', 'B', 'C', 'Y', 'Z'}};
    std::array<bool, 10> TF = {{true, false, true, true, false, false, true, false, false, true}};

    std::array<char, 10> CHARS = {(char)0, (char)1,  (char)-2, (char)3,  (char)-4,
                                  (char)5, (char)-6, (char)7,  (char)-8, (char)9};
    std::array<signed char, 10> SCHARS = {
        (signed char)0, (signed char)1,  (signed char)-2, (signed char)3,  (signed char)-4,
        (signed char)5, (signed char)-6, (signed char)7,  (signed char)-8, (signed char)9};
    std::array<unsigned char, 10> UCHARS = {(unsigned char)0,  (unsigned char)1,  (unsigned char)-2,
                                            (unsigned char)3,  (unsigned char)-4, (unsigned char)5,
                                            (unsigned char)-6, (unsigned char)7,  (unsigned char)-8,
                                            (unsigned char)9};
};

SmallTestData generateNewSmallTestData(SmallTestData in, size_t step, size_t rank, size_t size)
{
    size_t j = rank + 1 + step * size;
    std::for_each(in.I8.begin(), in.I8.end(), [&](int8_t &v) { v += static_cast<int8_t>(j); });
    std::for_each(in.I16.begin(), in.I16.end(), [&](int16_t &v) { v += static_cast<int16_t>(j); });
    std::for_each(in.I32.begin(), in.I32.end(), [&](int32_t &v) { v += static_cast<int32_t>(j); });
    std::for_each(in.I64.begin(), in.I64.end(), [&](int64_t &v) { v += static_cast<int64_t>(j); });
    std::for_each(in.U8.begin(), in.U8.end(), [&](uint8_t &v) { v += static_cast<uint8_t>(j); });
    std::for_each(in.U16.begin(), in.U16.end(),
                  [&](uint16_t &v) { v += static_cast<uint16_t>(j); });
    std::for_each(in.U32.begin(), in.U32.end(),
                  [&](uint32_t &v) { v += static_cast<uint32_t>(j); });
    std::for_each(in.U64.begin(), in.U64.end(),
                  [&](uint64_t &v) { v += static_cast<uint64_t>(j); });
    std::for_each(in.R32.begin(), in.R32.end(), [&](float &v) { v += static_cast<float>(j); });
    std::for_each(in.R64.begin(), in.R64.end(), [&](double &v) { v += static_cast<double>(j); });
    std::for_each(in.R128.begin(), in.R128.end(), [&](long double &v) { v += j; });

    std::for_each(in.CR32.begin(), in.CR32.end(), [&](std::complex<float> &v) {
        v.real(v.real() + static_cast<float>(j));
        v.imag(v.imag() + static_cast<float>(j));
    });
    std::for_each(in.CR64.begin(), in.CR64.end(), [&](std::complex<double> &v) {
        v.real(v.real() + static_cast<double>(j));
        v.imag(v.imag() + static_cast<double>(j));
    });

    std::for_each(in.CHAR.begin(), in.CHAR.end(), [&](char &v) {
        char jc = static_cast<char>(j);
        char inc = jc % ('z' - 'a');
        v += inc;
        if (v > 'z')
        {
            v = 'a' + (v - 'z') - 1;
        }
        else if (v > 'Z' && v < 'a')
        {
            v = 'A' + (v - 'Z') - 1;
        }
    });

    std::for_each(in.CHARS.begin(), in.CHARS.end(), [&](char &v) { v += static_cast<char>(j); });
    std::for_each(in.SCHARS.begin(), in.SCHARS.end(),
                  [&](signed char &v) { v += static_cast<signed char>(j); });
    std::for_each(in.UCHARS.begin(), in.UCHARS.end(),
                  [&](unsigned char &v) { v += static_cast<unsigned char>(j); });

    return in;
}

void UpdateSmallTestData(SmallTestData &in, int step, int rank, int size)
{
    int j = rank + 1 + step * size;
    std::for_each(in.I8.begin(), in.I8.end(), [&](int8_t &v) { v += j; });
    std::for_each(in.I16.begin(), in.I16.end(), [&](int16_t &v) { v += j; });
    std::for_each(in.I32.begin(), in.I32.end(), [&](int32_t &v) { v += j; });
    std::for_each(in.I64.begin(), in.I64.end(), [&](int64_t &v) { v += j; });
    std::for_each(in.U8.begin(), in.U8.end(), [&](uint8_t &v) { v += j; });
    std::for_each(in.U16.begin(), in.U16.end(), [&](uint16_t &v) { v += j; });
    std::for_each(in.U32.begin(), in.U32.end(), [&](uint32_t &v) { v += j; });
    std::for_each(in.U64.begin(), in.U64.end(), [&](uint64_t &v) { v += j; });
    std::for_each(in.R32.begin(), in.R32.end(), [&](float &v) { v += j; });
    std::for_each(in.R64.begin(), in.R64.end(), [&](double &v) { v += j; });
    std::for_each(in.R128.begin(), in.R128.end(), [&](long double &v) { v += j; });

    std::for_each(in.CR32.begin(), in.CR32.end(), [&](std::complex<float> &v) {
        v.real(v.real() + static_cast<float>(j));
        v.imag(v.imag() + static_cast<float>(j));
    });
    std::for_each(in.CR64.begin(), in.CR64.end(), [&](std::complex<double> &v) {
        v.real(v.real() + static_cast<double>(j));
        v.imag(v.imag() + static_cast<double>(j));
    });

    std::for_each(in.CHAR.begin(), in.CHAR.end(), [&](char &v) {
        char jc = static_cast<char>(j);
        char inc = jc % ('z' - 'a');
        v += inc;
        if (v > 'z')
        {
            v = 'a' + (v - 'z') - 1;
        }
        else if (v > 'Z' && v < 'a')
        {
            v = 'A' + (v - 'Z') - 1;
        }
    });

    std::for_each(in.CHARS.begin(), in.CHARS.end(), [&](char &v) { v += static_cast<char>(j); });
    std::for_each(in.SCHARS.begin(), in.SCHARS.end(),
                  [&](signed char &v) { v += static_cast<signed char>(j); });
    std::for_each(in.UCHARS.begin(), in.UCHARS.end(),
                  [&](unsigned char &v) { v += static_cast<unsigned char>(j); });
}

#endif // TESTING_ADIOS2_ENGINE_SMALLTESTDATA_H_

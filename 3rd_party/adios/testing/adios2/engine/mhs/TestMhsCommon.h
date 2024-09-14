#include <adios2.h>
#include <gtest/gtest.h>
#include <iostream>
#include <numeric>

using namespace adios2;

int printed_lines = 0;

template <class T>
void PrintData(const T *data, const size_t step, const Dims &start, const Dims &count,
               const std::string &var, const int to_print_lines)
{
    size_t size = std::accumulate(count.begin(), count.end(), static_cast<size_t>(1),
                                  std::multiplies<size_t>());
    std::cout << " Step: " << step << " Size:" << size << "\n";
    size_t printsize = 1024;

    if (size < printsize)
    {
        printsize = size;
    }
    size_t s = 0;
    for (size_t i = 0; i < printsize; ++i)
    {
        ++s;
        std::cout << data[i] << " ";
        if (s == count[1])
        {
            std::cout << " <--- Variable " << var << ", Step " << step << std::endl;
            s = 0;
        }
    }

    std::cout << "]" << std::endl;
}

template <class T>
void GenDataRecursive(std::vector<size_t> start, std::vector<size_t> count,
                      std::vector<size_t> shape, size_t n0, size_t y,
                      std::vector<std::complex<T>> &vec, const size_t step)
{
    for (size_t i = 0; i < count[0]; i++)
    {
        size_t i0 = n0 * count[0] + i;
        size_t z = y * shape[0] + (i + start[0]);

        auto start_next = start;
        auto count_next = count;
        auto shape_next = shape;
        start_next.erase(start_next.begin());
        count_next.erase(count_next.begin());
        shape_next.erase(shape_next.begin());

        if (start_next.size() == 1)
        {
            for (size_t j = 0; j < count_next[0]; j++)
            {
                vec[i0 * count_next[0] + j] = {
                    static_cast<T>(z * shape_next[0] + (j + start_next[0]) + step), 1};
            }
        }
        else
        {
            GenDataRecursive(start_next, count_next, shape_next, i0, z, vec, step);
        }
    }
}

template <class T>
void GenDataRecursive(std::vector<size_t> start, std::vector<size_t> count,
                      std::vector<size_t> shape, size_t n0, size_t y, std::vector<T> &vec,
                      const size_t step)
{
    for (size_t i = 0; i < count[0]; i++)
    {
        size_t i0 = n0 * count[0] + i;
        size_t z = y * shape[0] + (i + start[0]);

        auto start_next = start;
        auto count_next = count;
        auto shape_next = shape;
        start_next.erase(start_next.begin());
        count_next.erase(count_next.begin());
        shape_next.erase(shape_next.begin());

        if (start_next.size() == 1)
        {
            for (size_t j = 0; j < count_next[0]; j++)
            {
                vec[i0 * count_next[0] + j] =
                    static_cast<T>(z * shape_next[0] + (j + start_next[0]) + step);
            }
        }
        else
        {
            GenDataRecursive(start_next, count_next, shape_next, i0, z, vec, step);
        }
    }
}

template <class T>
void GenData(std::vector<T> &vec, const size_t step, const std::vector<size_t> &start,
             const std::vector<size_t> &count, const std::vector<size_t> &shape)
{
    size_t total_size = std::accumulate(count.begin(), count.end(), static_cast<size_t>(1),
                                        std::multiplies<size_t>());
    vec.resize(total_size);
    GenDataRecursive(start, count, shape, 0, 0, vec, step);
}

template <class T>
void VerifyData(const std::complex<T> *data, size_t step, const Dims &start, const Dims &count,
                const Dims &shape, const std::string &var, const int to_print_lines = 0,
                const int rank = 0)
{
    size_t size = std::accumulate(count.begin(), count.end(), static_cast<size_t>(1),
                                  std::multiplies<size_t>());
    std::vector<std::complex<T>> tmpdata(size);
    GenData(tmpdata, step, start, count, shape);
    for (size_t i = 0; i < size; ++i)
    {
        ASSERT_EQ(data[i], tmpdata[i])
            << "Step " << step << " Variable " << var << " at " << i << std::endl;
    }
    if (rank == 0 && printed_lines < to_print_lines)
    {
        PrintData(data, step, start, count, var, to_print_lines);
        ++printed_lines;
    }
}

template <class T>
void VerifyData(const T *data, size_t step, const Dims &start, const Dims &count, const Dims &shape,
                const std::string &var, const int to_print_lines = 0, const int rank = 0)
{
    size_t size = std::accumulate(count.begin(), count.end(), static_cast<size_t>(1),
                                  std::multiplies<size_t>());
    std::vector<T> tmpdata(size);
    if (rank == 0 && printed_lines < to_print_lines)
    {
        PrintData(data, step, start, count, var, to_print_lines);
        ++printed_lines;
    }
    GenData(tmpdata, step, start, count, shape);
    for (size_t i = 0; i < size; ++i)
    {
        ASSERT_EQ(data[i], tmpdata[i])
            << "Step " << step << " Variable " << var << " at " << i << std::endl;
    }
}

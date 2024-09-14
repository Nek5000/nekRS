
#include <adios2.h>

#include <gtest/gtest.h>

#include <array>
#include <initializer_list>

/**
 * MultiArray
 * Simple multi-d array class for use in testing
 */
template <class T, size_t N>
class MultiArray
{
    static_assert(N > 0, "MultiArray: number of dimensions must be >= 1");

public:
    using value_type = T;
    using Index = std::array<size_t, N>;

    MultiArray(const Index &dims) : m_Dims(dims), m_Data(new T[size()]()) {}

    const T &operator[](const Index &idx) const { return m_Data[offset(idx)]; }
    T &operator[](const Index &idx) { return m_Data[offset(idx)]; }

    const T *data() const { return m_Data.get(); }
    T *data() { return m_Data.get(); }

    Index dims() const { return m_Dims; }

    size_t size() const
    {
        size_t size = m_Dims[0];
        for (size_t d = 1; d < N; d++)
        {
            size *= m_Dims[d];
        }
        return size;
    }

    size_t offset(const Index &idx) const
    {
        size_t off = idx[0];
        for (size_t d = 1; d < N; d++)
        {
            off = off * m_Dims[d] + idx[d];
        }
        return off;
    }

private:
    const Index m_Dims;
    std::unique_ptr<T[]> m_Data;
};

// should be done in general (N dims)
template <class T>
bool operator==(const MultiArray<T, 2> &a, const MultiArray<T, 2> &b)
{
    if (a.dims() != b.dims())
    {
        return false;
    }

    auto dims = a.dims();
    for (size_t i = 0; i < dims[0]; i++)
    {
        for (size_t j = 0; j < dims[1]; j++)
        {
            if (a[{i, j}] != b[{i, j}])
            {
                return false;
            }
        }
    }

    return true;
}

// should be done in general (N dims)
template <class T>
std::ostream &operator<<(std::ostream &os, const MultiArray<T, 2> &arr)
{
    auto dims = arr.dims();
    os << "MultiArray(dims = {" << dims[0] << " ," << dims[1] << "})\n";
    os << "{";
    for (size_t i = 0; i < dims[0]; i++)
    {
        if (i == 0)
        {
            os << "{";
        }
        else
        {
            os << " {";
        }
        for (size_t j = 0; j < dims[1]; j++)
        {
            os << " " << std::setw(2) << arr[{i, j}];
        }
        if (i < dims[0] - 1)
        {
            os << "},\n";
        }
        else
        {
            os << "}";
        }
    }
    os << "}";
    return os;
}

template <class T>
MultiArray<T, 2> makeArray(std::initializer_list<std::initializer_list<T>> t)
{
    auto dims = typename MultiArray<T, 2>::Index{t.size(), t.begin()->size()};
    auto arr = MultiArray<T, 2>(dims);

    T *p = arr.data();
    for (auto &row : t)
    {
        for (auto val : row)
        {
            *p++ = val;
        }
    }
    return arr;
}

std::string engine = "BPfile"; // default if no argument

TEST(MultiArray, Constructor)
{
    using MultiArrayT = MultiArray<double, 4>;

    auto arr = MultiArrayT({2, 3, 4, 5});
}

TEST(MultiArray, Access)
{
    using MultiArrayT = MultiArray<double, 4>;
    using MultiIndex = MultiArrayT::Index;

    auto dims = MultiIndex{2, 3, 4, 5};
    auto arr = MultiArrayT(dims);

    for (size_t i = 0; i < dims[0]; i++)
    {
        for (size_t j = 0; j < dims[1]; j++)
        {
            for (size_t k = 0; k < dims[2]; k++)
            {
                for (size_t m = 0; m < dims[3]; m++)
                {
                    arr[{i, j, k, m}] = i * 1000. + j * 100. + k * 10. + m;
                }
            }
        }
    }

    for (size_t i = 0; i < dims[0]; i++)
    {
        for (size_t j = 0; j < dims[1]; j++)
        {
            for (size_t k = 0; k < dims[2]; k++)
            {
                for (size_t m = 0; m < dims[3]; m++)
                {
                    EXPECT_EQ((arr[{i, j, k, m}]), i * 1000. + j * 100. + k * 10. + m);
                }
            }
        }
    }
}

TEST(MultiArray, Data)
{
    using MultiArrayT = MultiArray<double, 4>;
    using MultiIndex = MultiArrayT::Index;

    auto dims = MultiIndex{2, 3, 4, 5};
    auto arr = MultiArrayT(dims);

    EXPECT_EQ(arr.data(), (&arr[{0, 0, 0, 0}]));
}

TEST(MultiArray, Make)
{
    auto arr = makeArray({{1, 2, 3}, {4, 5, 6}});
    EXPECT_EQ((arr[{0, 0}]), 1);
    EXPECT_EQ((arr[{1, 2}]), 6);
}

/**
 * ADIOS2_CXX11_API
 */
class ADIOS2_CXX11_API : public ::testing::Test
{
public:
    ADIOS2_CXX11_API() : m_Ad() {}

    adios2::ADIOS m_Ad;
};

/**
 * ADIOS2_CXX11_API_IO
 */
class ADIOS2_CXX11_API_IO : public ADIOS2_CXX11_API
{
public:
    ADIOS2_CXX11_API_IO() : m_Io(m_Ad.DeclareIO("CXX11_API_TestIO")) {}

    adios2::IO m_Io;
};

/**
 * Selection Tests
 */

class ADIOS2_CXX11_API_Selection : public ADIOS2_CXX11_API_IO
{
public:
    using DataType = double;
    using MultiArrayT = MultiArray<DataType, 2>;
    using MultiIndex = MultiArrayT::Index;

    ADIOS2_CXX11_API_Selection()
    : m_IOWriter(m_Ad.DeclareIO("CXX11_API_Writer")), m_IOReader(m_Ad.DeclareIO("CXX11_API_Reader"))
    {
    }

    adios2::IO m_IOWriter;
    adios2::IO m_IOReader;
};

TEST_F(ADIOS2_CXX11_API_Selection, SelectionNone)
{
    // clang-format off
    auto arr = makeArray({{ 0.,  1.,  2.,  3.},
	                  {10., 11., 12., 13.},
			  {20., 21., 22., 23.}});
    auto ref = makeArray({{ 0.,  1.,  2.,  3.},
	                  {10., 11., 12., 13.},
	                  {20., 21., 22., 23.}});
    // clang-format on

    // write
    m_IOWriter.SetEngine(engine);
    auto writer = m_IOWriter.Open("test_selection_none.bp", adios2::Mode::Write);
    auto var = m_IOWriter.DefineVariable<DataType>("var", {3, 4});
    var.SetSelection({{0, 0}, {3, 4}});
    writer.Put(var, arr.data());
    writer.Close();

    // read back
    auto arr_read = MultiArrayT(ref.dims());
    auto engine = m_IOReader.Open("test_selection_none.bp", adios2::Mode::Read);
    engine.BeginStep();
    var = m_IOReader.InquireVariable<DataType>("var");
    var.SetSelection({{0, 0}, {3, 4}});
    engine.Get(var, arr_read.data());
    engine.EndStep();
    engine.Close();

    EXPECT_EQ(arr_read, ref);
}

TEST_F(ADIOS2_CXX11_API_Selection, SelectionWrite)
{
    // clang-format off
    auto arr = makeArray({{ 0.,  1.},
	                  {10., 11.},
			  {20., 21.}});
    auto ref = makeArray({{ 0.,  1., 0., 0.},
			  {10., 11., 0., 0.},
			  {20., 21., 0., 0.}});
    // clang-format on

    // write
    m_IOWriter.SetEngine(engine);
    auto writer = m_IOWriter.Open("test_selection_write.bp", adios2::Mode::Write);
    auto var = m_IOWriter.DefineVariable<DataType>("var", {3, 4});
    var.SetSelection({{0, 0}, {3, 2}});
    writer.Put(var, arr.data());
    writer.Close();

    // read back
    auto arr_read = MultiArrayT(ref.dims());
    auto engine = m_IOReader.Open("test_selection_write.bp", adios2::Mode::Read);
    engine.BeginStep();
    var = m_IOReader.InquireVariable<DataType>("var");
    var.SetSelection({{0, 0}, {3, 4}});
    engine.Get(var, arr_read.data());
    engine.EndStep();
    engine.Close();

    EXPECT_EQ(arr_read, ref);
}

TEST_F(ADIOS2_CXX11_API_Selection, SelectionWriteStart)
{
    // clang-format off
    auto arr = makeArray({{ 0.,  1.},
	 		  {10., 11.},
			  {20., 21.}});
    auto ref = makeArray({{0., 0.,  0.,  1.},
			  {0., 0., 10., 11.},
			  {0., 0., 20., 21.}});
    // clang-format on

    // write
    m_IOWriter.SetEngine(engine);
    auto writer = m_IOWriter.Open("test_selection_write_start.bp", adios2::Mode::Write);
    auto var = m_IOWriter.DefineVariable<DataType>("var", {3, 4});
    var.SetSelection({{0, 2}, {3, 2}});
    writer.Put(var, arr.data());
    writer.Close();

    // read back
    auto arr_read = MultiArrayT(ref.dims());
    auto engine = m_IOReader.Open("test_selection_write_start.bp", adios2::Mode::Read);
    engine.BeginStep();
    var = m_IOReader.InquireVariable<DataType>("var");
    var.SetSelection({{0, 0}, {3, 4}});
    engine.Get(var, arr_read.data());
    engine.EndStep();
    engine.Close();

    EXPECT_EQ(arr_read, ref);
}

TEST_F(ADIOS2_CXX11_API_Selection, SelectionRead)
{
    // clang-format off
    auto arr = makeArray({{ 0.,  1.,  2.,  3.},
			  {10., 11., 12., 13.},
			  {20., 21., 22., 23.}});
    auto ref = makeArray({{ 0.,  1.},
			  {10., 11.},
			  {20., 21.}});
    // clang-format on

    m_IOWriter.SetEngine(engine);
    auto writer = m_IOWriter.Open("test_selection_read.bp", adios2::Mode::Write);
    auto var = m_IOWriter.DefineVariable<DataType>("var", {3, 4});
    var.SetSelection({{0, 0}, {3, 4}});
    writer.Put(var, arr.data());
    writer.Close();

    // read back
    auto arr_read = MultiArrayT(ref.dims());
    auto engine = m_IOReader.Open("test_selection_read.bp", adios2::Mode::Read);
    engine.BeginStep();
    var = m_IOReader.InquireVariable<DataType>("var");
    var.SetSelection({{0, 0}, {3, 2}});
    engine.Get(var, arr_read.data());
    engine.EndStep();
    engine.Close();

    EXPECT_EQ(arr_read, ref);
}

TEST_F(ADIOS2_CXX11_API_Selection, SelectionReadStart)
{
    // clang-format off
    auto arr = makeArray({{ 0.,  1.,  2.,  3.},
			  {10., 11., 12., 13.},
			  {20., 21., 22., 23.}});
    auto ref = makeArray({{ 2.,  3.},
			  {12., 13.},
			  {22., 23.}});
    // clang-format on

    // write
    m_IOWriter.SetEngine(engine);
    auto writer = m_IOWriter.Open("test_selection_read_start.bp", adios2::Mode::Write);
    auto var = m_IOWriter.DefineVariable<DataType>("var", {3, 4});
    var.SetSelection({{0, 0}, {3, 4}});
    writer.Put(var, arr.data());
    writer.Close();

    // read back
    auto arr_read = MultiArrayT(ref.dims());
    auto engine = m_IOReader.Open("test_selection_read_start.bp", adios2::Mode::Read);
    engine.BeginStep();
    var = m_IOReader.InquireVariable<DataType>("var");
    var.SetSelection({{0, 2}, {3, 2}});
    engine.Get(var, arr_read.data());
    engine.EndStep();
    engine.Close();

    EXPECT_EQ(arr_read, ref);
}

TEST_F(ADIOS2_CXX11_API_Selection, MemorySelectionNone)
{
    // clang-format off
    auto arr = makeArray({{0.,  0.,  0.,  0.,  0., 0.},
			  {0.,  0.,  1.,  2.,  3., 0.},
			  {0., 10., 11., 12., 13., 0.},
			  {0., 20., 21., 22., 23., 0.},
			  {0.,  0.,  0.,  0.,  0., 0.}});
    auto ref = makeArray({{ 0.,  1.,  2.,  3.},
			  {10., 11., 12., 13.},
			  {20., 21., 22., 23.}});
    // clang-format on

    // write
    m_IOWriter.SetEngine(engine);
    auto writer = m_IOWriter.Open("test_mem_selection_none.bp", adios2::Mode::Write);
    auto var = m_IOWriter.DefineVariable<DataType>("var", {3, 4});
    var.SetSelection({{0, 0}, {3, 4}});
    var.SetMemorySelection({{1, 1}, {5, 6}});
    writer.Put(var, arr.data());
    writer.Close();

    // read back
    auto arr_read = MultiArrayT(ref.dims());
    auto engine = m_IOReader.Open("test_mem_selection_none.bp", adios2::Mode::Read);
    engine.BeginStep();
    var = m_IOReader.InquireVariable<DataType>("var");
    var.SetSelection({{0, 0}, {3, 4}});
    engine.Get(var, arr_read.data());
    engine.EndStep();
    engine.Close();

    EXPECT_EQ(arr_read, ref);
}

TEST_F(ADIOS2_CXX11_API_Selection, MemorySelectionWrite)
{
    // clang-format off
    auto arr = makeArray({{0.,  0.,  0.,  0.,  0., 0.},
			  {0.,  0.,  1.,  2.,  3., 0.},
			  {0., 10., 11., 12., 13., 0.},
			  {0., 20., 21., 22., 23., 0.},
			  {0.,  0.,  0.,  0.,  0., 0.}});
    auto ref = makeArray({{ 0.,  1.},
			  {10., 11.},
			  {20., 21.}});
    // clang-format on

    // write
    m_IOWriter.SetEngine(engine);
    auto writer = m_IOWriter.Open("test_mem_selection_write.bp", adios2::Mode::Write);
    auto var = m_IOWriter.DefineVariable<DataType>("var", {3, 4});
    var.SetSelection({{0, 0}, {3, 2}});
    var.SetMemorySelection({{1, 1}, {5, 6}});
    writer.Put(var, arr.data());
    writer.Close();

    // read back
    auto arr_read = MultiArrayT(ref.dims());
    auto engine = m_IOReader.Open("test_mem_selection_write.bp", adios2::Mode::Read);
    engine.BeginStep();
    var = m_IOReader.InquireVariable<DataType>("var");
    var.SetSelection({{0, 0}, {3, 2}});
    engine.Get(var, arr_read.data());
    engine.EndStep();
    engine.Close();

    EXPECT_EQ(arr_read, ref);
}

TEST_F(ADIOS2_CXX11_API_Selection, MemorySelectionWriteStart)
{
    // clang-format off
    auto arr = makeArray({{0.,  0.,  0.,  0.,  0., 0.},
			  {0.,  0.,  1.,  2.,  3., 0.},
			  {0., 10., 11., 12., 13., 0.},
			  {0., 20., 21., 22., 23., 0.},
			  {0.,  0.,  0.,  0.,  0., 0.}});
    auto ref = makeArray({{ 2.,  3.},
			  {12., 13.},
			  {22., 23.}});
    // clang-format on

    // write
    m_IOWriter.SetEngine(engine);
    auto writer = m_IOWriter.Open("test_mem_selection_write_start.bp", adios2::Mode::Write);
    auto var = m_IOWriter.DefineVariable<DataType>("var", {3, 4});
    var.SetSelection({{0, 0}, {3, 2}});
    var.SetMemorySelection({{1, 3}, {5, 6}});
    writer.Put(var, arr.data());
    writer.Close();

    // read back
    auto arr_read = MultiArrayT(ref.dims());
    auto engine = m_IOReader.Open("test_mem_selection_write_start.bp", adios2::Mode::Read);
    engine.BeginStep();
    var = m_IOReader.InquireVariable<DataType>("var");
    var.SetSelection({{0, 0}, {3, 2}});
    engine.Get(var, arr_read.data());
    engine.EndStep();
    engine.Close();

    EXPECT_EQ(arr_read, ref);
}

TEST_F(ADIOS2_CXX11_API_Selection, MemorySelectionRead)
{
    // clang-format off
    auto arr = makeArray({{ 0.,  1.,  2.,  3.},
			  {10., 11., 12., 13.},
			  {20., 21., 22., 23.}});
    auto ref = makeArray({{0.,  0.,  0., 0.},
			  {0.,  0.,  1., 0.},
			  {0., 10., 11., 0.},
			  {0., 20., 21., 0.},
			  {0.,  0.,  0., 0.}});
    // clang-format on

    m_IOWriter.SetEngine(engine);
    auto writer = m_IOWriter.Open("test_mem_selection_read.bp", adios2::Mode::Write);
    auto var = m_IOWriter.DefineVariable<DataType>("var", {3, 4});
    var.SetSelection({{0, 0}, {3, 4}});
    writer.Put(var, arr.data());
    writer.Close();

    // read back
    auto arr_read = MultiArrayT(ref.dims());
    auto engine = m_IOReader.Open("test_mem_selection_read.bp", adios2::Mode::Read);
    engine.BeginStep();
    var = m_IOReader.InquireVariable<DataType>("var");
    var.SetSelection({{0, 0}, {3, 2}});
    var.SetMemorySelection({{1, 1}, {5, 4}});
    engine.Get(var, arr_read.data());
    engine.EndStep();
    engine.Close();

    EXPECT_EQ(arr_read, ref);
}

TEST_F(ADIOS2_CXX11_API_Selection, MemorySelectionReadStart)
{
    // clang-format off
    auto arr = makeArray({{ 0.,  1.,  2.,  3.},
			  {10., 11., 12., 13.},
			  {20., 21., 22., 23.}});
    auto ref = makeArray({{0.,  0.,  0., 0.},
			  {0.,  2.,  3., 0.},
			  {0., 12., 13., 0.},
			  {0., 22., 23., 0.},
			  {0.,  0.,  0., 0.}});
    // clang-format on

    m_IOWriter.SetEngine(engine);
    auto writer = m_IOWriter.Open("test_mem_selection_read_start.bp", adios2::Mode::Write);
    auto var = m_IOWriter.DefineVariable<DataType>("var", {3, 4});
    var.SetSelection({{0, 0}, {3, 4}});
    writer.Put(var, arr.data());
    writer.Close();

    // read back
    auto arr_read = MultiArrayT(ref.dims());
    auto engine =
        m_IOReader.Open("test_mem_selection_read_start.bp", adios2::Mode::ReadRandomAccess);
    var = m_IOReader.InquireVariable<DataType>("var");
    var.SetSelection({{0, 2}, {3, 2}});
    var.SetMemorySelection({{1, 1}, {5, 4}});
    engine.Get(var, arr_read.data());
    engine.Close();

    EXPECT_EQ(arr_read, ref);
}

TEST_F(ADIOS2_CXX11_API_Selection, MemorySelectionComplex)
{
    // clang-format off
    auto arr = makeArray({{0.,  0.,  0.,  0.,  0., 0.},
			  {0.,  0.,  1.,  2.,  3., 0.},
			  {0., 10., 11., 12., 13., 0.},
			  {0., 20., 21., 22., 23., 0.},
			  {0., 30., 31., 32., 33., 0.},
			  {0.,  0.,  0.,  0.,  0., 0.}});
    auto ref = makeArray({{0.,  0.,  0., 0.},
			  {0., 11., 12., 0.},
			  {0., 21., 22., 0.},
			  {0.,  0.,  0., 0.}});
    // clang-format on

    m_IOWriter.SetEngine(engine);
    auto writer = m_IOWriter.Open("test_mem_selection_complex.bp", adios2::Mode::Write);
    auto var = m_IOWriter.DefineVariable<DataType>("var", {4, 4});
    // write in 4 quarters
    var.SetSelection({{0, 0}, {2, 2}});
    var.SetMemorySelection({{1, 1}, {6, 6}});
    writer.Put(var, arr.data());
    var.SetSelection({{2, 0}, {2, 2}});
    var.SetMemorySelection({{3, 1}, {6, 6}});
    writer.Put(var, arr.data());
    var.SetSelection({{0, 2}, {2, 2}});
    var.SetMemorySelection({{1, 3}, {6, 6}});
    writer.Put(var, arr.data());
    var.SetSelection({{2, 2}, {2, 2}});
    var.SetMemorySelection({{3, 3}, {6, 6}});
    writer.Put(var, arr.data());
    writer.Close();

    // read back center block, with bits from every block written
    auto arr_read = MultiArrayT(ref.dims());
    auto engine = m_IOReader.Open("test_mem_selection_complex.bp", adios2::Mode::ReadRandomAccess);
    var = m_IOReader.InquireVariable<DataType>("var");
    var.SetSelection({{1, 1}, {2, 2}});
    var.SetMemorySelection({{1, 1}, {4, 4}});
    engine.Get(var, arr_read.data());
    engine.Close();

    EXPECT_EQ(arr_read, ref);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc > 1)
    {
        engine = argv[1];
        std::cout << "Running with engine " << engine << std::endl;
    }
    return RUN_ALL_TESTS();
}

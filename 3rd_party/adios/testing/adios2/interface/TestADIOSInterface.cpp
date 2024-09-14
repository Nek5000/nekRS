#include <cstdint>

#include <iostream>
#include <numeric>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

#if ADIOS2_USE_MPI

TEST(ADIOSInterface, MPICommRemoved)
{
    MPI_Comm myComm;
    MPI_Comm_dup(MPI_COMM_WORLD, &myComm);
    adios2::ADIOS adios(myComm);
    adios2::IO io = adios.DeclareIO("TestIO");
    MPI_Comm_free(&myComm);

    adios2::Engine engine = io.Open("test.bp", adios2::Mode::Write);
    (void)engine;
}

#endif

TEST(ADIOSInterface, BADConfigFile)
{
    EXPECT_THROW(adios2::ADIOS adios("notthere.xml"); adios2::IO io = adios.DeclareIO("TestIO");

                 io.Open("test.bp", adios2::Mode::Write);, std::logic_error);
}

/** ADIOS2_CXX11_API
 */
class ADIOS2_CXX11_API : public ::testing::Test
{
public:
    ADIOS2_CXX11_API()
#if ADIOS2_USE_MPI
    : m_Ad(MPI_COMM_WORLD)
#else
    : m_Ad()
#endif
    {
#if ADIOS2_USE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &m_MpiRank);
        MPI_Comm_size(MPI_COMM_WORLD, &m_MpiSize);
#endif
    }

    adios2::ADIOS m_Ad;
    int m_MpiRank = 0;
    int m_MpiSize = 1;
};

TEST_F(ADIOS2_CXX11_API, ToString)
{
    EXPECT_EQ(ToString(adios2::ShapeID::Unknown), "ShapeID::Unknown");
    EXPECT_EQ(ToString(adios2::ShapeID::GlobalValue), "ShapeID::GlobalValue");
    EXPECT_EQ(ToString(adios2::ShapeID::GlobalArray), "ShapeID::GlobalArray");
    EXPECT_EQ(ToString(adios2::ShapeID::JoinedArray), "ShapeID::JoinedArray");
    EXPECT_EQ(ToString(adios2::ShapeID::LocalValue), "ShapeID::LocalValue");
    EXPECT_EQ(ToString(adios2::ShapeID::LocalArray), "ShapeID::LocalArray");

    EXPECT_EQ(ToString(adios2::IOMode::Independent), "IOMode::Independent");
    EXPECT_EQ(ToString(adios2::IOMode::Collective), "IOMode::Collective");

    EXPECT_EQ(ToString(adios2::Mode::Undefined), "Mode::Undefined");
    EXPECT_EQ(ToString(adios2::Mode::Write), "Mode::Write");
    EXPECT_EQ(ToString(adios2::Mode::Read), "Mode::Read");
    EXPECT_EQ(ToString(adios2::Mode::Append), "Mode::Append");
    EXPECT_EQ(ToString(adios2::Mode::Sync), "Mode::Sync");
    EXPECT_EQ(ToString(adios2::Mode::Deferred), "Mode::Deferred");

    EXPECT_EQ(ToString(adios2::ReadMultiplexPattern::GlobalReaders),
              "ReadMultiplexPattern::GlobalReaders");
    EXPECT_EQ(ToString(adios2::ReadMultiplexPattern::RoundRobin),
              "ReadMultiplexPattern::RoundRobin");
    EXPECT_EQ(ToString(adios2::ReadMultiplexPattern::FirstInFirstOut),
              "ReadMultiplexPattern::FirstInFirstOut");
    EXPECT_EQ(ToString(adios2::ReadMultiplexPattern::OpenAllSteps),
              "ReadMultiplexPattern::OpenAllSteps");

    EXPECT_EQ(ToString(adios2::StreamOpenMode::Wait), "StreamOpenMode::Wait");
    EXPECT_EQ(ToString(adios2::StreamOpenMode::NoWait), "StreamOpenMode::NoWait");

    EXPECT_EQ(ToString(adios2::ReadMode::NonBlocking), "ReadMode::NonBlocking");
    EXPECT_EQ(ToString(adios2::ReadMode::Blocking), "ReadMode::Blocking");

    EXPECT_EQ(ToString(adios2::StepMode::Append), "StepMode::Append");
    EXPECT_EQ(ToString(adios2::StepMode::Update), "StepMode::Update");
    EXPECT_EQ(ToString(adios2::StepMode::Read), "StepMode::Read");

    EXPECT_EQ(ToString(adios2::StepStatus::OK), "StepStatus::OK");
    EXPECT_EQ(ToString(adios2::StepStatus::NotReady), "StepStatus::NotReady");
    EXPECT_EQ(ToString(adios2::StepStatus::EndOfStream), "StepStatus::EndOfStream");
    EXPECT_EQ(ToString(adios2::StepStatus::OtherError), "StepStatus::OtherError");

    EXPECT_EQ(ToString(adios2::TimeUnit::Microseconds), "TimeUnit::Microseconds");
    EXPECT_EQ(ToString(adios2::TimeUnit::Milliseconds), "TimeUnit::Milliseconds");
    EXPECT_EQ(ToString(adios2::TimeUnit::Seconds), "TimeUnit::Seconds");
    EXPECT_EQ(ToString(adios2::TimeUnit::Minutes), "TimeUnit::Minutes");
    EXPECT_EQ(ToString(adios2::TimeUnit::Hours), "TimeUnit::Hours");

    EXPECT_EQ(ToString(adios2::SelectionType::BoundingBox), "SelectionType::BoundingBox");
    EXPECT_EQ(ToString(adios2::SelectionType::WriteBlock), "SelectionType::WriteBlock");
}

TEST_F(ADIOS2_CXX11_API, APIToString)
{
    auto io = m_Ad.DeclareIO("CXX11_API_TestIO");
    EXPECT_EQ(ToString(io), "IO(Name: \"CXX11_API_TestIO\")");

    auto variable = io.DefineVariable<double>("var_double");
    EXPECT_EQ(ToString(variable), "Variable<double>(Name: \"var_double\")");

    auto attribute = io.DefineAttribute<float>("attr_float", 1.f);
    EXPECT_EQ(ToString(attribute), "Attribute<float>(Name: \"attr_float\")");

    auto engine = io.Open("test.bp", adios2::Mode::Write);

    EXPECT_EQ(ToString(engine), "Engine(Name: \"test.bp\", Type: \"BP5Writer\")");
}

TEST_F(ADIOS2_CXX11_API, operatorLL)
{
    std::stringstream result;
    result << adios2::Mode::Write;
    EXPECT_EQ(result.str(), "Mode::Write");
}

/** ADIOS2_CXX11_API_IO
 */
class ADIOS2_CXX11_API_IO : public ADIOS2_CXX11_API
{
public:
    ADIOS2_CXX11_API_IO() : m_Io(m_Ad.DeclareIO("CXX11_API_TestIO")) {}

    adios2::IO m_Io;
};

TEST_F(ADIOS2_CXX11_API_IO, Engine)
{
    m_Io.SetEngine("file");
    EXPECT_EQ(m_Io.EngineType(), "file");

    adios2::Engine engine = m_Io.Open("types.bp", adios2::Mode::Write);
    EXPECT_EQ(engine.Name(), "types.bp");
    EXPECT_EQ(engine.Type(), "BP5Writer");

    engine.Close();
}

TEST_F(ADIOS2_CXX11_API_IO, EngineDefault)
{
    m_Io.SetEngine("");
    EXPECT_EQ(m_Io.EngineType(), "");

    adios2::Engine engine = m_Io.Open("types.bp", adios2::Mode::Write);
    EXPECT_EQ(engine.Name(), "types.bp");
    EXPECT_EQ(engine.Type(), "BP5Writer");
    engine.Close();
}

template <class T>
struct MyData
{
    using Block = std::vector<T>;
    using Box = adios2::Box<adios2::Dims>;

    explicit MyData(const std::vector<Box> &selections)
    : m_Blocks(selections.size()), m_Selections(selections)
    {
        for (int b = 0; b < static_cast<int>(NBlocks()); ++b)
        {
            m_Blocks[b].resize(selections[b].second[0]);
        }
    }

    size_t NBlocks() const { return m_Selections.size(); }
    size_t Start(int b) const { return m_Selections[b].first[0]; }
    size_t Count(int b) const { return m_Selections[b].second[0]; }
    const Box &Selection(int b) const { return m_Selections[b]; }
    Block &operator[](int b) { return m_Blocks[b]; }
    T *Begin(int b) { return &(*m_Blocks[b].begin()); }
    T *End(int b) { return Begin(b) + Count(b); }

private:
    std::vector<Block> m_Blocks;
    std::vector<Box> m_Selections;
};

template <class T>
struct MyDataView
{
    using Block = T *;
    using Box = adios2::Box<adios2::Dims>;

    explicit MyDataView(const std::vector<Box> &selections)
    : m_Blocks(selections.size()), m_Selections(selections)
    {
    }

    void Place(int b, T *arr) { m_Blocks[b] = arr; }

    size_t NBlocks() const { return m_Selections.size(); }
    size_t Start(int b) const { return m_Selections[b].first[0]; }
    size_t Count(int b) const { return m_Selections[b].second[0]; }
    const Box &Selection(int b) const { return m_Selections[b]; }
    Block operator[](int b) { return m_Blocks[b]; }
    T *Begin(int b) { return m_Blocks[b]; }
    T *End(int b) { return m_Blocks[b] + Count(b); }

private:
    std::vector<Block> m_Blocks;
    std::vector<Box> m_Selections;
};

// ----------------------------------------------------------------------
// ADIOS2_CXX11_API_MultiBlock

template <typename _DataType, adios2::Mode _PutMode>
struct CaseMultiBlock
{
    using DataType = _DataType;
    static const adios2::Mode PutMode = _PutMode;
};

using MultiBlockTypes = ::testing::Types<CaseMultiBlock<double, adios2::Mode::Sync>,
                                         CaseMultiBlock<int64_t, adios2::Mode::Sync>,
                                         CaseMultiBlock<double, adios2::Mode::Deferred>,
                                         CaseMultiBlock<int64_t, adios2::Mode::Deferred>>;

template <typename TypeParam>
class ADIOS2_CXX11_API_MultiBlock : public ADIOS2_CXX11_API_IO
{
public:
    using DataType = typename TypeParam::DataType;
    using Box = adios2::Box<adios2::Dims>;

    using ADIOS2_CXX11_API_IO::ADIOS2_CXX11_API_IO;

    void SetUp() override
    {
        m_Shape = {m_MpiSize * m_Nx};
        m_Selections = {{{m_MpiRank * m_Nx}, {m_Nx / 2}},
                        {{m_MpiRank * m_Nx + m_Nx / 2}, {m_Nx / 2}}};
    }

    template <class MyData>
    void PopulateBlock(MyData &myData, int b)
    {
        std::iota(myData.Begin(b), myData.End(b), DataType(myData.Start(b)));
    }

    void GenerateOutput(std::string filename, std::string enginename)
    {
        auto io = m_Ad.DeclareIO("CXX11_API_GenerateOutput");
        io.SetEngine(enginename);
        auto engine = io.Open(filename, adios2::Mode::Write);
        auto var = io.template DefineVariable<DataType>("var", m_Shape);

        MyData<DataType> myData(m_Selections);

        for (int b = 0; b < static_cast<int>(myData.NBlocks()); ++b)
        {
            PopulateBlock(myData, b);

            var.SetSelection(myData.Selection(b));
            engine.Put(var, &myData[b][0], adios2::Mode::Sync);
        }
        engine.Close();

        m_Ad.RemoveIO("CXX11_API_GenerateOutput");
    }

    void CheckOutput(std::string filename)
    {
        if (m_MpiRank != 0)
        {
            return;
        }
        auto io = m_Ad.DeclareIO("CXX11_API_CheckOutput");
#if ADIOS2_USE_MPI
        auto engine = io.Open(filename, adios2::Mode::Read, MPI_COMM_SELF);
#else
        auto engine = io.Open(filename, adios2::Mode::Read);
#endif
        engine.BeginStep();
        auto var = io.template InquireVariable<DataType>("var");
        auto shape = var.Shape();
        std::vector<DataType> data(shape[0]);
        engine.Get(var, data, adios2::Mode::Sync);
        engine.EndStep();
        engine.Close();
        m_Ad.RemoveIO("CXX11_API_CheckOutput");

        std::vector<DataType> ref(shape[0]);
        std::iota(ref.begin(), ref.end(), DataType());
        EXPECT_EQ(data, ref);
    }

    size_t m_Nx = 10;
    adios2::Dims m_Shape;
    std::vector<Box> m_Selections;
};

TYPED_TEST_SUITE(ADIOS2_CXX11_API_MultiBlock, MultiBlockTypes);

TYPED_TEST(ADIOS2_CXX11_API_MultiBlock, Put)
{
    using T = typename TypeParam::DataType;

    std::string filename = "multi_put.bp";
    auto writer = this->m_Io.Open(filename, adios2::Mode::Write);
    auto var = this->m_Io.template DefineVariable<T>("var", this->m_Shape);
    MyData<T> myData(this->m_Selections);

    for (int b = 0; b < static_cast<int>(myData.NBlocks()); ++b)
    {
        this->PopulateBlock(myData, b);
        var.SetSelection(myData.Selection(b));
        writer.Put(var, &myData[b][0], TypeParam::PutMode);
    }

    writer.Close();
    this->CheckOutput(filename);
}

TYPED_TEST(ADIOS2_CXX11_API_MultiBlock, PutMixed)
{
    using T = typename TestFixture::DataType;

    std::string filename = "multi_putmixed.bp";
    this->m_Io.SetEngine("BP3");
    auto writer = this->m_Io.Open(filename, adios2::Mode::Write);
    auto var = this->m_Io.template DefineVariable<T>("var", this->m_Shape);
    MyData<T> myData(this->m_Selections);

    for (int b = 0; b < static_cast<int>(myData.NBlocks()); ++b)
    {
        this->PopulateBlock(myData, b);

        var.SetSelection(myData.Selection(b));
        if (b == 0)
        {
            writer.Put(var, &myData[b][0], adios2::Mode::Deferred);
        }
        else
        {
            writer.Put(var, &myData[b][0], adios2::Mode::Sync);
        }
    }

    writer.Close();
    this->CheckOutput(filename);
}

// this one doesn't actually use the PutMode param
TYPED_TEST(ADIOS2_CXX11_API_MultiBlock, PutZeroCopy)
{
    using T = typename TestFixture::DataType;

    std::string filename = "multi_putzerocopy.bp";
    this->m_Io.SetEngine("BP3");
    auto writer = this->m_Io.Open(filename, adios2::Mode::Write);
    auto var = this->m_Io.template DefineVariable<T>("var", this->m_Shape);
    MyDataView<T> myData(this->m_Selections);

#if 0
    // keeping this as an example on how not to use Span
    for (int b = 0; b < myData.NBlocks(); ++b)
    {
        var.SetSelection(myData.Selection(b));
        auto span = writer.Put(var);
        myData.Place(b, span.data());
    }
#else
    std::vector<typename adios2::Variable<T>::Span> spans;
    for (int b = 0; b < static_cast<int>(myData.NBlocks()); ++b)
    {
        var.SetSelection(myData.Selection(b));
        spans.push_back(writer.Put(var));
    }
    for (int b = 0; b < static_cast<int>(myData.NBlocks()); ++b)
    {
        myData.Place(b, spans[b].data());
    }
#endif
    for (int b = 0; b < static_cast<int>(myData.NBlocks()); ++b)
    {
        this->PopulateBlock(myData, b);
    }

    writer.Close();
    this->CheckOutput(filename);
}

// writes last block using regular Put(Sync/Deferred)
TYPED_TEST(ADIOS2_CXX11_API_MultiBlock, PutZeroCopyMixed)
{
    using T = typename TestFixture::DataType;

    std::string filename = "multi_putzerocopymixed.bp";
    this->m_Io.SetEngine("BP3");
    auto writer = this->m_Io.Open(filename, adios2::Mode::Write);
    auto var = this->m_Io.template DefineVariable<T>("var", this->m_Shape);

    MyDataView<T> myData(this->m_Selections);
    std::unique_ptr<T[]> lastBlock;
    for (int b = 0; b < static_cast<int>(myData.NBlocks()); ++b)
    {
        if (b < 1)
        {
            var.SetSelection(myData.Selection(b));
            auto span = writer.Put(var);
            myData.Place(b, span.data());
        }
        else
        {
            lastBlock.reset(new T[myData.Count(b)]);
            myData.Place(b, lastBlock.get());
        }
    }

    for (int b = 0; b < static_cast<int>(myData.NBlocks()); ++b)
    {
        this->PopulateBlock(myData, b);
    }
    var.SetSelection(myData.Selection(1));
    writer.Put(var, lastBlock.get(), TypeParam::PutMode);

    writer.Close();
    this->CheckOutput(filename);
}

TYPED_TEST(ADIOS2_CXX11_API_MultiBlock, Put2File)
{
    using T = typename TypeParam::DataType;

    this->GenerateOutput("multi_2f_input.bp", "BP3");

    std::string filename = "multi_put2file.bp";
    this->m_Io.SetEngine("BP3");
    auto writer = this->m_Io.Open(filename, adios2::Mode::Write);
    auto reader = this->m_Io.Open("multi_2f_input.bp", adios2::Mode::Read);
    auto var = this->m_Io.template InquireVariable<T>("var");

    MyData<T> myData(this->m_Selections);

    for (int b = 0; b < static_cast<int>(myData.NBlocks()); ++b)
    {
        var.SetSelection(myData.Selection(b));
        reader.Get(var, &myData[b][0], adios2::Mode::Sync);
        writer.Put(var, &myData[b][0], TypeParam::PutMode);
    }

    // Close the writer before the reader because the var goes away when the
    // reader does
    writer.Close();
    reader.Close();
    this->CheckOutput(filename);
}

// write two files simultaneously (in a more realistic use case, one
// would write different data into the two files, but this is enough to
// show a problem)
#if 0
TYPED_TEST(ADIOS2_CXX11_API_MultiBlock, Put2Writers)
{
    using T = typename TypeParam::DataType;

    std::string filename = "multi_put2writers.bp";
    auto writer = this->m_Io.Open(filename, adios2::Mode::Write);
    auto writer2 =
        this->m_Io.Open("multi_put2writers2.bp", adios2::Mode::Write);
    auto var = this->m_Io.template DefineVariable<T>("var", this->m_Shape);

    MyData<T> myData(this->m_Selections);

    for (int b = 0; b < myData.NBlocks(); ++b)
    {
        this->PopulateBlock(myData, b);
        var.SetSelection(myData.Selection(b));
        writer.Put(var, &myData[b][0], TypeParam::PutMode);
        writer2.Put(var, &myData[b][0], TypeParam::PutMode);
    }
    writer2.Close();
    writer.Close();

    this->CheckOutput(filename);
    this->CheckOutput("multi_put2writers2.bp");
}
#endif

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

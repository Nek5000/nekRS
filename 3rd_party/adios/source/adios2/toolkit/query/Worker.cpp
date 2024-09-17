#include "Worker.h"
#include "adios2/helper/adiosFunctions.h"

namespace adios2
{
namespace query
{
bool EndsWith(const std::string &hostStr, const std::string &fileTag)
{
    if (hostStr.size() >= fileTag.size() &&
        hostStr.compare(hostStr.size() - fileTag.size(), fileTag.size(), fileTag) == 0)
        return true;
    else
        return false;
}

bool IsFileNameXML(const std::string &filename) { return EndsWith(filename, ".xml"); }

bool IsFileNameJSON(const std::string &filename) { return EndsWith(filename, ".json"); }

Worker::Worker(const std::string &queryFile, adios2::core::Engine *adiosEngine)
: m_QueryFile(queryFile), m_SourceReader(adiosEngine)
{
}

Worker::~Worker()
{
    if (m_Query != nullptr)
        delete m_Query;
}

Worker *GetWorker(const std::string &configFile, adios2::core::Engine *adiosEngine)
{
    std::ifstream fileStream(configFile);

    if (!fileStream)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "query::Worker", "GetWorker",
                                              "file " + configFile + " not found");
    }

    if (adios2::query::IsFileNameXML(configFile))
    {
        return new XmlWorker(configFile, adiosEngine);
    }

#ifdef ADIOS2_HAVE_DATAMAN // so json is included
    if (adios2::query::IsFileNameJSON(configFile))
    {
        return new JsonWorker(configFile, adiosEngine);
    }
#endif
    helper::Throw<std::invalid_argument>("Toolkit", "query::Worker", "GetWorker",
                                         "Unable to construct xml query");
    return nullptr;
}

QueryVar *Worker::GetBasicVarQuery(adios2::core::IO &currentIO, const std::string &variableName)
{
    const DataType varType = currentIO.InquireVariableType(variableName);
    if (varType == DataType::None)
    {
        helper::Log("Query", "Worker", "GetBasicVarQuery", "No such variable: " + variableName,
                    helper::FATALERROR);
        return nullptr;
    }
#define declare_type(T)                                                                            \
    if (varType == helper::GetDataType<T>())                                                       \
    {                                                                                              \
        core::Variable<T> *var = currentIO.InquireVariable<T>(variableName);                       \
        if (var)                                                                                   \
        {                                                                                          \
            QueryVar *q = new QueryVar(variableName);                                              \
            adios2::Dims zero(var->Shape().size(), 0);                                             \
            adios2::Dims shape = var->Shape();                                                     \
            q->SetSelection(zero, shape);                                                          \
            return q;                                                                              \
        }                                                                                          \
    }
    ADIOS2_FOREACH_ATTRIBUTE_PRIMITIVE_STDTYPE_1ARG(declare_type)
#undef declare_type
    return nullptr;
}

void Worker::GetResultCoverage(std::vector<size_t> &touchedBlockIDs)
{
    touchedBlockIDs.clear();

    std::vector<BlockHit> blockHits;
    if (m_Query && m_SourceReader)
    {
        m_Query->BlockIndexEvaluate(m_SourceReader->m_IO, *m_SourceReader, blockHits);
    }

    for (auto blk : blockHits)
        touchedBlockIDs.push_back(blk.m_ID);
}

void Worker::GetResultCoverage(const adios2::Box<adios2::Dims> &outputRegion,
                               std::vector<Box<Dims>> &touchedBlocks)
{
    touchedBlocks.clear();

    if (!m_Query->UseOutputRegion(outputRegion))
    {
        helper::Throw<std::invalid_argument>("Toolkit", "query::Worker", "GetResultCoverage",
                                             "Unable to use the output region");
    }

    if (m_Query && m_SourceReader)
    {
        std::vector<BlockHit> blockHits;
        m_Query->BlockIndexEvaluate(m_SourceReader->m_IO, *m_SourceReader, blockHits);

        for (auto blk : blockHits)
            touchedBlocks.insert(touchedBlocks.end(), blk.m_Regions.begin(), blk.m_Regions.end());
    }
}
} // namespace query
} // namespace adios2

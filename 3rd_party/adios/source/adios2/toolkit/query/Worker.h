#ifndef QUERY_WORKER_H
#define QUERY_WORKER_H

#include <ios>      //std::ios_base::failure
#include <iostream> //std::cout

#include <stdexcept>
#include <string>
#include <vector>

#include "adios2/common/ADIOSTypes.h"
#include <fstream>

#include "Index.h"
#include "Query.h"
#include "Util.h"

// forward declaring to make pugixml a private dependency
namespace pugi
{
class xml_node;
}

namespace adios2
{
namespace query
{
class Worker
{
public:
    Worker(const Worker &other) = delete;

    Worker(Worker &&other)
    {
        this->m_QueryFile = other.m_QueryFile;
        this->m_SourceReader = other.m_SourceReader;
        this->m_Query = other.m_Query;
        other.m_Query = nullptr;
    }

    virtual ~Worker();

    adios2::core::Engine *GetSourceReader() { return m_SourceReader; }

    void GetResultCoverage(std::vector<size_t> &);
    void GetResultCoverage(const adios2::Box<adios2::Dims> &, std::vector<Box<adios2::Dims>> &);

protected:
    Worker(const std::string &configFile, adios2::core::Engine *adiosEngine);

    QueryVar *GetBasicVarQuery(adios2::core::IO &currentIO, const std::string &variableName);

    std::string m_QueryFile; // e.g. xml file

    adios2::core::Engine *m_SourceReader = nullptr;
    adios2::query::QueryBase *m_Query = nullptr;

private:
}; // worker

#ifdef ADIOS2_HAVE_DATAMAN
class JsonWorker : public Worker
{
public:
    JsonWorker(const std::string &configFile, adios2::core::Engine *adiosEngine)
    : Worker(configFile, adiosEngine)
    {
        ParseJson();
    }

private:
    void ParseJson();
};
#endif

class XmlWorker : public Worker
{
public:
    XmlWorker(const std::string &configFile, adios2::core::Engine *adiosEngine)
    : Worker(configFile, adiosEngine)
    {
        ParseMe();
    }

    void ParseMe();
    void Close() { std::cout << " .. closing in xml " << std::endl; }

private:
    void ConstructTree(RangeTree &host, const pugi::xml_node &node);
    void ConstructQuery(QueryVar &q, const pugi::xml_node &node);

    void ParseIONode(const pugi::xml_node &ioNode);
    adios2::query::QueryVar *ParseVarNode(const pugi::xml_node &node, adios2::core::IO &currentIO,
                                          adios2::core::Engine &reader);

}; // XmlWorker

Worker *GetWorker(const std::string &configFile, adios2::core::Engine *adiosEngine);

} // namespace query
} // name space adios2

#endif // QUERY_WORKER_H

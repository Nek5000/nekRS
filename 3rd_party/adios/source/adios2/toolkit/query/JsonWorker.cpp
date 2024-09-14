#include "Worker.h"
#include "adios2/helper/adiosLog.h"

#include <nlohmann_json.hpp>

namespace adios2
{
namespace query
{
namespace JsonUtil
{
bool HasEntry(nlohmann::json &jsonO, const char *name)
{
    int countMe = jsonO.count(name);
    if (countMe == 0)
        return false;
    return true;
}

void ConstructTree(adios2::query::RangeTree &host, nlohmann::json &opO)
{
    if (!HasEntry(opO, "value"))
        return;
    auto relationStr = opO["value"];
    host.SetRelation(adios2::query::strToRelation(relationStr));

    if (HasEntry(opO, "range"))
    {
        auto const rangeOs = opO["range"];
        for (auto r : rangeOs)
        {
            std::string valStr = r["value"];
            std::string opStr = r["compare"];
            host.AddLeaf(adios2::query::strToQueryOp(opStr), valStr);
        }
    }
    if (HasEntry(opO, "op"))
    {
        auto subOpO = opO["op"];
        if (subOpO.is_array())
        {
            for (auto sub : subOpO)
            {
                adios2::query::RangeTree subNode;
                ConstructTree(subNode, sub);
                host.AddNode(subNode);
            }
        }
        else
        {
            adios2::query::RangeTree subNode;
            ConstructTree(subNode, subOpO);
            host.AddNode(subNode);
        }
    }
} // construct tree

void LoadVarQuery(QueryVar *q, nlohmann::json &varO)
{
    if (!adios2::query::JsonUtil::HasEntry(varO, "op"))
        helper::Throw<std::ios_base::failure>("Toolkit", "query::JsonWorker", "LoadVarQuery",
                                              "No op entry specified for var:" + q->m_VarName);

    if (adios2::query::JsonUtil::HasEntry(varO, "boundingbox"))
    {
        auto bbO = varO["boundingbox"];
        q->LoadSelection(bbO["start"], bbO["count"]);
    }
    if (adios2::query::JsonUtil::HasEntry(varO, "op"))
    {
        auto opO = varO["op"];
        adios2::query::JsonUtil::ConstructTree(q->m_RangeTree, opO);
    }
} // LoadVarQuery
}
}
}

namespace adios2
{
namespace query
{
void JsonWorker::ParseJson()
{
    // local functions:

    auto lf_assertArray = [&](nlohmann::json &jsonO, const std::string &name) -> void {
        if (!jsonO.is_array())
            helper::Throw<std::ios_base::failure>("Toolkit", "query::JsonWorker", "ParseJson",
                                                  "Expecting Array for node:" + name);
    }; // lf assert

    auto lf_parseVar = [&](nlohmann::json &varO) -> QueryVar * {
        if (!adios2::query::JsonUtil::HasEntry(varO, "name"))
            helper::Throw<std::ios_base::failure>("Toolkit", "query::JsonWorker", "ParseJson",
                                                  "No var name specified!!");
        auto varName = (varO)["name"];
        adios2::core::IO &currIO = m_SourceReader->m_IO;
        const DataType varType = currIO.InquireVariableType(varName);
        if (varType == DataType::None)
        {
            helper::Log("Query", "JsonWorker", "ParseJson", "No such variable: " + varName.dump(),
                        helper::FATALERROR);
            return nullptr;
        }

        QueryVar *simpleQ = GetBasicVarQuery(currIO, varName);
        if (simpleQ)
            adios2::query::JsonUtil::LoadVarQuery(simpleQ, varO);
        return simpleQ;
    }; // local function to Parse var

    auto lf_parseTag = [&](nlohmann::json &tagO) -> QueryBase * {
        if (adios2::query::JsonUtil::HasEntry(tagO, "var"))
            return lf_parseVar(tagO["var"]);
        return nullptr;
    };

    std::ifstream fileStream(m_QueryFile);
    nlohmann::json jsonObj = nlohmann::json::parse(fileStream);

    if (!adios2::query::JsonUtil::HasEntry(jsonObj, "io"))
    {
        helper::Throw<std::ios_base::failure>(
            "Toolkit", "query::JsonWorker", "ParseJson",
            "No io node in json query file. Expecting the io node: " + m_SourceReader->m_IO.m_Name);
    }

    auto ioO = jsonObj.find("io");
    std::string const ioName = (*ioO)["name"];
    if (m_SourceReader->m_IO.m_Name.compare(ioName) != 0)
        helper::Throw<std::ios_base::failure>("Toolkit", "query::JsonWorker", "ParseJson",
                                              "invalid query io. Expecting io name = " +
                                                  m_SourceReader->m_IO.m_Name);
    if (adios2::query::JsonUtil::HasEntry(*ioO, "var"))
    {
        auto varO = (*ioO).find("var");
        m_Query = lf_parseVar(*varO);
        m_Query->Print();
        return;
    }
    if (!adios2::query::JsonUtil::HasEntry(*ioO, "query"))
        helper::Throw<std::ios_base::failure>("Toolkit", "query::JsonWorker", "ParseJson",
                                              "no query entry was defined for composite query");
    auto queryO = (*ioO)["query"];
    auto relationO = queryO["op"];
    QueryComposite *result = new QueryComposite(adios2::query::strToRelation(relationO));

    auto tagO = (*ioO)["tag"];
    lf_assertArray(tagO, "tag");
    std::map<std::string, QueryBase *> subqueries;

    for (auto tag : tagO)
    {
        QueryBase *q = lf_parseTag(tag);
        subqueries[tag["name"]] = q;
    }
    // for (nlohmann::json::iterator it = queryO.begin(); it != queryO.end();
    // ++it) {

    auto compO = queryO["comp"];
    lf_assertArray(compO, "comp");

    for (auto qname : compO)
    {
        std::cout << qname << std::endl;
        result->AddNode(subqueries[qname]);
    }

    m_Query = result;
    return;
} // parse
}
}

#include "Worker.h"

#include "adios2/helper/adiosLog.h"
#include "adios2/helper/adiosXMLUtil.h"

#include <pugixml.hpp>

namespace adios2
{
namespace query
{
void XmlWorker::ParseMe()
{
    auto lf_FileContents = [&](const std::string &configXML) -> std::string {
        std::ifstream fileStream(configXML);
        if (!fileStream)
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "query::XmlWorker", "ParseMe",
                                                  "file " + configXML + " not found");
        }
        std::ostringstream fileSS;
        fileSS << fileStream.rdbuf();
        fileStream.close();

        if (fileSS.str().empty())
        {
            helper::Throw<std::invalid_argument>("Toolkit", "query::XmlWorker", "ParseMe",
                                                 "config xml file is empty");
        }

        return fileSS.str();
    }; // local  function lf_FileContents

    const std::string fileContents = lf_FileContents(m_QueryFile);
    const std::unique_ptr<pugi::xml_document> document =
        adios2::helper::XMLDocument(fileContents, "in Query XMLWorker");

    const std::unique_ptr<pugi::xml_node> config =
        adios2::helper::XMLNode("adios-query", *document, "in adios-query", true);

    const pugi::xml_node ioNode = config->child("io");
    ParseIONode(ioNode);

} // Parse()

void XmlWorker::ParseIONode(const pugi::xml_node &ioNode)
{
#ifdef PARSE_IO
    const std::unique_ptr<pugi::xml_attribute> ioName =
        adios2::helper::XMLAttribute("name", ioNode, "in query");
    const std::unique_ptr<pugi::xml_attribute> fileName =
        adios2::helper::XMLAttribute("file", ioNode, "in query");

    // must be unique per io
    const std::unique_ptr<pugi::xml_node> &engineNode =
        adios2::helper::XMLNode("engine", ioNode, "in query", false, true);
    m_IO = &(m_Adios2.DeclareIO(ioName->value()));

    if (engineNode)
    {
        const std::unique_ptr<pugi::xml_attribute> type =
            adios2::query::XmlUtil::XMLAttribute("type", engineNode, "in query");
        m_IO->SetEngine(type->value());

        const adios2::Params parameters = helper::GetParameters(engineNode, "in query");
        m_IO->SetParameters(parameters);
    }
    else
    {
        m_IO->SetEngine("BPFile");
    }
    // adios2::Engine reader =  currIO.Open(fileName.value(),
    // adios2::Mode::Read, m_Comm);
    m_SourceReader = &(m_IO->Open(fileName.value(), adios2::Mode::Read, m_Comm));
#else
    const std::unique_ptr<pugi::xml_attribute> ioName =
        adios2::helper::XMLAttribute("name", ioNode, "in query");
    if (m_SourceReader->m_IO.m_Name.compare(ioName->value()) != 0)
    {
        helper::Throw<std::ios_base::failure>(
            "Toolkit", "query::XmlWorker", "ParseIONode",
            "invalid query io. Expecting io name = " + m_SourceReader->m_IO.m_Name +
                " found:" + ioName->value());
    }
#endif

    std::map<std::string, QueryBase *> subqueries;

    adios2::Box<adios2::Dims> ref;
    for (const pugi::xml_node &qTagNode : ioNode.children("tag"))
    {
        const std::unique_ptr<pugi::xml_attribute> name =
            adios2::helper::XMLAttribute("name", qTagNode, "in query");
        const pugi::xml_node &variable = qTagNode.child("var");
        QueryVar *q = ParseVarNode(variable, m_SourceReader->m_IO, *m_SourceReader);
        if (!q)
            continue;

        if (ref.first.size() == 0)
        {
            ref = q->m_Selection;
        }
        else if (!q->IsCompatible(ref))
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "query::XmlWorker", "ParseIONode",
                                                  "impactible query found on var:" +
                                                      q->GetVarName());
        }
        subqueries[name->value()] = q;
    }

    const pugi::xml_node &qNode = ioNode.child("query");
    if (qNode == nullptr)
    {
        const pugi::xml_node &variable = ioNode.child("var");
        m_Query = ParseVarNode(variable, m_SourceReader->m_IO, *m_SourceReader);
    }
    else
    {
        const std::unique_ptr<pugi::xml_attribute> op =
            adios2::helper::XMLAttribute("op", qNode, "in query");
        QueryComposite *q = new QueryComposite(adios2::query::strToRelation(op->value()));
        for (const pugi::xml_node &sub : qNode.children())
        {
            q->AddNode(subqueries[sub.name()]);
        }
        m_Query = q;
    }
} // parse_io_node

// node is the variable node
QueryVar *XmlWorker::ParseVarNode(const pugi::xml_node &node, adios2::core::IO &currentIO,
                                  adios2::core::Engine &reader)

{
    const std::string variableName =
        std::string(adios2::helper::XMLAttribute("name", node, "in query")->value());

    // const std::string varType = currentIO.VariableType(variableName);
    const DataType varType = currentIO.InquireVariableType(variableName);
    if (varType == DataType::None)
    {
        helper::Log("Query", "XmlWorker", "ParseVarNode", "No such variable: " + variableName,
                    helper::FATALERROR);
        helper::Throw<std::ios_base::failure>("Toolkit", "query::XmlWorker", "ParseVarNode",
                                              "variable: " + variableName + " not found");
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
            ConstructQuery(*q, node);                                                              \
            return q;                                                                              \
        }                                                                                          \
    }
    ADIOS2_FOREACH_ATTRIBUTE_PRIMITIVE_STDTYPE_1ARG(declare_type)
#undef declare_type
    return nullptr;
} //  parse_var_node

void XmlWorker::ConstructTree(RangeTree &host, const pugi::xml_node &node)
{
    std::string relationStr = adios2::helper::XMLAttribute("value", node, "in query")->value();
    host.SetRelation(adios2::query::strToRelation(relationStr));
    for (const pugi::xml_node &rangeNode : node.children("range"))
    {
        std::string valStr = adios2::helper::XMLAttribute("value", rangeNode, "in query")->value();
        std::string opStr = adios2::helper::XMLAttribute("compare", rangeNode, "in query")->value();

        host.AddLeaf(adios2::query::strToQueryOp(opStr), valStr);
    }

    for (const pugi::xml_node &opNode : node.children("op"))
    {
        adios2::query::RangeTree subNode;
        ConstructTree(subNode, opNode);
        host.AddNode(subNode);
    }
}

void XmlWorker::ConstructQuery(QueryVar &simpleQ, const pugi::xml_node &node)
{
    // QueryVar* simpleQ = new QueryVar(variableName);
    pugi::xml_node bbNode = node.child("boundingbox");

    if (bbNode)
    {
        std::string startStr = adios2::helper::XMLAttribute("start", bbNode, "in query")->value();
        std::string countStr = adios2::helper::XMLAttribute("count", bbNode, "in query")->value();

        adios2::Dims start = split(startStr, ',');
        adios2::Dims count = split(countStr, ',');

        if (start.size() != count.size())
        {
            helper::Throw<std::ios_base::failure>(
                "Toolkit", "query::XmlWorker", "ConstructQuery",
                "dim of startcount does match in the bounding box definition");
        }

        // simpleQ.setSelection(box.first, box.second);
        adios2::Dims shape = simpleQ.m_Selection.second; // set at the creation for default
        simpleQ.SetSelection(start, count);
        if (!simpleQ.IsSelectionValid(shape))
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "query::XmlWorker", "ConstructQuery",
                                                  "invalid selections for selection of var: " +
                                                      simpleQ.GetVarName());
        }
    }

#ifdef NEVER // don't know whether this is useful.
    pugi::xml_node tsNode = node.child("tstep");
    if (tsNode)
    {
        std::string startStr = adios2::helper::XMLAttribute("start", tsNode, "in query").value();
        std::string countStr = adios2::helper::XMLAttribute("count", tsNode, "in query").value();

        if ((startStr.size() > 0) && (countStr.size() > 0))
        {
            std::stringstream ss(startStr), cc(countStr);
            ss >> simpleQ.m_TimestepStart;
            cc >> simpleQ.m_TimestepCount;
        }
    }
#endif
    pugi::xml_node relationNode = node.child("op");
    ConstructTree(simpleQ.m_RangeTree, relationNode);
}

} // namespace query
} // namespace adios2

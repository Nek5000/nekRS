#include "VariableDerived.h"
#include "adios2/helper/adiosType.h"

namespace adios2
{
namespace core
{

VariableDerived::VariableDerived(const std::string &name, adios2::derived::Expression expr,
                                 const DataType exprType, const bool isConstant,
                                 const DerivedVarType varType)
: VariableBase(name, exprType, helper::GetDataTypeSize(exprType), expr.GetShape(), expr.GetStart(),
               expr.GetCount(), isConstant),
  m_DerivedType(varType), m_Expr(expr)
{
}

DerivedVarType VariableDerived::GetDerivedType() { return m_DerivedType; }

std::vector<std::string> VariableDerived::VariableNameList() { return m_Expr.VariableNameList(); }
void VariableDerived::UpdateExprDim(std::map<std::string, std::tuple<Dims, Dims, Dims>> NameToDims)
{
    m_Expr.SetDims(NameToDims);
    m_Shape = m_Expr.GetShape();
    m_Start = m_Expr.GetStart();
    m_Count = m_Expr.GetCount();
}

std::vector<std::tuple<void *, Dims, Dims>>
VariableDerived::ApplyExpression(std::map<std::string, std::unique_ptr<MinVarInfo>> &NameToMVI)
{
    size_t numBlocks = 0;
    // check that all variables have the same number of blocks
    for (const auto &variable : NameToMVI)
    {
        if (numBlocks == 0)
            numBlocks = variable.second->BlocksInfo.size();
        if (numBlocks != variable.second->BlocksInfo.size())
            helper::Throw<std::invalid_argument>("Core", "VariableDerived", "ApplyExpression",
                                                 " variables do not have the same number of blocks "
                                                 " in computing the derived variable " +
                                                     m_Name);
    }

    std::map<std::string, std::vector<adios2::derived::DerivedData>> inputData;
    // create the map between variable name and DerivedData object
    for (const auto &variable : NameToMVI)
    {
        // add the dimensions of all blocks into a vector
        std::vector<adios2::derived::DerivedData> varData;
        for (size_t i = 0; i < numBlocks; i++)
        {
            Dims start;
            Dims count;
            for (int d = 0; d < variable.second->Dims; d++)
            {
                start.push_back(variable.second->BlocksInfo[i].Start[d]);
                count.push_back(variable.second->BlocksInfo[i].Count[d]);
            }
            varData.push_back(adios2::derived::DerivedData(
                {variable.second->BlocksInfo[i].BufferP, start, count}));
        }
        inputData.insert({variable.first, varData});
    }
    // TODO check that the dimensions are still corrects
    std::vector<adios2::derived::DerivedData> outputData =
        m_Expr.ApplyExpression(m_Type, numBlocks, inputData);

    std::vector<std::tuple<void *, Dims, Dims>> blockData;
    for (size_t i = 0; i < numBlocks; i++)
    {
        blockData.push_back({outputData[i].Data, outputData[i].Start, outputData[i].Count});
    }

    return blockData;
}

std::vector<void *>
VariableDerived::ApplyExpression(std::map<std::string, std::vector<void *>> NameToData,
                                 std::map<std::string, std::tuple<Dims, Dims, Dims>> NameToDims)
{
    size_t numBlocks = 0;
    std::map<std::string, std::vector<adios2::derived::DerivedData>> inputData;
    // check that all variables have the same number of blocks
    for (auto variable : NameToData)
    {
        if (numBlocks == 0)
            numBlocks = variable.second.size();
        if (numBlocks != variable.second.size())
            helper::Throw<std::invalid_argument>("Core", "VariableDerived", "ApplyExpression",
                                                 " variables do not have the same number of blocks "
                                                 " in computing the derived variable " +
                                                     m_Name);
    }
    std::cout << "Derived variable " << m_Name
              << ": PASS : variables have written the same num of blocks" << std::endl;
    // create the map between variable name and DerivedData object
    for (auto variable : NameToData)
    {
        // add the dimensions of all blocks into a vector
        std::vector<adios2::derived::DerivedData> varData;
        for (size_t i = 0; i < numBlocks; i++)
        {
            varData.push_back(adios2::derived::DerivedData(
                {variable.second[i], std::get<0>(NameToDims[variable.first]),
                 std::get<1>(NameToDims[variable.first])}));
        }
        inputData.insert({variable.first, varData});
    }
    std::vector<adios2::derived::DerivedData> outputData =
        m_Expr.ApplyExpression(m_Type, numBlocks, inputData);
    std::vector<void *> blockData;
    for (size_t i = 0; i < numBlocks; i++)
    {
        blockData.push_back(outputData[i].Data);
    }
    return blockData;
}

} // end namespace core
} // end namespace adios2
